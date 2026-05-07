// ---------------------------------------------------------------------
//
// Copyright (C) 2026 by Luca Heltai
//
// This file is part of the deal.II Testbench for NMOPT Laboratories,
// deal.II library.
//
// The deal.II Testbench for NMOPT Laboratories is free software; you can use
// it, redistribute it, and/or modify it under the terms of the Apache-2.0
// License WITH LLVM-exception as published by the Free Software Foundation;
// either version 3.0 of the License, or (at your option) any later version. The
// full text of the license can be found in the file LICENSE.md at the top level
// of the deal.II Testbench for NMOPT Laboratories distribution.
//
// ---------------------------------------------------------------------

#include "parabolic_projected_gradient.h"

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace dealii;


template <int dim>
ParabolicProjectedGradient<dim>::ParabolicProjectedGradient(
  const std::string &subsection_name)
  : ParameterAcceptor(subsection_name)
  , fe(std::make_unique<FESystem<dim>>(FE_Q<dim>(1),
                                       1,
                                       FE_Q<dim>(1),
                                       1,
                                       FE_DGQ<dim>(0),
                                       1))
  , dof_handler(triangulation)
  , dofs_per_block(n_blocks)
  , state_degree(1)
  , control_degree(0)
  , continuous_control(false)
  , global_refinements(5)
  , final_time(1.0)
  , n_time_steps(20)
  , time_step(final_time / n_time_steps)
  , regularization(1e-3)
  , terminal_tracking_weight(1.0)
  , max_iterations(20)
  , gradient_tolerance(1e-8)
  , projected_gradient_tolerance(1e-8)
  , fixed_step_length(1.0)
  , output_name("parabolic_projected_gradient")
  , grid_generator_function("hyper_cube")
  , grid_generator_arguments("0.0 : 1.0 : false")
  , diffusion_coefficient(1)
  , reaction_coefficient(1)
  , forcing_term(1)
  , initial_state(1)
  , initial_control(1)
  , desired_state(1)
  , terminal_target_function(1)
  , lower_bound_function(1)
  , upper_bound_function(1)
  , state_boundary_values(1)
  , adjoint_boundary_values(1)
{}



template <int dim>
std::string
ParabolicProjectedGradient<dim>::default_variable_names()
{
  const std::vector<std::string> spatial_names = {"x", "y", "z"};

  std::string names = spatial_names[0];
  for (unsigned int d = 1; d < dim; ++d)
    names += "," + spatial_names[d];
  names += ",t";

  return names;
}



template <int dim>
void
ParabolicProjectedGradient<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.add_parameter("State degree", state_degree);
  prm.add_parameter("Control degree", control_degree);
  prm.add_parameter("Continuous control", continuous_control);
  prm.add_parameter("Global refinements", global_refinements);
  prm.add_parameter("Final time", final_time);
  prm.add_parameter("Number of time steps", n_time_steps);
  prm.add_parameter("Regularization", regularization);
  prm.add_parameter("Terminal tracking weight", terminal_tracking_weight);
  prm.add_parameter("Max optimization iterations", max_iterations);
  prm.add_parameter("Gradient tolerance", gradient_tolerance);
  prm.add_parameter("Projected gradient tolerance",
                    projected_gradient_tolerance);
  prm.add_parameter("Fixed step length", fixed_step_length);
  prm.add_parameter("Output name", output_name);
  prm.add_parameter("Grid generator function", grid_generator_function);
  prm.add_parameter("Grid generator arguments", grid_generator_arguments);

  const std::string variables = default_variable_names();

  prm.enter_subsection("Armijo exercise parameters");
  armijo_parameters.add_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Diffusion coefficient");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "1.0");
  prm.leave_subsection();

  prm.enter_subsection("Reaction coefficient");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "0.0");
  prm.leave_subsection();

  prm.enter_subsection("Forcing term");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "0.0");
  prm.leave_subsection();

  prm.enter_subsection("Initial state");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "0.0");
  prm.leave_subsection();

  prm.enter_subsection("Initial control");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "0.0");
  prm.leave_subsection();

  prm.enter_subsection("Desired state");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set(
    "Function expression",
    "x > 0.35 ? (x < 0.65 ? (t > 0.25 ? (t < 0.75 ? 1 : 0) : 0) : 0) : 0");
  prm.leave_subsection();

  prm.enter_subsection("Terminal target");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "0.0");
  prm.leave_subsection();

  prm.enter_subsection("Lower bound");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "-2.0");
  prm.leave_subsection();

  prm.enter_subsection("Upper bound");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "2.0");
  prm.leave_subsection();

  prm.enter_subsection("State boundary values");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "0.0");
  prm.leave_subsection();

  prm.enter_subsection("Adjoint boundary values");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Variable names", variables);
  prm.set("Function expression", "0.0");
  prm.leave_subsection();
}



template <int dim>
void
ParabolicProjectedGradient<dim>::parse_parameters(ParameterHandler &prm)
{
  AssertThrow(n_time_steps > 0,
              ExcMessage("The number of time steps must be positive."));
  AssertThrow(final_time > 0.0, ExcMessage("Final time must be positive."));
  AssertThrow(regularization > 0.0,
              ExcMessage("Regularization must be positive."));
  AssertThrow(fixed_step_length > 0.0,
              ExcMessage("Fixed step length must be positive."));

  time_step = final_time / static_cast<double>(n_time_steps);
  fe        = create_fe_system();

  prm.enter_subsection("Diffusion coefficient");
  diffusion_coefficient.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Reaction coefficient");
  reaction_coefficient.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Forcing term");
  forcing_term.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Initial state");
  initial_state.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Initial control");
  initial_control.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Desired state");
  desired_state.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Terminal target");
  terminal_target_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Lower bound");
  lower_bound_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Upper bound");
  upper_bound_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("State boundary values");
  state_boundary_values.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Adjoint boundary values");
  adjoint_boundary_values.parse_parameters(prm);
  prm.leave_subsection();
}



template <int dim>
std::unique_ptr<FESystem<dim>>
ParabolicProjectedGradient<dim>::create_fe_system() const
{
  if (continuous_control)
    return std::make_unique<FESystem<dim>>(FE_Q<dim>(state_degree),
                                           1,
                                           FE_Q<dim>(state_degree),
                                           1,
                                           FE_Q<dim>(control_degree),
                                           1);

  return std::make_unique<FESystem<dim>>(FE_Q<dim>(state_degree),
                                         1,
                                         FE_Q<dim>(state_degree),
                                         1,
                                         FE_DGQ<dim>(control_degree),
                                         1);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::make_grid()
{
  triangulation.clear();
  GridGenerator::generate_from_name_and_arguments(triangulation,
                                                  grid_generator_function,
                                                  grid_generator_arguments);

  for (const auto &cell : triangulation.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        face->set_boundary_id(0);

  triangulation.refine_global(global_refinements);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::initialize()
{
  make_grid();
  setup_system();
  assemble_matrices();
  initialize_time_step_solvers();
  initialize_time_data();
}



template <int dim>
void
ParabolicProjectedGradient<dim>::setup_system()
{
  dof_handler.distribute_dofs(*fe);

  const std::vector<unsigned int> block_component = {state_block,
                                                     adjoint_block,
                                                     control_block};
  DoFRenumbering::component_wise(dof_handler, block_component);

  dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  state_boundary_values.set_time(0.0);
  adjoint_boundary_values.set_time(0.0);

  const VectorFunctionFromScalarFunctionObject<dim> state_bc(
    [this](const Point<dim> &p) { return state_boundary_values.value(p); },
    state_block,
    n_blocks);
  const VectorFunctionFromScalarFunctionObject<dim> adjoint_bc(
    [this](const Point<dim> &p) { return adjoint_boundary_values.value(p); },
    adjoint_block,
    n_blocks);

  VectorTools::interpolate_boundary_values(
    dof_handler,
    0,
    state_bc,
    constraints,
    fe->component_mask(FEValuesExtractors::Scalar(state_block)));
  VectorTools::interpolate_boundary_values(
    dof_handler,
    0,
    adjoint_bc,
    constraints,
    fe->component_mask(FEValuesExtractors::Scalar(adjoint_block)));
  constraints.close();

  BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs = */ false);
  sparsity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparsity_pattern);
  stiffness_matrix.reinit(sparsity_pattern);
  coupling_matrix.reinit(sparsity_pattern);
  tracking_matrix.reinit(sparsity_pattern);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::assemble_matrices()
{
  mass_matrix      = 0;
  stiffness_matrix = 0;
  coupling_matrix  = 0;
  tracking_matrix  = 0;

  const QGauss<dim> quadrature_formula(std::max(state_degree, control_degree) +
                                       1);
  FEValues<dim>     fe_values(*fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Scalar state(state_block);
  const FEValuesExtractors::Scalar adjoint(adjoint_block);
  const FEValuesExtractors::Scalar control(control_block);

  const unsigned int dofs_per_cell = fe->n_dofs_per_cell();

  FullMatrix<double> cell_mass(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_stiffness(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_coupling(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_tracking(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_mass      = 0;
      cell_stiffness = 0;
      cell_coupling  = 0;
      cell_tracking  = 0;

      fe_values.reinit(cell);

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          const Point<dim> &x_q       = fe_values.quadrature_point(q_index);
          const double      diffusion = diffusion_coefficient.value(x_q);
          const double      reaction  = reaction_coefficient.value(x_q);

          for (const unsigned int i : fe_values.dof_indices())
            {
              const double phi_state_i = fe_values[state].value(i, q_index);
              const Tensor<1, dim> grad_state_i =
                fe_values[state].gradient(i, q_index);
              const double phi_adjoint_i = fe_values[adjoint].value(i, q_index);
              const Tensor<1, dim> grad_adjoint_i =
                fe_values[adjoint].gradient(i, q_index);
              const double phi_control_i = fe_values[control].value(i, q_index);

              for (const unsigned int j : fe_values.dof_indices())
                {
                  const double phi_state_j = fe_values[state].value(j, q_index);
                  const Tensor<1, dim> grad_state_j =
                    fe_values[state].gradient(j, q_index);
                  const double phi_adjoint_j =
                    fe_values[adjoint].value(j, q_index);
                  const Tensor<1, dim> grad_adjoint_j =
                    fe_values[adjoint].gradient(j, q_index);
                  const double phi_control_j =
                    fe_values[control].value(j, q_index);

                  const double dx = fe_values.JxW(q_index);

                  cell_mass(i, j) +=
                    (phi_state_i * phi_state_j + phi_adjoint_i * phi_adjoint_j +
                     phi_control_i * phi_control_j) *
                    dx;

                  cell_stiffness(i, j) +=
                    (diffusion * grad_state_i * grad_state_j +
                     reaction * phi_state_i * phi_state_j +
                     diffusion * grad_adjoint_i * grad_adjoint_j +
                     reaction * phi_adjoint_i * phi_adjoint_j) *
                    dx;

                  cell_coupling(i, j) += (phi_state_i * phi_control_j +
                                          phi_adjoint_i * phi_control_j) *
                                         dx;

                  cell_tracking(i, j) += phi_adjoint_i * phi_state_j * dx;
                }
            }
        }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_mass,
                                             local_dof_indices,
                                             mass_matrix);
      constraints.distribute_local_to_global(cell_stiffness,
                                             local_dof_indices,
                                             stiffness_matrix);
      constraints.distribute_local_to_global(cell_coupling,
                                             local_dof_indices,
                                             coupling_matrix);
      constraints.distribute_local_to_global(cell_tracking,
                                             local_dof_indices,
                                             tracking_matrix);
    }
}



template <int dim>
void
ParabolicProjectedGradient<dim>::initialize_time_step_solvers()
{
  state_time_step_matrix.reinit(
    sparsity_pattern.block(state_block, state_block));
  state_time_step_matrix.copy_from(mass_matrix.block(state_block, state_block));
  state_time_step_matrix.add(time_step,
                             stiffness_matrix.block(state_block, state_block));

  adjoint_time_step_matrix.reinit(
    sparsity_pattern.block(adjoint_block, adjoint_block));
  adjoint_time_step_matrix.copy_from(
    mass_matrix.block(adjoint_block, adjoint_block));
  adjoint_time_step_matrix.add(time_step,
                               stiffness_matrix.block(adjoint_block,
                                                      adjoint_block));

  state_time_step_inverse.initialize(state_time_step_matrix);
  adjoint_time_step_inverse.initialize(adjoint_time_step_matrix);
  control_mass_inverse.initialize(
    mass_matrix.block(control_block, control_block));
}



template <int dim>
void
ParabolicProjectedGradient<dim>::reinit_time_vector(
  TimeVector        &values,
  const unsigned int block,
  const bool         include_terminal_slot) const
{
  values.resize(n_time_steps + (include_terminal_slot ? 2 : 1));
  for (auto &value : values)
    {
      value.reinit(dofs_per_block[block]);
      value = 0;
    }
}



template <int dim>
void
ParabolicProjectedGradient<dim>::interpolate_function_to_block(
  Functions::ParsedFunction<dim> &function,
  const unsigned int              block,
  const double                    time,
  Vector<double>                 &vector) const
{
  BlockVector<double> tmp(dofs_per_block);
  tmp = 0;

  function.set_time(time);
  const VectorFunctionFromScalarFunctionObject<dim> vector_function(
    [&function](const Point<dim> &point) { return function.value(point); },
    block,
    n_blocks);

  VectorTools::interpolate(dof_handler, vector_function, tmp);
  constraints.distribute(tmp);
  vector = tmp.block(block);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::distribute_block(const unsigned int block,
                                                  Vector<double> &vector) const
{
  BlockVector<double> tmp(dofs_per_block);
  tmp              = 0;
  tmp.block(block) = vector;
  constraints.distribute(tmp);
  vector = tmp.block(block);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::initialize_time_data()
{
  reinit_time_vector(lower_bound, control_block);
  reinit_time_vector(upper_bound, control_block);
  reinit_time_vector(target_state, state_block);

  terminal_target.reinit(dofs_per_block[state_block]);

  for (unsigned int k = 0; k <= n_time_steps; ++k)
    {
      const double t = k * time_step;

      interpolate_function_to_block(lower_bound_function,
                                    control_block,
                                    t,
                                    lower_bound[k]);
      interpolate_function_to_block(upper_bound_function,
                                    control_block,
                                    t,
                                    upper_bound[k]);
      interpolate_function_to_block(desired_state,
                                    state_block,
                                    t,
                                    target_state[k]);

      for (unsigned int i = 0; i < lower_bound[k].size(); ++i)
        AssertThrow(lower_bound[k][i] <= upper_bound[k][i],
                    ExcMessage("Lower control bound exceeds upper bound."));
    }

  interpolate_function_to_block(terminal_target_function,
                                state_block,
                                final_time,
                                terminal_target);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::initialize_control(TimeVector &control) const
{
  reinit_time_vector(control, control_block);

  for (unsigned int k = 0; k <= n_time_steps; ++k)
    const_cast<ParabolicProjectedGradient<dim> *>(this)
      ->interpolate_function_to_block(
        const_cast<Functions::ParsedFunction<dim> &>(initial_control),
        control_block,
        k * time_step,
        control[k]);

  project_onto_box(control);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::project_onto_box(TimeVector &control) const
{
  for (unsigned int k = 0; k <= n_time_steps; ++k)
    project_time_level(k, control[k]);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::project_time_level(
  const unsigned int time_index,
  Vector<double>    &control) const
{
  for (unsigned int i = 0; i < control.size(); ++i)
    control[i] = std::clamp(control[i],
                            lower_bound[time_index][i],
                            upper_bound[time_index][i]);

  distribute_block(control_block, control);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::solve_state(const TimeVector &control,
                                             TimeVector       &state) const
{
  reinit_time_vector(state, state_block);

  const_cast<ParabolicProjectedGradient<dim> *>(this)
    ->interpolate_function_to_block(
      const_cast<Functions::ParsedFunction<dim> &>(initial_state),
      state_block,
      0.0,
      state[0]);

  Vector<double> rhs(dofs_per_block[state_block]);
  Vector<double> tmp(dofs_per_block[state_block]);
  Vector<double> forcing(dofs_per_block[state_block]);
  Vector<double> forcing_coefficients(dofs_per_block[state_block]);

  for (unsigned int k = 1; k <= n_time_steps; ++k)
    {
      rhs = 0;

      mass_matrix.block(state_block, state_block).vmult(rhs, state[k - 1]);

      coupling_matrix.block(state_block, control_block).vmult(tmp, control[k]);
      rhs.add(time_step, tmp);

      const_cast<ParabolicProjectedGradient<dim> *>(this)
        ->interpolate_function_to_block(
          const_cast<Functions::ParsedFunction<dim> &>(forcing_term),
          state_block,
          k * time_step,
          forcing_coefficients);
      mass_matrix.block(state_block, state_block)
        .vmult(forcing, forcing_coefficients);
      rhs.add(time_step, forcing);

      state_time_step_inverse.vmult(state[k], rhs);
      distribute_block(state_block, state[k]);
    }
}



template <int dim>
void
ParabolicProjectedGradient<dim>::solve_adjoint(const TimeVector &state,
                                               TimeVector       &adjoint) const
{
  reinit_time_vector(adjoint, adjoint_block, true);

  Vector<double> rhs(dofs_per_block[adjoint_block]);
  Vector<double> tmp(dofs_per_block[adjoint_block]);
  Vector<double> mismatch(dofs_per_block[state_block]);

  for (unsigned int k = n_time_steps; k > 0; --k)
    {
      rhs = 0;

      if (k < n_time_steps)
        {
          mass_matrix.block(adjoint_block, adjoint_block)
            .vmult(rhs, adjoint[k + 1]);
        }

      mismatch = state[k];
      mismatch -= target_state[k];
      tracking_matrix.block(adjoint_block, state_block).vmult(tmp, mismatch);
      rhs.add(time_step, tmp);

      if (k == n_time_steps && terminal_tracking_weight != 0.0)
        {
          mismatch = state[k];
          mismatch -= terminal_target;
          tracking_matrix.block(adjoint_block, state_block)
            .vmult(tmp, mismatch);
          rhs.add(terminal_tracking_weight, tmp);
        }

      adjoint_time_step_inverse.vmult(adjoint[k], rhs);
      distribute_block(adjoint_block, adjoint[k]);
    }
}



template <int dim>
void
ParabolicProjectedGradient<dim>::compute_gradient(const TimeVector &control,
                                                  const TimeVector &adjoint,
                                                  TimeVector &gradient) const
{
  reinit_time_vector(gradient, control_block);

  Vector<double> dual_gradient(dofs_per_block[control_block]);
  Vector<double> tmp(dofs_per_block[control_block]);

  for (unsigned int k = 1; k <= n_time_steps; ++k)
    {
      mass_matrix.block(control_block, control_block)
        .vmult(dual_gradient, control[k]);
      dual_gradient *= regularization;

      coupling_matrix.block(adjoint_block, control_block)
        .Tvmult(tmp, adjoint[k]);
      dual_gradient += tmp;

      control_mass_inverse.vmult(gradient[k], dual_gradient);
      distribute_block(control_block, gradient[k]);
    }
}



template <int dim>
void
ParabolicProjectedGradient<dim>::projected_gradient_step(
  TimeVector       &trial_control,
  const TimeVector &gradient,
  const double      step_length) const
{
  for (unsigned int k = 1; k <= n_time_steps; ++k)
    {
      trial_control[k].add(-step_length, gradient[k]);
      project_time_level(k, trial_control[k]);
    }

  trial_control[0] = trial_control[1];
}



template <int dim>
double
ParabolicProjectedGradient<dim>::control_inner_product(
  const TimeVector &left,
  const TimeVector &right) const
{
  Vector<double> weighted(dofs_per_block[control_block]);
  double         result = 0.0;

  for (unsigned int k = 1; k <= n_time_steps; ++k)
    {
      mass_matrix.block(control_block, control_block).vmult(weighted, right[k]);
      result += time_step * (left[k] * weighted);
    }

  return result;
}



template <int dim>
double
ParabolicProjectedGradient<dim>::control_norm(const TimeVector &control) const
{
  return std::sqrt(std::max(0.0, control_inner_product(control, control)));
}



template <int dim>
typename ParabolicProjectedGradient<dim>::Diagnostics
ParabolicProjectedGradient<dim>::compute_diagnostics(
  const TimeVector &state,
  const TimeVector &control,
  const TimeVector &gradient,
  const double      step_length) const
{
  Diagnostics diagnostics;
  diagnostics.step_length = step_length;

  Vector<double> mismatch(dofs_per_block[state_block]);
  Vector<double> weighted_state(dofs_per_block[state_block]);
  Vector<double> weighted_control(dofs_per_block[control_block]);

  for (unsigned int k = 1; k <= n_time_steps; ++k)
    {
      mismatch = state[k];
      mismatch -= target_state[k];
      mass_matrix.block(state_block, state_block)
        .vmult(weighted_state, mismatch);
      diagnostics.tracking_cost +=
        0.5 * time_step * (mismatch * weighted_state);

      mass_matrix.block(control_block, control_block)
        .vmult(weighted_control, control[k]);
      diagnostics.control_cost +=
        0.5 * regularization * time_step * (control[k] * weighted_control);
    }

  mismatch = state[n_time_steps];
  mismatch -= terminal_target;
  mass_matrix.block(state_block, state_block).vmult(weighted_state, mismatch);
  diagnostics.terminal_tracking_cost =
    0.5 * terminal_tracking_weight * (mismatch * weighted_state);

  diagnostics.cost = diagnostics.tracking_cost +
                     diagnostics.terminal_tracking_cost +
                     diagnostics.control_cost;
  diagnostics.gradient_norm = control_norm(gradient);
  diagnostics.projected_gradient_norm =
    compute_projected_gradient_norm(control, gradient, step_length);
  diagnostics.control_norm = control_norm(control);

  return diagnostics;
}



template <int dim>
double
ParabolicProjectedGradient<dim>::compute_projected_gradient_norm(
  const TimeVector &control,
  const TimeVector &gradient,
  const double      step_length) const
{
  TimeVector trial_control = control;
  projected_gradient_step(trial_control, gradient, step_length);

  TimeVector residual = control;
  for (unsigned int k = 1; k <= n_time_steps; ++k)
    {
      residual[k] -= trial_control[k];
      residual[k] /= step_length;
    }

  return control_norm(residual);
}



template <int dim>
double
ParabolicProjectedGradient<dim>::choose_step_length(
  const TimeVector & /*control*/,
  const TimeVector & /*gradient*/,
  const Diagnostics & /*diagnostics*/) const
{
  // Laboratory exercise: replace this fixed value by a projected Armijo
  // backtracking loop.  The ingredients are already available in this class:
  // projected_gradient_step(), solve_state(), compute_diagnostics(), and
  // control_inner_product().
  return fixed_step_length;
}



template <int dim>
void
ParabolicProjectedGradient<dim>::output_time_sequence(
  const unsigned int iteration,
  const TimeVector  &state,
  const TimeVector  &adjoint,
  const TimeVector  &control) const
{
  std::vector<std::pair<double, std::string>> times_and_names;

  for (unsigned int k = 0; k <= n_time_steps; ++k)
    {
      BlockVector<double> solution(dofs_per_block);
      BlockVector<double> target(dofs_per_block);
      BlockVector<double> terminal_target_output(dofs_per_block);
      BlockVector<double> lower_bound_output(dofs_per_block);
      BlockVector<double> upper_bound_output(dofs_per_block);

      solution               = 0;
      target                 = 0;
      terminal_target_output = 0;
      lower_bound_output     = 0;
      upper_bound_output     = 0;

      solution.block(state_block) = state[k];
      if (k > 0)
        solution.block(adjoint_block) = adjoint[k];
      solution.block(control_block) = control[k];

      target.block(state_block)                 = target_state[k];
      terminal_target_output.block(state_block) = terminal_target;
      lower_bound_output.block(control_block)   = lower_bound[k];
      upper_bound_output.block(control_block)   = upper_bound[k];

      std::vector<std::string> solution_names = {"state", "adjoint", "control"};
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation(n_blocks,
                       DataComponentInterpretation::component_is_scalar);

      DataOut<dim> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(dof_handler,
                               solution,
                               solution_names,
                               interpretation);
      data_out.add_data_vector(dof_handler,
                               target,
                               {"desired_state",
                                "zero_adjoint_target",
                                "zero_control_target"},
                               interpretation);
      data_out.add_data_vector(dof_handler,
                               terminal_target_output,
                               {"terminal_target",
                                "zero_adjoint_terminal",
                                "zero_control_terminal"},
                               interpretation);
      data_out.add_data_vector(dof_handler,
                               lower_bound_output,
                               {"zero_state_lower",
                                "zero_adjoint_lower",
                                "lower_bound"},
                               interpretation);
      data_out.add_data_vector(dof_handler,
                               upper_bound_output,
                               {"zero_state_upper",
                                "zero_adjoint_upper",
                                "upper_bound"},
                               interpretation);
      data_out.build_patches();

      const std::string filename =
        output_name + "-opt" + Utilities::int_to_string(iteration, 4) +
        "-time" + Utilities::int_to_string(k, 4) + ".vtu";
      std::ofstream output(filename);
      data_out.write_vtu(output);

      times_and_names.emplace_back(k * time_step, filename);
    }

  const std::string pvd_filename =
    output_name + "-opt" + Utilities::int_to_string(iteration, 4) + ".pvd";
  std::ofstream pvd_output(pvd_filename);
  DataOutBase::write_pvd_record(pvd_output, times_and_names);
}



template <int dim>
void
ParabolicProjectedGradient<dim>::write_csv_header(std::ofstream &csv) const
{
  csv << "iteration,cost,tracking_cost,terminal_tracking_cost,control_cost,"
      << "gradient_norm,projected_gradient_norm,control_norm,step_length\n";
}



template <int dim>
void
ParabolicProjectedGradient<dim>::write_csv_line(
  std::ofstream     &csv,
  const unsigned int iteration,
  const Diagnostics &diagnostics) const
{
  csv << iteration << ',' << diagnostics.cost << ','
      << diagnostics.tracking_cost << ',' << diagnostics.terminal_tracking_cost
      << ',' << diagnostics.control_cost << ',' << diagnostics.gradient_norm
      << ',' << diagnostics.projected_gradient_norm << ','
      << diagnostics.control_norm << ',' << diagnostics.step_length << '\n';
}



template <int dim>
void
ParabolicProjectedGradient<dim>::run()
{
  initialize();

  std::cout << "Number of active cells:       "
            << triangulation.n_active_cells() << std::endl;
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  std::cout << "State dofs:                   " << dofs_per_block[state_block]
            << std::endl;
  std::cout << "Adjoint dofs:                 " << dofs_per_block[adjoint_block]
            << std::endl;
  std::cout << "Control dofs:                 " << dofs_per_block[control_block]
            << std::endl;
  std::cout << "Time step:                    " << time_step << std::endl;

  TimeVector control;
  TimeVector state;
  TimeVector adjoint;
  TimeVector gradient;

  initialize_control(control);

  std::ofstream csv(output_name + "-history.csv");
  write_csv_header(csv);

  for (unsigned int iteration = 0; iteration <= max_iterations; ++iteration)
    {
      solve_state(control, state);
      solve_adjoint(state, adjoint);
      compute_gradient(control, adjoint, gradient);

      Diagnostics diagnostics =
        compute_diagnostics(state, control, gradient, fixed_step_length);

      const double step_length =
        choose_step_length(control, gradient, diagnostics);
      diagnostics.step_length = step_length;
      diagnostics.projected_gradient_norm =
        compute_projected_gradient_norm(control, gradient, step_length);

      write_csv_line(csv, iteration, diagnostics);
      output_time_sequence(iteration, state, adjoint, control);

      std::cout << "Projected GD iteration " << iteration
                << ": J=" << diagnostics.cost
                << ", ||g||=" << diagnostics.gradient_norm
                << ", ||g_P||=" << diagnostics.projected_gradient_norm
                << ", step=" << diagnostics.step_length << std::endl;

      if (diagnostics.gradient_norm < gradient_tolerance ||
          diagnostics.projected_gradient_norm < projected_gradient_tolerance)
        break;

      TimeVector next_control = control;
      projected_gradient_step(next_control, gradient, step_length);
      control = next_control;
    }
}


template class ParabolicProjectedGradient<DEAL_DIMENSION>;

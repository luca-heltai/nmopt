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
// either version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at the top
// level of the deal.II Testbench for NMOPT Laboratories distribution.
//
// ---------------------------------------------------------------------

#include "inverse_poisson_kkt.h"

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>

using namespace dealii;


template <int dim>
InversePoissonKKT<dim>::InversePoissonKKT(const std::string &subsection_name)
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
  , global_refinements(5)
  , max_pdas_iterations(30)
  , regularization(1e-4)
  , active_set_parameter(1.0)
  , newton_tolerance(1e-10)
  , grid_generator_function("hyper_cube")
  , grid_generator_arguments("0.0 : 1.0 : false")
  , output_name("inverse_poisson_kkt")
  , forcing_term(1)
  , desired_state_function(1)
  , reference_coefficient_function(1)
  , initial_coefficient_function(1)
  , exact_coefficient_function(1)
  , lower_bound_function(1)
  , upper_bound_function(1)
{}



template <int dim>
void
InversePoissonKKT<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.add_parameter("State degree", state_degree);
  prm.add_parameter("Control degree", control_degree);
  prm.add_parameter("Global refinements", global_refinements);
  prm.add_parameter("Regularization", regularization);
  prm.add_parameter("Max PDAS iterations", max_pdas_iterations);
  prm.add_parameter("Active-set parameter", active_set_parameter);
  prm.add_parameter("Newton tolerance", newton_tolerance);
  prm.add_parameter("Grid generator function", grid_generator_function);
  prm.add_parameter("Grid generator arguments", grid_generator_arguments);
  prm.add_parameter("Output name", output_name);

  prm.enter_subsection("Forcing term");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "1.0");
  prm.leave_subsection();

  prm.enter_subsection("Desired state");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "sin(pi*x)*sin(pi*y)");
  prm.leave_subsection();

  prm.enter_subsection("Reference coefficient");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "1.0");
  prm.leave_subsection();

  prm.enter_subsection("Initial coefficient");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "1.0");
  prm.leave_subsection();

  prm.enter_subsection("Exact coefficient");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "1.0");
  prm.leave_subsection();

  prm.enter_subsection("Lower bound");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "0.2");
  prm.leave_subsection();

  prm.enter_subsection("Upper bound");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "3.0");
  prm.leave_subsection();
}



template <int dim>
void
InversePoissonKKT<dim>::parse_parameters(ParameterHandler &prm)
{
  fe = std::make_unique<FESystem<dim>>(FE_Q<dim>(state_degree),
                                       1,
                                       FE_Q<dim>(state_degree),
                                       1,
                                       FE_DGQ<dim>(control_degree),
                                       1);

  prm.enter_subsection("Forcing term");
  forcing_term.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Desired state");
  desired_state_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Reference coefficient");
  reference_coefficient_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Initial coefficient");
  initial_coefficient_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Exact coefficient");
  exact_coefficient_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Lower bound");
  lower_bound_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Upper bound");
  upper_bound_function.parse_parameters(prm);
  prm.leave_subsection();
}



template <int dim>
void
InversePoissonKKT<dim>::make_grid()
{
  triangulation.clear();
  GridGenerator::generate_from_name_and_arguments(triangulation,
                                                  grid_generator_function,
                                                  grid_generator_arguments);

  if constexpr (dim == 1)
    if (triangulation.n_active_cells() == 1)
      triangulation.begin_active()->face(1)->set_boundary_id(0);

  triangulation.refine_global(global_refinements);
}



template <int dim>
void
InversePoissonKKT<dim>::setup_system()
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
  VectorTools::interpolate_boundary_values(
    dof_handler,
    0,
    Functions::ZeroFunction<dim>(n_blocks),
    constraints,
    fe->component_mask(FEValuesExtractors::Scalar(state_block)));
  VectorTools::interpolate_boundary_values(
    dof_handler,
    0,
    Functions::ZeroFunction<dim>(n_blocks),
    constraints,
    fe->component_mask(FEValuesExtractors::Scalar(adjoint_block)));
  constraints.close();

  BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs=*/false);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);
}



template <int dim>
void
InversePoissonKKT<dim>::initialize_vectors()
{
  const std::vector<types::global_dof_index> block_sizes(dofs_per_block.begin(),
                                                         dofs_per_block.end());

  current_solution.reinit(block_sizes);
  system_rhs.reinit(block_sizes);
  residual.reinit(block_sizes);
  lower_bound.reinit(block_sizes);
  upper_bound.reinit(block_sizes);
  reference_coefficient_vector.reinit(block_sizes);
  desired_state_vector.reinit(block_sizes);
  exact_coefficient_vector.reinit(block_sizes);
}



template <int dim>
void
InversePoissonKKT<dim>::interpolate_component_function(
  const Function<dim> &function,
  const unsigned int   block,
  BlockVector<double> &vector) const
{
  vector = 0;

  const VectorFunctionFromScalarFunctionObject<dim> block_function(
    [&function](const Point<dim> &point) { return function.value(point); },
    block,
    n_blocks);

  VectorTools::interpolate(dof_handler, block_function, vector);
  constraints.distribute(vector);
}



template <int dim>
types::global_dof_index
InversePoissonKKT<dim>::control_global_offset() const
{
  return dofs_per_block[state_block] + dofs_per_block[adjoint_block];
}



template <int dim>
void
InversePoissonKKT<dim>::initialize_control_mass_matrix()
{
  AssertThrow(regularization > 0.0,
              ExcMessage("Regularization must be strictly positive."));

  DynamicSparsityPattern               dsp(dofs_per_block[control_block],
                             dofs_per_block[control_block]);
  std::vector<types::global_dof_index> local_dof_indices(fe->n_dofs_per_cell());
  const auto                           offset = control_global_offset();

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
        if (fe->system_to_block_index(i).first == control_block)
          for (unsigned int j = 0; j < fe->n_dofs_per_cell(); ++j)
            if (fe->system_to_block_index(j).first == control_block)
              dsp.add(local_dof_indices[i] - offset,
                      local_dof_indices[j] - offset);
    }

  SparsityPattern sparsity;
  sparsity.copy_from(dsp);
  control_mass_matrix.reinit(sparsity);

  const QGauss<dim> quadrature_formula(std::max(state_degree, control_degree) +
                                       1);
  FEValues<dim>     fe_values(*fe,
                          quadrature_formula,
                          update_values | update_JxW_values);
  const FEValuesExtractors::Scalar control(control_block);

  FullMatrix<double> cell_matrix(fe->n_dofs_per_cell(), fe->n_dofs_per_cell());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0;
      fe_values.reinit(cell);

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        for (const unsigned int i : fe_values.dof_indices())
          {
            const double phi_control_i = fe_values[control].value(i, q_index);

            if (phi_control_i == 0.0)
              continue;

            for (const unsigned int j : fe_values.dof_indices())
              {
                const double phi_control_j =
                  fe_values[control].value(j, q_index);

                if (phi_control_j == 0.0)
                  continue;

                cell_matrix(i, j) +=
                  phi_control_i * phi_control_j * fe_values.JxW(q_index);
              }
          }

      cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values.dof_indices())
        if (fe->system_to_block_index(i).first == control_block)
          for (const unsigned int j : fe_values.dof_indices())
            if (fe->system_to_block_index(j).first == control_block)
              control_mass_matrix.add(local_dof_indices[i] - offset,
                                      local_dof_indices[j] - offset,
                                      cell_matrix(i, j));
    }

  control_mass_inverse.initialize(control_mass_matrix);
}



template <int dim>
void
InversePoissonKKT<dim>::initialize_solution()
{
  current_solution = 0;
  interpolate_component_function(initial_coefficient_function,
                                 control_block,
                                 current_solution);
}



template <int dim>
void
InversePoissonKKT<dim>::assemble_system()
{
  system_matrix = 0;
  system_rhs    = 0;
  residual      = 0;

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

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_residual(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<double>         state_values(quadrature_formula.size());
  std::vector<Tensor<1, dim>> state_gradients(quadrature_formula.size());
  std::vector<double>         adjoint_values(quadrature_formula.size());
  std::vector<Tensor<1, dim>> adjoint_gradients(quadrature_formula.size());
  std::vector<double>         control_values(quadrature_formula.size());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix   = 0;
      cell_residual = 0;

      fe_values.reinit(cell);
      fe_values[state].get_function_values(current_solution, state_values);
      fe_values[state].get_function_gradients(current_solution,
                                              state_gradients);
      fe_values[adjoint].get_function_values(current_solution, adjoint_values);
      fe_values[adjoint].get_function_gradients(current_solution,
                                                adjoint_gradients);
      fe_values[control].get_function_values(current_solution, control_values);

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          const Point<dim> &x_q           = fe_values.quadrature_point(q_index);
          const double      forcing_value = forcing_term.value(x_q);
          const double      desired_value = desired_state_function.value(x_q);
          const double      reference_value =
            reference_coefficient_function.value(x_q);

          const double          coefficient_value = control_values[q_index];
          const double          state_value       = state_values[q_index];
          const double          adjoint_value     = adjoint_values[q_index];
          const Tensor<1, dim> &state_gradient    = state_gradients[q_index];
          const Tensor<1, dim> &adjoint_gradient  = adjoint_gradients[q_index];

          (void)adjoint_value;

          for (const unsigned int i : fe_values.dof_indices())
            {
              const double phi_state_i = fe_values[state].value(i, q_index);
              const Tensor<1, dim> grad_state_i =
                fe_values[state].gradient(i, q_index);
              const double phi_adjoint_i = fe_values[adjoint].value(i, q_index);
              const Tensor<1, dim> grad_adjoint_i =
                fe_values[adjoint].gradient(i, q_index);
              const double phi_control_i = fe_values[control].value(i, q_index);

              cell_residual(i) +=
                (coefficient_value * (state_gradient * grad_state_i) -
                 forcing_value * phi_state_i +
                 coefficient_value * (adjoint_gradient * grad_adjoint_i) +
                 (state_value - desired_value) * phi_adjoint_i +
                 (regularization * (coefficient_value - reference_value) -
                  state_gradient * adjoint_gradient) *
                   phi_control_i) *
                fe_values.JxW(q_index);

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

                  cell_matrix(i, j) +=
                    (coefficient_value * (grad_state_j * grad_state_i) +
                     phi_control_j * (state_gradient * grad_state_i) +
                     phi_state_j * phi_adjoint_i +
                     coefficient_value * (grad_adjoint_j * grad_adjoint_i) +
                     phi_control_j * (adjoint_gradient * grad_adjoint_i) +
                     regularization * phi_control_j * phi_control_i -
                     phi_control_i * (grad_state_j * adjoint_gradient) -
                     phi_control_i * (state_gradient * grad_adjoint_j)) *
                    fe_values.JxW(q_index);

                  (void)phi_adjoint_j;
                }
            }
        }

      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(
        cell_matrix, cell_residual, local_dof_indices, system_matrix, residual);
    }

  system_rhs = residual;
  system_rhs *= -1.0;
}



template <int dim>
BlockVector<double>
InversePoissonKKT<dim>::assemble_residual_only(
  const BlockVector<double> &iterate) const
{
  BlockVector<double> local_residual(iterate);
  local_residual = 0;

  const QGauss<dim> quadrature_formula(std::max(state_degree, control_degree) +
                                       1);
  FEValues<dim>     fe_values(*fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Scalar state(state_block);
  const FEValuesExtractors::Scalar adjoint(adjoint_block);
  const FEValuesExtractors::Scalar control(control_block);

  Vector<double>                       cell_residual(fe->n_dofs_per_cell());
  std::vector<types::global_dof_index> local_dof_indices(fe->n_dofs_per_cell());

  std::vector<double>         state_values(quadrature_formula.size());
  std::vector<Tensor<1, dim>> state_gradients(quadrature_formula.size());
  std::vector<double>         adjoint_values(quadrature_formula.size());
  std::vector<Tensor<1, dim>> adjoint_gradients(quadrature_formula.size());
  std::vector<double>         control_values(quadrature_formula.size());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_residual = 0;

      fe_values.reinit(cell);
      fe_values[state].get_function_values(iterate, state_values);
      fe_values[state].get_function_gradients(iterate, state_gradients);
      fe_values[adjoint].get_function_values(iterate, adjoint_values);
      fe_values[adjoint].get_function_gradients(iterate, adjoint_gradients);
      fe_values[control].get_function_values(iterate, control_values);

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          const Point<dim> &x_q           = fe_values.quadrature_point(q_index);
          const double      forcing_value = forcing_term.value(x_q);
          const double      desired_value = desired_state_function.value(x_q);
          const double      reference_value =
            reference_coefficient_function.value(x_q);

          const double          coefficient_value = control_values[q_index];
          const double          state_value       = state_values[q_index];
          const Tensor<1, dim> &state_gradient    = state_gradients[q_index];
          const Tensor<1, dim> &adjoint_gradient  = adjoint_gradients[q_index];

          (void)adjoint_values;

          for (const unsigned int i : fe_values.dof_indices())
            {
              const double phi_state_i = fe_values[state].value(i, q_index);
              const Tensor<1, dim> grad_state_i =
                fe_values[state].gradient(i, q_index);
              const double phi_adjoint_i = fe_values[adjoint].value(i, q_index);
              const Tensor<1, dim> grad_adjoint_i =
                fe_values[adjoint].gradient(i, q_index);
              const double phi_control_i = fe_values[control].value(i, q_index);

              cell_residual(i) +=
                (coefficient_value * (state_gradient * grad_state_i) -
                 forcing_value * phi_state_i +
                 coefficient_value * (adjoint_gradient * grad_adjoint_i) +
                 (state_value - desired_value) * phi_adjoint_i +
                 (regularization * (coefficient_value - reference_value) -
                  state_gradient * adjoint_gradient) *
                   phi_control_i) *
                fe_values.JxW(q_index);
            }
        }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_residual,
                                             local_dof_indices,
                                             local_residual);
    }

  return local_residual;
}



template <int dim>
Vector<double>
InversePoissonKKT<dim>::compute_stationarity() const
{
  Vector<double> dual_stationarity(dofs_per_block[control_block]);
  dual_stationarity = residual.block(control_block);

  Vector<double> stationarity(dofs_per_block[control_block]);
  control_mass_inverse.vmult(stationarity, dual_stationarity);

  return stationarity;
}



template <int dim>
std::vector<bool>
InversePoissonKKT<dim>::compute_lower_active_set() const
{
  const auto &control_lower = lower_bound.block(control_block);
  const auto &control_value = current_solution.block(control_block);
  const auto  stationarity  = compute_stationarity();

  std::vector<bool> lower_active(control_value.size(), false);

  for (unsigned int i = 0; i < control_value.size(); ++i)
    {
      const double indicator =
        stationarity[i] +
        active_set_parameter * (control_value[i] - control_lower[i]);
      if (indicator < -newton_tolerance)
        lower_active[i] = true;
    }

  return lower_active;
}



template <int dim>
std::vector<bool>
InversePoissonKKT<dim>::compute_upper_active_set() const
{
  const auto &control_upper = upper_bound.block(control_block);
  const auto &control_value = current_solution.block(control_block);
  const auto  stationarity  = compute_stationarity();

  std::vector<bool> upper_active(control_value.size(), false);

  for (unsigned int i = 0; i < control_value.size(); ++i)
    {
      const double indicator =
        stationarity[i] +
        active_set_parameter * (control_value[i] - control_upper[i]);
      if (indicator > newton_tolerance)
        upper_active[i] = true;
    }

  return upper_active;
}



template <int dim>
AffineConstraints<double>
InversePoissonKKT<dim>::make_control_constraints(
  const std::vector<bool> &lower_active,
  const std::vector<bool> &upper_active) const
{
  AssertDimension(lower_active.size(), upper_active.size());

  AffineConstraints<double> control_constraints;
  const auto                offset = control_global_offset();

  for (unsigned int i = 0; i < lower_active.size(); ++i)
    if (lower_active[i] || upper_active[i])
      {
        control_constraints.add_line(offset + i);
        control_constraints.set_inhomogeneity(
          offset + i,
          lower_active[i] ? lower_bound.block(control_block)[i] -
                              current_solution.block(control_block)[i] :
                            upper_bound.block(control_block)[i] -
                              current_solution.block(control_block)[i]);
      }

  control_constraints.close();

  AffineConstraints<double> combined_constraints;
  combined_constraints.copy_from(constraints);
  combined_constraints.merge(control_constraints);
  return combined_constraints;
}



template <int dim>
void
InversePoissonKKT<dim>::solve_linearized_system(
  const std::vector<bool> &lower_active,
  const std::vector<bool> &upper_active,
  BlockVector<double>     &update) const
{
  const AffineConstraints<double> combined_constraints =
    make_control_constraints(lower_active, upper_active);

  BlockSparseMatrix<double> modified_matrix;
  modified_matrix.reinit(sparsity_pattern);
  modified_matrix.copy_from(system_matrix);

  BlockVector<double> modified_rhs(system_rhs);

  combined_constraints.condense(modified_matrix, modified_rhs);

  SparseDirectUMFPACK direct_solver;
  direct_solver.initialize(modified_matrix);
  direct_solver.vmult(update, modified_rhs);
  combined_constraints.distribute(update);
}



template <int dim>
void
InversePoissonKKT<dim>::enforce_box_constraints(
  BlockVector<double>     &iterate,
  const std::vector<bool> &lower_active,
  const std::vector<bool> &upper_active) const
{
  constraints.distribute(iterate);

  for (unsigned int i = 0; i < lower_active.size(); ++i)
    if (lower_active[i])
      iterate.block(control_block)[i] = lower_bound.block(control_block)[i];
    else if (upper_active[i])
      iterate.block(control_block)[i] = upper_bound.block(control_block)[i];

  for (unsigned int i = 0; i < lower_active.size(); ++i)
    if (!lower_active[i] && !upper_active[i])
      iterate.block(control_block)[i] =
        std::clamp(iterate.block(control_block)[i],
                   lower_bound.block(control_block)[i],
                   upper_bound.block(control_block)[i]);
}



template <int dim>
std::string
InversePoissonKKT<dim>::output_iteration(
  const unsigned int         iteration,
  const BlockVector<double> &iterate,
  const std::vector<bool>   &lower_active,
  const std::vector<bool>   &upper_active) const
{
  const std::string filename =
    output_name + "-" + Utilities::int_to_string(iteration, 4);
  output_results(iterate, lower_active, upper_active, filename);
  return filename + ".vtu";
}



template <int dim>
void
InversePoissonKKT<dim>::output_results(const BlockVector<double> &iterate,
                                       const std::vector<bool>   &lower_active,
                                       const std::vector<bool>   &upper_active,
                                       const std::string &filename) const
{
  BlockVector<double> lower_active_output(iterate);
  BlockVector<double> upper_active_output(iterate);
  BlockVector<double> inactive_output(iterate);
  lower_active_output = 0;
  upper_active_output = 0;
  inactive_output     = 0;

  for (unsigned int i = 0; i < lower_active.size(); ++i)
    if (lower_active[i])
      lower_active_output.block(control_block)[i] = 1.0;
    else if (upper_active[i])
      upper_active_output.block(control_block)[i] = 1.0;
    else
      inactive_output.block(control_block)[i] = 1.0;

  std::vector<std::string> solution_names = {"state", "adjoint", "diffusion"};
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(n_blocks, DataComponentInterpretation::component_is_scalar);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(dof_handler,
                           iterate,
                           solution_names,
                           interpretation);
  data_out.add_data_vector(dof_handler,
                           desired_state_vector,
                           {"desired_state_projection",
                            "desired_state_dummy_adjoint",
                            "desired_state_dummy_control"},
                           interpretation);
  data_out.add_data_vector(dof_handler,
                           reference_coefficient_vector,
                           {"reference_dummy_state",
                            "reference_dummy_adjoint",
                            "reference_diffusion"},
                           interpretation);
  data_out.add_data_vector(dof_handler,
                           exact_coefficient_vector,
                           {"exact_dummy_state",
                            "exact_dummy_adjoint",
                            "exact_diffusion"},
                           interpretation);
  data_out.add_data_vector(dof_handler,
                           lower_bound,
                           {"lower_dummy_state",
                            "lower_dummy_adjoint",
                            "lower_bound"},
                           interpretation);
  data_out.add_data_vector(dof_handler,
                           upper_bound,
                           {"upper_dummy_state",
                            "upper_dummy_adjoint",
                            "upper_bound"},
                           interpretation);
  data_out.add_data_vector(dof_handler,
                           lower_active_output,
                           {"lower_active_dummy_state",
                            "lower_active_dummy_adjoint",
                            "lower_active"},
                           interpretation);
  data_out.add_data_vector(dof_handler,
                           upper_active_output,
                           {"upper_active_dummy_state",
                            "upper_active_dummy_adjoint",
                            "upper_active"},
                           interpretation);
  data_out.add_data_vector(dof_handler,
                           inactive_output,
                           {"inactive_dummy_state",
                            "inactive_dummy_adjoint",
                            "inactive"},
                           interpretation);
  data_out.build_patches();

  std::ofstream output(filename + ".vtu");
  data_out.write_vtu(output);
}



template <int dim>
void
InversePoissonKKT<dim>::run()
{
  make_grid();
  setup_system();
  initialize_vectors();

  interpolate_component_function(lower_bound_function,
                                 control_block,
                                 lower_bound);
  interpolate_component_function(upper_bound_function,
                                 control_block,
                                 upper_bound);
  interpolate_component_function(reference_coefficient_function,
                                 control_block,
                                 reference_coefficient_vector);
  interpolate_component_function(desired_state_function,
                                 state_block,
                                 desired_state_vector);
  interpolate_component_function(exact_coefficient_function,
                                 control_block,
                                 exact_coefficient_vector);

  const auto &control_lower = lower_bound.block(control_block);
  const auto &control_upper = upper_bound.block(control_block);
  for (unsigned int i = 0; i < control_lower.size(); ++i)
    AssertThrow(control_lower[i] <= control_upper[i],
                ExcMessage("Lower bound exceeds upper bound."));

  initialize_solution();
  enforce_box_constraints(current_solution,
                          std::vector<bool>(control_lower.size(), false),
                          std::vector<bool>(control_lower.size(), false));
  initialize_control_mass_matrix();

  std::vector<bool> lower_active(control_lower.size(), false);
  std::vector<bool> upper_active(control_lower.size(), false);
  std::vector<std::pair<double, std::string>> times_and_names;

  for (unsigned int iteration = 0; iteration < max_pdas_iterations; ++iteration)
    {
      assemble_system();
      const double current_residual_norm = residual.l2_norm();

      const auto next_lower_active = compute_lower_active_set();
      const auto next_upper_active = compute_upper_active_set();

      BlockVector<double> newton_update(current_solution);
      newton_update = 0;
      solve_linearized_system(next_lower_active,
                              next_upper_active,
                              newton_update);

      BlockVector<double> next_solution(current_solution);
      BlockVector<double> next_residual(current_solution);
      next_residual        = residual;
      double accepted_step = 0.0;
      double residual_norm = std::numeric_limits<double>::infinity();
      double best_norm     = current_residual_norm;
      double best_step     = 0.0;
      auto   best_solution = current_solution;
      auto   best_residual = residual;

      for (double step_length = 1.0; step_length >= 1.0 / 1024.0;
           step_length *= 0.5)
        {
          BlockVector<double> trial_solution(current_solution);
          trial_solution.add(step_length, newton_update);
          enforce_box_constraints(trial_solution,
                                  next_lower_active,
                                  next_upper_active);

          const BlockVector<double> trial_residual =
            assemble_residual_only(trial_solution);
          const double trial_norm = trial_residual.l2_norm();

          if (trial_norm < best_norm)
            {
              best_norm     = trial_norm;
              best_step     = step_length;
              best_solution = trial_solution;
              best_residual = trial_residual;
            }

          if (trial_norm < current_residual_norm || trial_norm < 1e-12)
            {
              next_solution = trial_solution;
              next_residual = trial_residual;
              residual_norm = trial_norm;
              accepted_step = step_length;
              break;
            }
        }

      if (accepted_step == 0.0)
        {
          next_solution = best_solution;
          next_residual = best_residual;
          residual_norm = best_norm;
          accepted_step = best_step;
        }

      BlockVector<double> update(next_solution);
      update -= current_solution;
      const double update_norm = update.l2_norm();

      const unsigned int n_lower =
        std::count(next_lower_active.begin(), next_lower_active.end(), true);
      const unsigned int n_upper =
        std::count(next_upper_active.begin(), next_upper_active.end(), true);

      std::cout << "PDAS/Newton iteration " << iteration
                << ": |A_-|=" << n_lower << ", |A_+|=" << n_upper
                << ", step=" << accepted_step << ", ||delta||=" << update_norm
                << ", ||R||=" << residual_norm << std::endl;

      times_and_names.emplace_back(static_cast<double>(iteration),
                                   output_iteration(iteration,
                                                    next_solution,
                                                    next_lower_active,
                                                    next_upper_active));

      if (next_lower_active == lower_active &&
          next_upper_active == upper_active && update_norm < newton_tolerance &&
          residual_norm < newton_tolerance)
        {
          current_solution = next_solution;
          residual         = next_residual;
          std::cout << "Converged after " << iteration + 1
                    << " PDAS/Newton iterations." << std::endl;
          break;
        }

      if (next_lower_active == lower_active &&
          next_upper_active == upper_active && accepted_step == 0.0)
        {
          current_solution = next_solution;
          residual         = next_residual;
          std::cout << "Stopped because the damped Newton step no longer"
                    << " improves the KKT residual." << std::endl;
          break;
        }

      current_solution = next_solution;
      residual         = next_residual;
      lower_active     = next_lower_active;
      upper_active     = next_upper_active;

      if (iteration + 1 == max_pdas_iterations)
        std::cout << "Reached the maximum number of PDAS iterations."
                  << std::endl;
    }

  output_results(current_solution, lower_active, upper_active, output_name);

  std::ofstream pvd_output(output_name + ".pvd");
  DataOutBase::write_pvd_record(pvd_output, times_and_names);
}



template class InversePoissonKKT<DEAL_DIMENSION>;

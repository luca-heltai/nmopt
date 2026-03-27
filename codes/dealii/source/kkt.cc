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

#include "kkt.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <vector>

using namespace dealii;


template <int dim>
KKT<dim>::KKT(const std::string &subsection_name)
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
  , regularization(1e-4)
  , grid_generator_function("hyper_cube")
  , grid_generator_arguments("-1.0 : 1.0 : false")
  , diffusion_coefficient(1)
  , reaction_coefficient(1)
  , transport_field(dim)
  , forcing_term(1)
  , desired_state(1)
  , state_boundary_values(1)
  , adjoint_boundary_values(1)
{}



template <int dim>
void
KKT<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.add_parameter("State degree", state_degree);
  prm.add_parameter("Control degree", control_degree);
  prm.add_parameter("Continuous control", continuous_control);
  prm.add_parameter("Global refinements", global_refinements);
  prm.add_parameter("Regularization", regularization);
  prm.add_parameter("Grid generator function", grid_generator_function);
  prm.add_parameter("Grid generator arguments", grid_generator_arguments);

  prm.enter_subsection("Diffusion coefficient");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "1");
  prm.leave_subsection();

  prm.enter_subsection("Reaction coefficient");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();

  prm.enter_subsection("Transport field");
  Functions::ParsedFunction<dim>::declare_parameters(prm, dim);
  prm.leave_subsection();

  prm.enter_subsection("Forcing term");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();

  prm.enter_subsection("Desired state");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "x^2+y^2 < .5^2 ? 1 : 0");
  prm.leave_subsection();

  prm.enter_subsection("State boundary values");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();

  prm.enter_subsection("Adjoint boundary values");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();
}



template <int dim>
void
KKT<dim>::parse_parameters(ParameterHandler &prm)
{
  fe = create_fe_system();

  prm.enter_subsection("Diffusion coefficient");
  diffusion_coefficient.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Reaction coefficient");
  reaction_coefficient.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Transport field");
  transport_field.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Forcing term");
  forcing_term.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Desired state");
  desired_state.parse_parameters(prm);
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
KKT<dim>::create_fe_system() const
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
KKT<dim>::make_grid()
{
  triangulation.clear();
  GridGenerator::generate_from_name_and_arguments(triangulation,
                                                  grid_generator_function,
                                                  grid_generator_arguments);
  triangulation.refine_global(global_refinements);
}



template <int dim>
void
KKT<dim>::initialize()
{
  make_grid();
  setup_system();
  assemble_system();
}



template <int dim>
void
KKT<dim>::setup_system()
{
  dof_handler.distribute_dofs(*fe);

  const std::vector<unsigned int> block_component = {state_block,
                                                     adjoint_block,
                                                     control_block};
  DoFRenumbering::component_wise(dof_handler, block_component);

  dofs_per_block = DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
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

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(std::vector<types::global_dof_index>(dofs_per_block.begin(),
                                                       dofs_per_block.end()));
  system_rhs.reinit(std::vector<types::global_dof_index>(dofs_per_block.begin(),
                                                         dofs_per_block.end()));
}



template <int dim>
void
KKT<dim>::assemble_system()
{
  system_matrix = 0;
  system_rhs    = 0;

  const QGauss<dim> quadrature_formula(
    std::max(state_degree, control_degree) + 1);
  FEValues<dim>     fe_values(*fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Scalar state(state_block);
  const FEValuesExtractors::Scalar adjoint(adjoint_block);
  const FEValuesExtractors::Scalar control(control_block);

  const unsigned int dofs_per_cell = fe->n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  Vector<double>                       transport_values(dim);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0;
      cell_rhs    = 0;

      fe_values.reinit(cell);

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          const Point<dim> &x_q       = fe_values.quadrature_point(q_index);
          const double      diffusion = diffusion_coefficient.value(x_q);
          const double      reaction  = reaction_coefficient.value(x_q);
          const double      rhs       = forcing_term.value(x_q);
          const double      target    = desired_state.value(x_q);

          transport_field.vector_value(x_q, transport_values);
          Tensor<1, dim> transport;
          for (unsigned int d = 0; d < dim; ++d)
            transport[d] = transport_values[d];

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

                  // Tracking term on the state variable.
                  cell_matrix(i, j) +=
                    phi_state_i * phi_state_j * fe_values.JxW(q_index);

                  // State equation block A(y,p).
                  cell_matrix(i, j) +=
                    (diffusion * grad_state_i * grad_adjoint_j +
                     (transport * grad_state_i) * phi_adjoint_j +
                     reaction * phi_state_i * phi_adjoint_j) *
                    fe_values.JxW(q_index);

                  // Adjoint equation block A'(p,y).
                  cell_matrix(i, j) +=
                    (diffusion * grad_adjoint_i * grad_state_j +
                     (transport * grad_state_j) * phi_adjoint_i +
                     reaction * phi_state_j * phi_adjoint_i) *
                    fe_values.JxW(q_index);

                  // Optimality block coupling adjoint and control.
                  cell_matrix(i, j) +=
                    (-phi_adjoint_i * phi_control_j -
                     phi_control_i * phi_adjoint_j +
                     regularization * phi_control_i * phi_control_j) *
                    fe_values.JxW(q_index);
                }

              // Right-hand side for tracking and forcing.
              cell_rhs(i) += (target * phi_state_i + rhs * phi_adjoint_i) *
                             fe_values.JxW(q_index);
            }
        }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
}



template <int dim>
void
KKT<dim>::output_results(const std::string &filename) const
{
  std::vector<std::string> solution_names = {"state", "adjoint", "control"};
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(3, DataComponentInterpretation::component_is_scalar);
  BlockVector<double> desired_state_output;
  desired_state_output.reinit(std::vector<types::global_dof_index>(
    dofs_per_block.begin(), dofs_per_block.end()));
  const VectorFunctionFromScalarFunctionObject<dim> desired_state_function(
    [this](const Point<dim> &p) { return desired_state.value(p); },
    state_block,
    n_blocks);

  VectorTools::interpolate(dof_handler, desired_state_function, desired_state_output);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(dof_handler,
                           solution,
                           solution_names,
                           interpretation);
  data_out.add_data_vector(dof_handler,
                           desired_state_output,
                           {"desired_state", "zero_adjoint", "zero_control"},
                           interpretation);
  data_out.build_patches();

  std::ofstream output(filename + ".vtu");
  data_out.write_vtu(output);
}



template <int dim>
void
KKT<dim>::solve()
{
  SparseDirectUMFPACK direct_solver;
  direct_solver.initialize(system_matrix);
  direct_solver.vmult(solution, system_rhs);
  constraints.distribute(solution);
}



template <int dim>
void
KKT<dim>::run()
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

  solve();
}



template <int dim>
const BlockSparseMatrix<double> &
KKT<dim>::get_system_matrix() const
{
  return system_matrix;
}



template <int dim>
const SparseMatrix<double> &
KKT<dim>::get_system_block(const unsigned int row,
                           const unsigned int column) const
{
  AssertIndexRange(row, n_blocks);
  AssertIndexRange(column, n_blocks);
  return system_matrix.block(row, column);
}



template <int dim>
const AffineConstraints<double> &
KKT<dim>::get_constraints() const
{
  return constraints;
}



template <int dim>
const DoFHandler<dim> &
KKT<dim>::get_dof_handler() const
{
  return dof_handler;
}



template <int dim>
const BlockVector<double> &
KKT<dim>::get_system_rhs() const
{
  return system_rhs;
}



template <int dim>
const BlockVector<double> &
KKT<dim>::get_solution() const
{
  return solution;
}



template <int dim>
const std::vector<types::global_dof_index> &
KKT<dim>::get_dofs_per_block() const
{
  return dofs_per_block;
}


template class KKT<DEAL_DIMENSION>;

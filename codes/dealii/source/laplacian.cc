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

#include "laplacian.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>

using namespace dealii;


template <int dim>
Laplacian<dim>::Laplacian(const std::string &subsection_name)
  : ParameterAcceptor(subsection_name)
  , fe(std::make_unique<FE_Q<dim>>(1))
  , dof_handler(triangulation)
  , fe_degree(1)
  , global_refinements(5)
  , solver_max_iterations(1000)
  , solver_tolerance(1e-12)
  , preconditioner_relaxation(1.2)
  , grid_generator_function("hyper_cube")
  , grid_generator_arguments("0.0 : 1.0 : false")
  , diffusion_coefficient(1)
  , forcing_term(1)
  , boundary_values(1)
{}



template <int dim>
void
Laplacian<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.add_parameter("FE degree", fe_degree);
  prm.add_parameter("Global refinements", global_refinements);
  prm.add_parameter("Solver max iterations", solver_max_iterations);
  prm.add_parameter("Solver tolerance", solver_tolerance);
  prm.add_parameter("Preconditioner relaxation", preconditioner_relaxation);
  prm.add_parameter("Grid generator function", grid_generator_function);
  prm.add_parameter("Grid generator arguments", grid_generator_arguments);

  prm.enter_subsection("Diffusion coefficient");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();

  prm.enter_subsection("Forcing term");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();

  prm.enter_subsection("Boundary values");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();
}



template <int dim>
void
Laplacian<dim>::parse_parameters(ParameterHandler &prm)
{
  fe = std::make_unique<FE_Q<dim>>(fe_degree);

  prm.enter_subsection("Diffusion coefficient");
  diffusion_coefficient.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Forcing term");
  forcing_term.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Boundary values");
  boundary_values.parse_parameters(prm);
  prm.leave_subsection();
}



template <int dim>
void
Laplacian<dim>::make_grid()
{
  triangulation.clear();
  GridGenerator::generate_from_name_and_arguments(triangulation,
                                                  grid_generator_function,
                                                  grid_generator_arguments);
  triangulation.refine_global(global_refinements);
}



template <int dim>
void
Laplacian<dim>::initialize()
{
  make_grid();
  setup_system();
  assemble_system();
}



template <int dim>
void
Laplacian<dim>::setup_system()
{
  dof_handler.distribute_dofs(*fe);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           boundary_values,
                                           constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs = */ false);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  mass_matrix.reinit(sparsity_pattern);
}



template <int dim>
void
Laplacian<dim>::assemble_system()
{
  system_matrix = 0;
  mass_matrix   = 0;
  system_rhs    = 0;

  const QGauss<dim> quadrature_formula(fe->degree + 1);
  FEValues<dim>     fe_values(*fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe->n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix      = 0;
      cell_mass_matrix = 0;
      cell_rhs         = 0;

      fe_values.reinit(cell);

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          const double current_coefficient =
            diffusion_coefficient.value(fe_values.quadrature_point(q_index));
          const double current_rhs =
            forcing_term.value(fe_values.quadrature_point(q_index));

          for (const unsigned int i : fe_values.dof_indices())
            {
              for (const unsigned int j : fe_values.dof_indices())
                {
                  cell_matrix(i, j) +=
                    current_coefficient * fe_values.shape_grad(i, q_index) *
                    fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index);
                  cell_mass_matrix(i, j) += fe_values.shape_value(i, q_index) *
                                            fe_values.shape_value(j, q_index) *
                                            fe_values.JxW(q_index);
                }

              cell_rhs(i) += current_rhs * fe_values.shape_value(i, q_index) *
                             fe_values.JxW(q_index);
            }
        }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      constraints.distribute_local_to_global(cell_mass_matrix,
                                             local_dof_indices,
                                             mass_matrix);
    }
}



template <int dim>
void
Laplacian<dim>::solve()
{
  solve(system_rhs, solution);
}



template <int dim>
void
Laplacian<dim>::solve(const Vector<double> &rhs, Vector<double> &dst) const
{
  dst.reinit(solution);

  SolverControl solver_control(solver_max_iterations, solver_tolerance);
  SolverCG<Vector<double>> solver(solver_control);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, preconditioner_relaxation);

  solver.solve(system_matrix, dst, rhs, preconditioner);
  constraints.distribute(dst);
}



template <int dim>
void
Laplacian<dim>::output_results(const std::string &filename) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  std::ofstream output(filename + ".vtu");
  data_out.write_vtu(output);
}



template <int dim>
void
Laplacian<dim>::run()
{
  initialize();

  std::cout << "Number of active cells:       "
            << triangulation.n_active_cells() << std::endl;
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  solve();
  output_results();
}



template <int dim>
const SparseMatrix<double> &
Laplacian<dim>::get_mass_matrix() const
{
  return mass_matrix;
}



template <int dim>
const SparseMatrix<double> &
Laplacian<dim>::get_system_matrix() const
{
  return system_matrix;
}



template <int dim>
const AffineConstraints<double> &
Laplacian<dim>::get_constraints() const
{
  return constraints;
}



template <int dim>
const DoFHandler<dim> &
Laplacian<dim>::get_dof_handler() const
{
  return dof_handler;
}



template <int dim>
const Vector<double> &
Laplacian<dim>::get_system_rhs() const
{
  return system_rhs;
}



template <int dim>
const Vector<double> &
Laplacian<dim>::get_solution() const
{
  return solution;
}


template class Laplacian<DEAL_DIMENSION>;

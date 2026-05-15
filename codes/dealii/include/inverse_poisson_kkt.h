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

#ifndef dealii_inverse_poisson_kkt_h
#define dealii_inverse_poisson_kkt_h

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <memory>
#include <string>
#include <vector>

using namespace dealii;

/**
 * One-shot KKT solver for a coefficient-identification problem for Poisson's
 * equation with homogeneous Dirichlet boundary conditions. The control is the
 * diffusion coefficient, discretized with FE_DGQ on the same triangulation as
 * state and adjoint variables.
 *
 * The optimization problem is
 *
 *   min 1/2 ||y-y_d||^2 + alpha/2 ||a-a_ref||^2
 *
 * subject to
 *
 *   -div(a grad y) = f in Omega,
 *   y = 0 on boundary(Omega),
 *   a_min <= a <= a_max.
 *
 * Box constraints are enforced with a primal-dual active-set strategy. On each
 * PDAS iteration we assemble the Jacobian of the all-at-once KKT system and
 * solve the resulting linearized system with AffineConstraints<double>.
 */
template <int dim>
class InversePoissonKKT : public ParameterAcceptor
{
public:
  static constexpr unsigned int n_blocks      = 3;
  static constexpr unsigned int state_block   = 0;
  static constexpr unsigned int adjoint_block = 1;
  static constexpr unsigned int control_block = 2;

  InversePoissonKKT(const std::string &subsection_name =
                      "Inverse Poisson KKT problem");

  void
  run();

protected:
  void
  declare_parameters(ParameterHandler &prm) override;

  void
  parse_parameters(ParameterHandler &prm) override;

private:
  void
  make_grid();

  void
  setup_system();

  void
  initialize_vectors();

  void
  initialize_control_mass_matrix();

  void
  initialize_solution();

  void
  assemble_system();

  BlockVector<double>
  assemble_residual_only(const BlockVector<double> &iterate) const;

  std::vector<bool>
  compute_lower_active_set() const;

  std::vector<bool>
  compute_upper_active_set() const;

  Vector<double>
  compute_stationarity() const;

  AffineConstraints<double>
  make_control_constraints(const std::vector<bool> &lower_active,
                           const std::vector<bool> &upper_active) const;

  void
  solve_linearized_system(const std::vector<bool> &lower_active,
                          const std::vector<bool> &upper_active,
                          BlockVector<double>     &update) const;

  void
  enforce_box_constraints(BlockVector<double>     &iterate,
                          const std::vector<bool> &lower_active,
                          const std::vector<bool> &upper_active) const;

  void
  interpolate_component_function(const Function<dim> &function,
                                 const unsigned int   block,
                                 BlockVector<double> &vector) const;

  void
  output_results(const BlockVector<double>     &iterate,
                 const std::vector<bool>       &lower_active,
                 const std::vector<bool>       &upper_active,
                 const std::string             &filename) const;

  std::string
  output_iteration(const unsigned int          iteration,
                   const BlockVector<double>   &iterate,
                   const std::vector<bool>     &lower_active,
                   const std::vector<bool>     &upper_active) const;

  types::global_dof_index
  control_global_offset() const;

  Triangulation<dim> triangulation;

  std::unique_ptr<FESystem<dim>> fe;
  DoFHandler<dim>                dof_handler;

  AffineConstraints<double> constraints;

  BlockSparsityPattern      sparsity_pattern;
  BlockSparseMatrix<double> system_matrix;
  BlockVector<double>       current_solution;
  BlockVector<double>       system_rhs;
  BlockVector<double>       residual;
  BlockVector<double>       lower_bound;
  BlockVector<double>       upper_bound;
  BlockVector<double>       reference_coefficient_vector;
  BlockVector<double>       desired_state_vector;
  BlockVector<double>       exact_coefficient_vector;

  SparseMatrix<double> control_mass_matrix;
  SparseDirectUMFPACK  control_mass_inverse;

  std::vector<types::global_dof_index> dofs_per_block;

  unsigned int state_degree;
  unsigned int control_degree;
  unsigned int global_refinements;
  unsigned int max_pdas_iterations;
  double       regularization;
  double       active_set_parameter;
  double       newton_tolerance;
  std::string  grid_generator_function;
  std::string  grid_generator_arguments;
  std::string  output_name;

  Functions::ParsedFunction<dim> forcing_term;
  Functions::ParsedFunction<dim> desired_state_function;
  Functions::ParsedFunction<dim> reference_coefficient_function;
  Functions::ParsedFunction<dim> initial_coefficient_function;
  Functions::ParsedFunction<dim> exact_coefficient_function;
  Functions::ParsedFunction<dim> lower_bound_function;
  Functions::ParsedFunction<dim> upper_bound_function;
};

#endif

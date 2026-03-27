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

#ifndef dealii_kkt_h
#define dealii_kkt_h

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

#include <fstream>
#include <memory>

using namespace dealii;

/**
 * Assemble and solve the all-at-once KKT system for a linear-quadratic
 * PDE-constrained optimization problem.
 *
 * The code looks for a state y, an adjoint p, and a control u minimizing
 *
 *   1/2 ||y-y_d||^2 + alpha/2 ||u||^2
 *
 * subject to the diffusion-reaction-transport equation
 *
 *   -div(mu grad y) + beta . grad y + sigma y = u + f
 *
 * with Dirichlet boundary conditions on state and adjoint.
 *
 * The discrete unknown is ordered in three blocks:
 *
 * - block 0: state
 * - block 1: adjoint
 * - block 2: control
 *
 * and the class assembles the full block KKT matrix so that one can access
 * both the monolithic system and its individual blocks.
 */
template <int dim>
class KKT : public ParameterAcceptor
{
public:
  static constexpr unsigned int n_blocks      = 3;
  static constexpr unsigned int state_block   = 0;
  static constexpr unsigned int adjoint_block = 1;
  static constexpr unsigned int control_block = 2;

  KKT(const std::string &subsection_name = "KKT problem");

  void
  initialize();

  void
  solve();

  void
  output_results(const std::string &filename = "kkt_solution") const;

  void
  run();

  const BlockSparseMatrix<double> &
  get_system_matrix() const;

  const SparseMatrix<double> &
  get_system_block(const unsigned int row, const unsigned int column) const;

  const AffineConstraints<double> &
  get_constraints() const;

  const DoFHandler<dim> &
  get_dof_handler() const;

  const BlockVector<double> &
  get_system_rhs() const;

  const BlockVector<double> &
  get_solution() const;

  const std::vector<types::global_dof_index> &
  get_dofs_per_block() const;

protected:
  void
  declare_parameters(ParameterHandler &prm) override;

  void
  parse_parameters(ParameterHandler &prm) override;

  void
  make_grid();

  void
  setup_system();

  void
  assemble_system();

  std::unique_ptr<FESystem<dim>>
  create_fe_system() const;

  Triangulation<dim> triangulation;

  std::unique_ptr<FESystem<dim>> fe;
  DoFHandler<dim>                dof_handler;

  AffineConstraints<double> constraints;

  BlockSparsityPattern      sparsity_pattern;
  BlockSparseMatrix<double> system_matrix;
  BlockVector<double>       solution;
  BlockVector<double>       system_rhs;

  std::vector<types::global_dof_index> dofs_per_block;

  unsigned int state_degree;
  unsigned int control_degree;
  bool         continuous_control;
  unsigned int global_refinements;
  double       regularization;
  std::string  grid_generator_function;
  std::string  grid_generator_arguments;

  Functions::ParsedFunction<dim> diffusion_coefficient;
  Functions::ParsedFunction<dim> reaction_coefficient;
  Functions::ParsedFunction<dim> transport_field;
  Functions::ParsedFunction<dim> forcing_term;
  Functions::ParsedFunction<dim> desired_state;
  Functions::ParsedFunction<dim> state_boundary_values;
  Functions::ParsedFunction<dim> adjoint_boundary_values;
};

#endif

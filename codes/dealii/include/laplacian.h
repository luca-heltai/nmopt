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

#ifndef dealii_laplacian_h
#define dealii_laplacian_h

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <memory>

using namespace dealii;

template <int dim>
class Laplacian : public ParameterAcceptor
{
public:
  Laplacian(const std::string &subsection_name = "Poisson problem");

  void
  initialize();

  void
  solve();

  void
  solve(const Vector<double> &rhs, Vector<double> &dst) const;

  void
  output_results(const std::string &filename = "solution") const;

  void
  run();

  const SparseMatrix<double> &
  get_mass_matrix() const;

  const SparseMatrix<double> &
  get_system_matrix() const;

  const AffineConstraints<double> &
  get_constraints() const;

  const DoFHandler<dim> &
  get_dof_handler() const;

  const Vector<double> &
  get_system_rhs() const;

  const Vector<double> &
  get_solution() const;

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

  Triangulation<dim> triangulation;

  std::unique_ptr<FE_Q<dim>> fe;
  DoFHandler<dim>            dof_handler;


  AffineConstraints<double> constraints;

  SparseMatrix<double> system_matrix;
  SparseMatrix<double> mass_matrix;
  SparsityPattern      sparsity_pattern;

  Vector<double> solution;
  Vector<double> system_rhs;

  unsigned int fe_degree;
  unsigned int global_refinements;
  unsigned int solver_max_iterations;
  double       solver_tolerance;
  double       preconditioner_relaxation;
  std::string  grid_generator_function;
  std::string  grid_generator_arguments;

  Functions::ParsedFunction<dim> diffusion_coefficient;
  Functions::ParsedFunction<dim> forcing_term;
  Functions::ParsedFunction<dim> boundary_values;
};

#endif

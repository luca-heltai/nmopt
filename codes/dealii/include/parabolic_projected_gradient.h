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

#ifndef dealii_parabolic_projected_gradient_h
#define dealii_parabolic_projected_gradient_h

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
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "optimization_tools.h"

using namespace dealii;

/**
 * Reduced projected-gradient solver for a linear-quadratic parabolic optimal
 * control problem.
 *
 * The time discretization is implicit Euler.  For a fixed control sequence the
 * state is propagated forward in time, the adjoint is propagated backward in
 * time, and the control is updated by a box-projected gradient step.
 *
 * The finite element system has three components:
 *
 * - block 0: state space;
 * - block 1: adjoint space;
 * - block 2: control space.
 *
 * The state and adjoint spaces are continuous FE_Q spaces.  The control space
 * can be either continuous FE_Q or discontinuous FE_DGQ, as in the elliptic
 * laboratories.
 */
template <int dim>
class ParabolicProjectedGradient : public ParameterAcceptor
{
public:
  static constexpr unsigned int n_blocks      = 3;
  static constexpr unsigned int state_block   = 0;
  static constexpr unsigned int adjoint_block = 1;
  static constexpr unsigned int control_block = 2;

  struct Diagnostics
  {
    double cost                    = 0.0;
    double tracking_cost           = 0.0;
    double terminal_tracking_cost  = 0.0;
    double control_cost            = 0.0;
    double gradient_norm           = 0.0;
    double projected_gradient_norm = 0.0;
    double control_norm            = 0.0;
    double step_length             = 0.0;
  };

  ParabolicProjectedGradient(const std::string &subsection_name =
                               "Parabolic projected-gradient problem");

  void
  initialize();

  void
  run();

protected:
  void
  declare_parameters(ParameterHandler &prm) override;

  void
  parse_parameters(ParameterHandler &prm) override;

private:
  using TimeVector = std::vector<Vector<double>>;

  void
  make_grid();

  void
  setup_system();

  void
  assemble_matrices();

  void
  initialize_time_step_solvers();

  std::unique_ptr<FESystem<dim>>
  create_fe_system() const;

  void
  reinit_time_vector(TimeVector        &values,
                     const unsigned int block,
                     const bool         include_terminal_slot = false) const;

  void
  interpolate_function_to_block(Functions::ParsedFunction<dim> &function,
                                const unsigned int              block,
                                const double                    time,
                                Vector<double>                 &vector) const;

  void
  distribute_block(const unsigned int block, Vector<double> &vector) const;

  void
  initialize_time_data();

  void
  initialize_control(TimeVector &control) const;

  void
  project_onto_box(TimeVector &control) const;

  void
  project_time_level(const unsigned int time_index,
                     Vector<double>    &control) const;

  void
  solve_state(const TimeVector &control, TimeVector &state) const;

  void
  solve_adjoint(const TimeVector &state, TimeVector &adjoint) const;

  void
  compute_gradient(const TimeVector &control,
                   const TimeVector &adjoint,
                   TimeVector       &gradient) const;

  void
  projected_gradient_step(TimeVector       &trial_control,
                          const TimeVector &gradient,
                          const double      step_length) const;

  double
  control_inner_product(const TimeVector &left, const TimeVector &right) const;

  double
  control_norm(const TimeVector &control) const;

  Diagnostics
  compute_diagnostics(const TimeVector &state,
                      const TimeVector &control,
                      const TimeVector &gradient,
                      const double      step_length) const;

  double
  compute_projected_gradient_norm(const TimeVector &control,
                                  const TimeVector &gradient,
                                  const double      step_length) const;

  double
  choose_step_length(const TimeVector  &control,
                     const TimeVector  &gradient,
                     const Diagnostics &diagnostics) const;

  void
  output_time_sequence(const unsigned int iteration,
                       const TimeVector  &state,
                       const TimeVector  &adjoint,
                       const TimeVector  &control) const;

  void
  write_csv_header(std::ofstream &csv) const;

  void
  write_csv_line(std::ofstream     &csv,
                 const unsigned int iteration,
                 const Diagnostics &diagnostics) const;

  static std::string
  default_variable_names();

  Triangulation<dim> triangulation;

  std::unique_ptr<FESystem<dim>> fe;
  DoFHandler<dim>                dof_handler;

  AffineConstraints<double> constraints;

  BlockSparsityPattern      sparsity_pattern;
  BlockSparseMatrix<double> mass_matrix;
  BlockSparseMatrix<double> stiffness_matrix;
  BlockSparseMatrix<double> coupling_matrix;
  BlockSparseMatrix<double> tracking_matrix;

  SparseMatrix<double> state_time_step_matrix;
  SparseMatrix<double> adjoint_time_step_matrix;

  SparseDirectUMFPACK state_time_step_inverse;
  SparseDirectUMFPACK adjoint_time_step_inverse;
  SparseDirectUMFPACK control_mass_inverse;

  std::vector<types::global_dof_index> dofs_per_block;

  TimeVector     lower_bound;
  TimeVector     upper_bound;
  TimeVector     target_state;
  Vector<double> terminal_target;

  unsigned int state_degree;
  unsigned int control_degree;
  bool         continuous_control;
  unsigned int global_refinements;

  double       final_time;
  unsigned int n_time_steps;
  double       time_step;

  double       regularization;
  double       terminal_tracking_weight;
  unsigned int max_iterations;
  double       gradient_tolerance;
  double       projected_gradient_tolerance;
  double       fixed_step_length;
  std::string  output_name;
  std::string  grid_generator_function;
  std::string  grid_generator_arguments;

  OptimizationTools::ArmijoParameters<double> armijo_parameters;

  Functions::ParsedFunction<dim> diffusion_coefficient;
  Functions::ParsedFunction<dim> reaction_coefficient;
  Functions::ParsedFunction<dim> forcing_term;
  Functions::ParsedFunction<dim> initial_state;
  Functions::ParsedFunction<dim> initial_control;
  Functions::ParsedFunction<dim> desired_state;
  Functions::ParsedFunction<dim> terminal_target_function;
  Functions::ParsedFunction<dim> lower_bound_function;
  Functions::ParsedFunction<dim> upper_bound_function;
  Functions::ParsedFunction<dim> state_boundary_values;
  Functions::ParsedFunction<dim> adjoint_boundary_values;
};

#endif

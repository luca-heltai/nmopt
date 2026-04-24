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

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "kkt.h"

using namespace dealii;


template <int dim>
class BoxConstrainedKKT : public ParameterAcceptor
{
public:
  BoxConstrainedKKT();

  void
  run();

protected:
  void
  declare_parameters(ParameterHandler &prm) override;

  void
  parse_parameters(ParameterHandler &prm) override;

private:
  void
  initialize_control_mass_matrix();

  std::string
  output_iteration(const unsigned int      iteration,
                   const BlockVector<double> &solution,
                   const std::vector<bool> &lower_active,
                   const std::vector<bool> &upper_active) const;

  void
  output_results(const BlockVector<double> &solution,
                 const std::vector<bool>   &lower_active,
                 const std::vector<bool>   &upper_active) const;

  types::global_dof_index
  control_global_offset() const;

  void
  interpolate_control_function(const Function<dim> &function,
                               BlockVector<double> &vector) const;

  Vector<double>
  compute_stationarity(const BlockVector<double> &solution);

  AffineConstraints<double>
  make_control_constraints(const std::vector<bool> &lower_active,
                           const std::vector<bool> &upper_active) const;

  void
  solve_active_set_system(const std::vector<bool> &lower_active,
                          const std::vector<bool> &upper_active,
                          BlockVector<double>     &solution) const;

  KKT<dim>                       kkt_problem;
  Functions::ParsedFunction<dim> lower_bound_function;
  Functions::ParsedFunction<dim> upper_bound_function;

  unsigned int max_iterations;
  double       active_set_parameter;
  double       active_set_tolerance;
  std::string  output_name;

  BlockVector<double> lower_bound;
  BlockVector<double> upper_bound;

  SparseMatrix<double> control_mass_matrix;
  SparseDirectUMFPACK  control_mass_inverse;
};



template <int dim>
BoxConstrainedKKT<dim>::BoxConstrainedKKT()
  : ParameterAcceptor("Box-constrained KKT problem")
  , kkt_problem("KKT problem")
  , lower_bound_function(1)
  , upper_bound_function(1)
  , max_iterations(25)
  , active_set_parameter(1.0)
  , active_set_tolerance(1e-10)
  , output_name("kkt_box_constraints")
{}



template <int dim>
void
BoxConstrainedKKT<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.add_parameter("Max active-set iterations", max_iterations);
  prm.add_parameter("Active-set parameter", active_set_parameter);
  prm.add_parameter("Active-set tolerance", active_set_tolerance);
  prm.add_parameter("Output name", output_name);

  prm.enter_subsection("Lower bound");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "-0.25");
  prm.leave_subsection();

  prm.enter_subsection("Upper bound");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.set("Function expression", "0.25");
  prm.leave_subsection();
}



template <int dim>
void
BoxConstrainedKKT<dim>::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Lower bound");
  lower_bound_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Upper bound");
  upper_bound_function.parse_parameters(prm);
  prm.leave_subsection();
}



template <int dim>
void
BoxConstrainedKKT<dim>::interpolate_control_function(
  const Function<dim> &function,
  BlockVector<double> &vector) const
{
  vector.reinit(kkt_problem.get_solution());
  vector = 0;

  const VectorFunctionFromScalarFunctionObject<dim> control_function(
    [&function](const Point<dim> &point) { return function.value(point); },
    KKT<dim>::control_block,
    KKT<dim>::n_blocks);

  VectorTools::interpolate(kkt_problem.get_dof_handler(),
                           control_function,
                           vector);
  kkt_problem.get_constraints().distribute(vector);
}



template <int dim>
void
BoxConstrainedKKT<dim>::initialize_control_mass_matrix()
{
  const auto &control_block =
    kkt_problem.get_system_block(KKT<dim>::control_block,
                                 KKT<dim>::control_block);
  const double regularization = kkt_problem.get_regularization();

  AssertThrow(regularization > 0.0,
              ExcMessage("Regularization must be strictly positive."));

  DynamicSparsityPattern dsp(control_block.m(), control_block.n());
  for (unsigned int row = 0; row < control_block.m(); ++row)
    for (auto entry = control_block.begin(row); entry != control_block.end(row);
         ++entry)
      dsp.add(row, entry->column());

  SparsityPattern sparsity;
  sparsity.copy_from(dsp);
  control_mass_matrix.reinit(sparsity);

  for (unsigned int row = 0; row < control_block.m(); ++row)
    for (auto entry = control_block.begin(row); entry != control_block.end(row);
         ++entry)
      control_mass_matrix.set(row,
                              entry->column(),
                              entry->value() / regularization);

  control_mass_inverse.initialize(control_mass_matrix);
}



template <int dim>
std::string
BoxConstrainedKKT<dim>::output_iteration(
  const unsigned int        iteration,
  const BlockVector<double> &solution,
  const std::vector<bool>   &lower_active,
  const std::vector<bool>   &upper_active) const
{
  BlockVector<double> lower_active_output(solution);
  BlockVector<double> upper_active_output(solution);
  BlockVector<double> inactive_output(solution);

  lower_active_output = 0;
  upper_active_output = 0;
  inactive_output     = 0;

  for (unsigned int i = 0; i < lower_active.size(); ++i)
    if (lower_active[i])
      lower_active_output.block(KKT<dim>::control_block)[i] = 1.0;
    else if (upper_active[i])
      upper_active_output.block(KKT<dim>::control_block)[i] = 1.0;
    else
      inactive_output.block(KKT<dim>::control_block)[i] = 1.0;

  std::vector<std::string> solution_names = {"state", "adjoint", "control"};
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(KKT<dim>::n_blocks,
                   DataComponentInterpretation::component_is_scalar);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(kkt_problem.get_dof_handler());
  data_out.add_data_vector(solution,
                           solution_names,
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(lower_bound,
                           {"lower_bound_state",
                            "lower_bound_adjoint",
                            "lower_bound"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(upper_bound,
                           {"upper_bound_state",
                            "upper_bound_adjoint",
                            "upper_bound"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(lower_active_output,
                           {"lower_active_state",
                            "lower_active_adjoint",
                            "lower_active"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(upper_active_output,
                           {"upper_active_state",
                            "upper_active_adjoint",
                            "upper_active"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(inactive_output,
                           {"inactive_state", "inactive_adjoint", "inactive"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.build_patches();

  const std::string filename =
    output_name + "-" + Utilities::int_to_string(iteration, 4) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);

  return filename;
}



template <int dim>
void
BoxConstrainedKKT<dim>::output_results(
  const BlockVector<double> &solution,
  const std::vector<bool>   &lower_active,
  const std::vector<bool>   &upper_active) const
{
  BlockVector<double> lower_active_output(solution);
  BlockVector<double> upper_active_output(solution);
  BlockVector<double> inactive_output(solution);

  lower_active_output = 0;
  upper_active_output = 0;
  inactive_output     = 0;

  for (unsigned int i = 0; i < lower_active.size(); ++i)
    if (lower_active[i])
      lower_active_output.block(KKT<dim>::control_block)[i] = 1.0;
    else if (upper_active[i])
      upper_active_output.block(KKT<dim>::control_block)[i] = 1.0;
    else
      inactive_output.block(KKT<dim>::control_block)[i] = 1.0;

  std::vector<std::string> solution_names = {"state", "adjoint", "control"};
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(KKT<dim>::n_blocks,
                   DataComponentInterpretation::component_is_scalar);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(kkt_problem.get_dof_handler());
  data_out.add_data_vector(solution,
                           solution_names,
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(lower_bound,
                           {"lower_bound_state",
                            "lower_bound_adjoint",
                            "lower_bound"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(upper_bound,
                           {"upper_bound_state",
                            "upper_bound_adjoint",
                            "upper_bound"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(lower_active_output,
                           {"lower_active_state",
                            "lower_active_adjoint",
                            "lower_active"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(upper_active_output,
                           {"upper_active_state",
                            "upper_active_adjoint",
                            "upper_active"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.add_data_vector(inactive_output,
                           {"inactive_state", "inactive_adjoint", "inactive"},
                           DataOut<dim>::type_dof_data,
                           interpretation);
  data_out.build_patches();

  std::ofstream output(output_name + ".vtu");
  data_out.write_vtu(output);
}



template <int dim>
types::global_dof_index
BoxConstrainedKKT<dim>::control_global_offset() const
{
  const auto &dofs_per_block = kkt_problem.get_dofs_per_block();

  return dofs_per_block[KKT<dim>::state_block] +
         dofs_per_block[KKT<dim>::adjoint_block];
}



template <int dim>
Vector<double>
BoxConstrainedKKT<dim>::compute_stationarity(
  const BlockVector<double> &solution)
{
  BlockVector<double> residual(solution);
  kkt_problem.vmult_system(residual, solution);
  residual -= kkt_problem.get_system_rhs();

  Vector<double> dual_stationarity(solution.block(KKT<dim>::control_block).size());
  dual_stationarity = residual.block(KKT<dim>::control_block);

  Vector<double> stationarity(solution.block(KKT<dim>::control_block).size());
  control_mass_inverse.vmult(stationarity, dual_stationarity);

  return stationarity;
}



template <int dim>
AffineConstraints<double>
BoxConstrainedKKT<dim>::make_control_constraints(
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
          lower_active[i] ? lower_bound.block(KKT<dim>::control_block)[i] :
                            upper_bound.block(KKT<dim>::control_block)[i]);
      }

  control_constraints.close();

  AffineConstraints<double> combined_constraints;
  combined_constraints.copy_from(kkt_problem.get_constraints());
  combined_constraints.merge(control_constraints);

  return combined_constraints;
}



template <int dim>
void
BoxConstrainedKKT<dim>::solve_active_set_system(
  const std::vector<bool> &lower_active,
  const std::vector<bool> &upper_active,
  BlockVector<double>     &solution) const
{
  const AffineConstraints<double> combined_constraints =
    make_control_constraints(lower_active, upper_active);

  BlockSparseMatrix<double> modified_matrix;
  kkt_problem.copy_system_matrix(modified_matrix);
  BlockVector<double> modified_rhs = kkt_problem.get_system_rhs();

  combined_constraints.condense(modified_matrix, modified_rhs);

  SparseDirectUMFPACK direct_solver;
  direct_solver.initialize(modified_matrix);
  direct_solver.vmult(solution, modified_rhs);
  combined_constraints.distribute(solution);
}



template <int dim>
void
BoxConstrainedKKT<dim>::run()
{
  kkt_problem.initialize();
  initialize_control_mass_matrix();

  interpolate_control_function(lower_bound_function, lower_bound);
  interpolate_control_function(upper_bound_function, upper_bound);

  const auto &control_lower = lower_bound.block(KKT<dim>::control_block);
  const auto &control_upper = upper_bound.block(KKT<dim>::control_block);
  for (unsigned int i = 0; i < control_lower.size(); ++i)
    AssertThrow(control_lower[i] <= control_upper[i],
                ExcMessage("Lower bound exceeds upper bound."));

  kkt_problem.solve();
  BlockVector<double> current_solution = kkt_problem.get_solution();

  std::vector<bool> lower_active(control_lower.size(), false);
  std::vector<bool> upper_active(control_lower.size(), false);
  std::vector<std::pair<double, std::string>> times_and_names;

  for (unsigned int iteration = 0; iteration < max_iterations; ++iteration)
    {
      const Vector<double> stationarity =
        compute_stationarity(current_solution);
      std::vector<bool> next_lower_active(control_lower.size(), false);
      std::vector<bool> next_upper_active(control_lower.size(), false);

      for (unsigned int i = 0; i < control_lower.size(); ++i)
        {
          const double control_value =
            current_solution.block(KKT<dim>::control_block)[i];
          const double lower_indicator =
            -stationarity[i] +
            active_set_parameter * (control_value - control_lower[i]);
          const double upper_indicator =
            -stationarity[i] +
            active_set_parameter * (control_value - control_upper[i]);

          if (lower_indicator < -active_set_tolerance)
            next_lower_active[i] = true;
          else if (upper_indicator > active_set_tolerance)
            next_upper_active[i] = true;
        }

      BlockVector<double> next_solution(current_solution);
      solve_active_set_system(next_lower_active,
                              next_upper_active,
                              next_solution);

      Vector<double> control_update =
        next_solution.block(KKT<dim>::control_block);
      control_update -= current_solution.block(KKT<dim>::control_block);

      const unsigned int n_lower =
        std::count(next_lower_active.begin(), next_lower_active.end(), true);
      const unsigned int n_upper =
        std::count(next_upper_active.begin(), next_upper_active.end(), true);

      std::cout << "PDAS iteration " << iteration << ": |A_-|=" << n_lower
                << ", |A_+|=" << n_upper
                << ", ||delta_u||=" << control_update.l2_norm() << std::endl;

      times_and_names.emplace_back(
        static_cast<double>(iteration),
        output_iteration(
          iteration, next_solution, next_lower_active, next_upper_active));

      if (next_lower_active == lower_active &&
          next_upper_active == upper_active &&
          control_update.l2_norm() < active_set_tolerance)
        {
          lower_active     = std::move(next_lower_active);
          upper_active     = std::move(next_upper_active);
          current_solution = next_solution;
          break;
        }

      lower_active     = std::move(next_lower_active);
      upper_active     = std::move(next_upper_active);
      current_solution = next_solution;
    }

  kkt_problem.set_solution(current_solution);
  output_results(current_solution, lower_active, upper_active);

  std::ofstream pvd_output(output_name + ".pvd");
  DataOutBase::write_pvd_record(pvd_output, times_and_names);
}


int
main(int argc, char **argv)
{
  try
    {
      BoxConstrainedKKT<DEAL_II_DIMENSION> kkt_problem;

      const std::string parameter_file =
        (argc > 1 ? argv[1] : "parameters/kkt_box_constraints.prm");

      ParameterAcceptor::initialize(parameter_file, "last_used_parameters.prm");
      kkt_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}

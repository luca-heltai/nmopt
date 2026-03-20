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
#include "optimization_tools.h"

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/linear_operator_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/base/data_out_base.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;


template <int dim>
class UnconstrainedPoissonOptimization : public ParameterAcceptor
{
public:
  UnconstrainedPoissonOptimization();

  void
  run();

protected:
  void
  declare_parameters(ParameterHandler &prm) override;

  void
  parse_parameters(ParameterHandler &prm) override;

private:
  void
  interpolate_configuration(const Function<dim> &function,
                            Vector<double>      &vector) const;

  void
  solve_state(const Vector<double> &control, Vector<double> &state) const;

  void
  solve_adjoint(const Vector<double> &state,
                const Vector<double> &target_state,
                Vector<double>       &adjoint) const;

  void
  output_results(const Vector<double> &optimized_control,
                 const Vector<double> &target_state) const;

  std::string
  output_iteration(const unsigned int   iteration,
                   const Vector<double> &control,
                   const Vector<double> &target_state) const;

  Laplacian<dim>               poisson_problem;
  Functions::ParsedFunction<dim> initial_control_function;
  Functions::ParsedFunction<dim> target_state_function;

  std::string                                            optimization_method;
  OptimizationTools::OptimizationParameters<double>      optimization_parameters;
  OptimizationTools::NLCGParameters<double>              nlcg_parameters;
  OptimizationTools::LBFGSParameters<double>             lbfgs_parameters;
  double                                                 regularization;
  std::string                                            output_name;
};



template <int dim>
UnconstrainedPoissonOptimization<dim>::UnconstrainedPoissonOptimization()
  : ParameterAcceptor("Unconstrained Poisson optimization")
  , poisson_problem("Poisson problem")
  , initial_control_function(1)
  , target_state_function(1)
  , optimization_method("nlcg")
  , regularization(1e-4)
  , output_name("optimized_poisson_state")
{
  optimization_parameters.store_iterates = false;
  nlcg_parameters.store_iterates         = false;
  lbfgs_parameters.store_iterates        = false;
}



template <int dim>
void
UnconstrainedPoissonOptimization<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.add_parameter("Optimization method", optimization_method);
  prm.add_parameter("Regularization", regularization);
  prm.add_parameter("Output name", output_name);

  prm.enter_subsection("Optimization parameters");
  optimization_parameters.add_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("NLCG parameters");
  nlcg_parameters.add_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("LBFGS parameters");
  lbfgs_parameters.add_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Initial control");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1, "0.0");
  prm.leave_subsection();

  prm.enter_subsection("Target state");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1, "sin(pi*x)");
  prm.leave_subsection();
}



template <int dim>
void
UnconstrainedPoissonOptimization<dim>::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Initial control");
  initial_control_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Target state");
  target_state_function.parse_parameters(prm);
  prm.leave_subsection();
}



template <int dim>
void
UnconstrainedPoissonOptimization<dim>::interpolate_configuration(
  const Function<dim> &function,
  Vector<double>      &vector) const
{
  vector.reinit(poisson_problem.get_solution());
  VectorTools::interpolate(poisson_problem.get_dof_handler(), function, vector);
  poisson_problem.get_constraints().distribute(vector);
}



template <int dim>
void
UnconstrainedPoissonOptimization<dim>::solve_state(
  const Vector<double> &control,
  Vector<double>       &state) const
{
  const auto M = linear_operator<Vector<double>>(poisson_problem.get_mass_matrix());
  const Vector<double> rhs = poisson_problem.get_system_rhs();
  const Vector<double> state_rhs = rhs + M * control;

  poisson_problem.solve(state_rhs, state);
}



template <int dim>
void
UnconstrainedPoissonOptimization<dim>::solve_adjoint(
  const Vector<double> &state,
  const Vector<double> &target_state,
  Vector<double>       &adjoint) const
{
  const auto M = linear_operator<Vector<double>>(poisson_problem.get_mass_matrix());
  const Vector<double> mismatch = state - target_state;
  const Vector<double> adjoint_rhs = M * mismatch;

  poisson_problem.solve(adjoint_rhs, adjoint);
}



template <int dim>
void
UnconstrainedPoissonOptimization<dim>::output_results(
  const Vector<double> &optimized_control,
  const Vector<double> &target_state) const
{
  Vector<double> state;
  Vector<double> adjoint;

  solve_state(optimized_control, state);
  solve_adjoint(state, target_state, adjoint);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(poisson_problem.get_dof_handler());
  data_out.add_data_vector(optimized_control, "optimized_control");
  data_out.add_data_vector(state, "state");
  data_out.add_data_vector(adjoint, "adjoint");
  data_out.add_data_vector(target_state, "target_state");
  data_out.add_data_vector(poisson_problem.get_system_rhs(), "poisson_rhs");
  data_out.build_patches();

  std::ofstream output(output_name + ".vtu");
  data_out.write_vtu(output);
}



template <int dim>
std::string
UnconstrainedPoissonOptimization<dim>::output_iteration(
  const unsigned int   iteration,
  const Vector<double> &control,
  const Vector<double> &target_state) const
{
  Vector<double> state;
  Vector<double> adjoint;

  solve_state(control, state);
  solve_adjoint(state, target_state, adjoint);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(poisson_problem.get_dof_handler());
  data_out.add_data_vector(control, "control");
  data_out.add_data_vector(state, "state");
  data_out.add_data_vector(adjoint, "adjoint");
  data_out.add_data_vector(target_state, "target_state");
  data_out.add_data_vector(poisson_problem.get_system_rhs(), "poisson_rhs");
  data_out.build_patches();

  const std::string filename =
    output_name + "-" + Utilities::int_to_string(iteration, 4) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);

  return filename;
}



template <int dim>
void
UnconstrainedPoissonOptimization<dim>::run()
{
  poisson_problem.initialize();

  Vector<double> initial_control;
  Vector<double> target_state;

  interpolate_configuration(initial_control_function, initial_control);
  interpolate_configuration(target_state_function, target_state);

  const auto M = linear_operator<Vector<double>>(poisson_problem.get_mass_matrix());

  const auto value = [&](const Vector<double> &control) {
    Vector<double> state;
    solve_state(control, state);

    const Vector<double> mismatch = state - target_state;
    const Vector<double> weighted_mismatch = M * mismatch;
    const Vector<double> weighted_control = M * control;

    return 0.5 * (mismatch * weighted_mismatch) +
           0.5 * regularization * (control * weighted_control);
  };

  const auto gradient = [&](const Vector<double> &control) {
    Vector<double> state;
    Vector<double> adjoint;

    solve_state(control, state);
    solve_adjoint(state, target_state, adjoint);

    const Vector<double> result = M * (regularization * control + adjoint);
    return result;
  };

  std::vector<std::pair<double, std::string>> times_and_names;
  const OptimizationTools::IterationCallback<Vector<double>> callback =
    [&](const unsigned int iteration,
        const Vector<double> &control,
        const double /*fx*/,
        const double /*gnorm*/) {
      const std::string filename =
        output_iteration(iteration, control, target_state);
      times_and_names.emplace_back(static_cast<double>(iteration), filename);
    };

  Vector<double> optimized_control;

  if (optimization_method == "gd")
    {
      const auto result =
        OptimizationTools::optimize_gd(value,
                                       gradient,
                                       initial_control,
                                       optimization_parameters,
                                       callback);
      optimized_control = result.x;
    }
  else if (optimization_method == "nlcg")
    {
      const auto result = OptimizationTools::optimize_nlcg(
        value, gradient, initial_control, nlcg_parameters, callback);
      optimized_control = result.x;
    }
  else if (optimization_method == "bfgs")
    {
      const auto result = OptimizationTools::optimize_bfgs(
        value, gradient, initial_control, lbfgs_parameters, callback);
      optimized_control = result.x;
    }
  else
    AssertThrow(false, ExcMessage("Unknown optimization method."));

  Vector<double>       optimized_state;
  Vector<double>       optimized_adjoint;
  solve_state(optimized_control, optimized_state);
  solve_adjoint(optimized_state, target_state, optimized_adjoint);

  const Vector<double> final_gradient = gradient(optimized_control);
  const Vector<double> final_mismatch = optimized_state - target_state;
  const Vector<double> weighted_final_mismatch = M * final_mismatch;
  const double         final_l2_error =
    std::sqrt(final_mismatch * weighted_final_mismatch);

  std::cout << "Final objective value: " << value(optimized_control) << std::endl;
  std::cout << "Final gradient norm:   " << final_gradient.l2_norm()
            << std::endl;
  std::cout << "Final L2 error:        " << final_l2_error << std::endl;

  std::ofstream pvd_output(output_name + ".pvd");
  DataOutBase::write_pvd_record(pvd_output, times_and_names);
  output_results(optimized_control, target_state);
}


int
main(int argc, char **argv)
{
  try
    {
      UnconstrainedPoissonOptimization<DEAL_II_DIMENSION> optimization_problem;

      const std::string parameter_file =
        (argc > 1 ?
           argv[1] :
           "unconstrained_poisson_optimization.prm");

      ParameterAcceptor::initialize(parameter_file);
      optimization_problem.run();
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

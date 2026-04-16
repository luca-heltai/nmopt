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

#include "laplacian.h"
#include "optimization_tools.h"

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/linear_operator_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <fstream>
#include <iostream>

using namespace dealii;


template <int dim>
class BoxConstrainedPoissonOptimization : public ParameterAcceptor
{
public:
  BoxConstrainedPoissonOptimization();

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
  project_onto_box(Vector<double> &vector) const;

  void
  solve_state(const Vector<double> &control, Vector<double> &state) const;

  void
  solve_adjoint(const Vector<double> &state,
                const Vector<double> &target_state,
                Vector<double>       &adjoint) const;

  double
  value(const Vector<double> &control, const Vector<double> &target_state) const;

  Vector<double>
  gradient(const Vector<double> &control, const Vector<double> &target_state) const;

  Vector<double>
  projected_gradient_residual(const Vector<double> &control,
                              const Vector<double> &gradient) const;

  std::string
  output_iteration(const unsigned int    iteration,
                   const Vector<double> &control,
                   const Vector<double> &target_state) const;

  void
  output_results(const Vector<double> &control,
                 const Vector<double> &target_state) const;

  Laplacian<dim>                 poisson_problem;
  Functions::ParsedFunction<dim> initial_control_function;
  Functions::ParsedFunction<dim> target_state_function;
  Functions::ParsedFunction<dim> lower_bound_function;
  Functions::ParsedFunction<dim> upper_bound_function;

  OptimizationTools::OptimizationParameters<double> optimization_parameters;
  double                                            regularization;
  std::string                                       output_name;

  Vector<double> lower_bound;
  Vector<double> upper_bound;
};



template <int dim>
BoxConstrainedPoissonOptimization<dim>::BoxConstrainedPoissonOptimization()
  : ParameterAcceptor("Box-constrained Poisson optimization")
  , poisson_problem("Poisson problem")
  , initial_control_function(1)
  , target_state_function(1)
  , lower_bound_function(1)
  , upper_bound_function(1)
  , regularization(1e-4)
  , output_name("laplacian_box_constraints")
{
  optimization_parameters.store_iterates = false;
}



template <int dim>
void
BoxConstrainedPoissonOptimization<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.add_parameter("Regularization", regularization);
  prm.add_parameter("Output name", output_name);

  prm.enter_subsection("Optimization parameters");
  optimization_parameters.add_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Initial control");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();

  prm.enter_subsection("Target state");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();

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
BoxConstrainedPoissonOptimization<dim>::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Initial control");
  initial_control_function.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Target state");
  target_state_function.parse_parameters(prm);
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
BoxConstrainedPoissonOptimization<dim>::interpolate_configuration(
  const Function<dim> &function,
  Vector<double>      &vector) const
{
  vector.reinit(poisson_problem.get_solution());
  VectorTools::interpolate(poisson_problem.get_dof_handler(), function, vector);
  poisson_problem.get_constraints().distribute(vector);
}



template <int dim>
void
BoxConstrainedPoissonOptimization<dim>::project_onto_box(Vector<double> &vector) const
{
  for (unsigned int i = 0; i < vector.size(); ++i)
    vector[i] = std::clamp(vector[i], lower_bound[i], upper_bound[i]);

  poisson_problem.get_constraints().distribute(vector);
}



template <int dim>
void
BoxConstrainedPoissonOptimization<dim>::solve_state(
  const Vector<double> &control,
  Vector<double>       &state) const
{
  const auto M = linear_operator<Vector<double>>(poisson_problem.get_mass_matrix());
  const Vector<double> rhs       = poisson_problem.get_system_rhs();
  const Vector<double> state_rhs = rhs + M * control;

  poisson_problem.solve(state_rhs, state);
}



template <int dim>
void
BoxConstrainedPoissonOptimization<dim>::solve_adjoint(
  const Vector<double> &state,
  const Vector<double> &target_state,
  Vector<double>       &adjoint) const
{
  const auto M = linear_operator<Vector<double>>(poisson_problem.get_mass_matrix());
  const Vector<double> mismatch    = state - target_state;
  const Vector<double> adjoint_rhs = M * mismatch;

  poisson_problem.solve(adjoint_rhs, adjoint);
}



template <int dim>
double
BoxConstrainedPoissonOptimization<dim>::value(
  const Vector<double> &control,
  const Vector<double> &target_state) const
{
  Vector<double> state;
  solve_state(control, state);

  const auto M = linear_operator<Vector<double>>(poisson_problem.get_mass_matrix());
  const Vector<double> mismatch          = state - target_state;
  const Vector<double> weighted_mismatch = M * mismatch;
  const Vector<double> weighted_control  = M * control;

  return 0.5 * (mismatch * weighted_mismatch) +
         0.5 * regularization * (control * weighted_control);
}



template <int dim>
Vector<double>
BoxConstrainedPoissonOptimization<dim>::gradient(
  const Vector<double> &control,
  const Vector<double> &target_state) const
{
  Vector<double> state;
  Vector<double> adjoint;

  solve_state(control, state);
  solve_adjoint(state, target_state, adjoint);

  const auto M = linear_operator<Vector<double>>(poisson_problem.get_mass_matrix());
  return M * (regularization * control + adjoint);
}



template <int dim>
Vector<double>
BoxConstrainedPoissonOptimization<dim>::projected_gradient_residual(
  const Vector<double> &control,
  const Vector<double> &gradient_vector) const
{
  Vector<double> projected = control;
  projected.add(-1.0, gradient_vector);
  project_onto_box(projected);

  projected -= control;
  projected *= -1.0;
  return projected;
}



template <int dim>
std::string
BoxConstrainedPoissonOptimization<dim>::output_iteration(
  const unsigned int    iteration,
  const Vector<double> &control,
  const Vector<double> &target_state) const
{
  Vector<double> state;
  Vector<double> adjoint;
  Vector<double> lower_active(control);
  Vector<double> upper_active(control);
  Vector<double> inactive(control);

  solve_state(control, state);
  solve_adjoint(state, target_state, adjoint);

  lower_active = 0;
  upper_active = 0;
  inactive     = 0;
  for (unsigned int i = 0; i < control.size(); ++i)
    {
      if (std::abs(control[i] - lower_bound[i]) < 1e-10)
        lower_active[i] = 1.0;
      else if (std::abs(control[i] - upper_bound[i]) < 1e-10)
        upper_active[i] = 1.0;
      else
        inactive[i] = 1.0;
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler(poisson_problem.get_dof_handler());
  data_out.add_data_vector(control, "control");
  data_out.add_data_vector(state, "state");
  data_out.add_data_vector(adjoint, "adjoint");
  data_out.add_data_vector(target_state, "target_state");
  data_out.add_data_vector(lower_bound, "lower_bound");
  data_out.add_data_vector(upper_bound, "upper_bound");
  data_out.add_data_vector(lower_active, "lower_active");
  data_out.add_data_vector(upper_active, "upper_active");
  data_out.add_data_vector(inactive, "inactive");
  data_out.build_patches();

  const std::string filename =
    output_name + "-" + Utilities::int_to_string(iteration, 4) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);

  return filename;
}



template <int dim>
void
BoxConstrainedPoissonOptimization<dim>::output_results(
  const Vector<double> &control,
  const Vector<double> &target_state) const
{
  Vector<double> state;
  Vector<double> adjoint;
  Vector<double> lower_active(control);
  Vector<double> upper_active(control);
  Vector<double> inactive(control);

  solve_state(control, state);
  solve_adjoint(state, target_state, adjoint);

  lower_active = 0;
  upper_active = 0;
  inactive     = 0;
  for (unsigned int i = 0; i < control.size(); ++i)
    {
      if (std::abs(control[i] - lower_bound[i]) < 1e-10)
        lower_active[i] = 1.0;
      else if (std::abs(control[i] - upper_bound[i]) < 1e-10)
        upper_active[i] = 1.0;
      else
        inactive[i] = 1.0;
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler(poisson_problem.get_dof_handler());
  data_out.add_data_vector(control, "optimized_control");
  data_out.add_data_vector(state, "state");
  data_out.add_data_vector(adjoint, "adjoint");
  data_out.add_data_vector(target_state, "target_state");
  data_out.add_data_vector(lower_bound, "lower_bound");
  data_out.add_data_vector(upper_bound, "upper_bound");
  data_out.add_data_vector(lower_active, "lower_active");
  data_out.add_data_vector(upper_active, "upper_active");
  data_out.add_data_vector(inactive, "inactive");
  data_out.build_patches();

  std::ofstream output(output_name + ".vtu");
  data_out.write_vtu(output);
}



template <int dim>
void
BoxConstrainedPoissonOptimization<dim>::run()
{
  poisson_problem.initialize();

  Vector<double> control;
  Vector<double> target_state;

  interpolate_configuration(initial_control_function, control);
  interpolate_configuration(target_state_function, target_state);
  interpolate_configuration(lower_bound_function, lower_bound);
  interpolate_configuration(upper_bound_function, upper_bound);

  for (unsigned int i = 0; i < control.size(); ++i)
    AssertThrow(lower_bound[i] <= upper_bound[i],
                ExcMessage("Lower bound exceeds upper bound."));

  project_onto_box(control);

  std::vector<std::pair<double, std::string>> times_and_names;

  for (unsigned int iteration = 0;
       iteration < optimization_parameters.max_iterations;
       ++iteration)
    {
      const double         fx          = value(control, target_state);
      const Vector<double> grad        = gradient(control, target_state);
      const Vector<double> projected_g = projected_gradient_residual(control, grad);
      const double         pg_norm     = projected_g.l2_norm();

      std::cout << "Projected GD iteration " << iteration << ": f=" << fx
                << ", ||pg||=" << pg_norm;

      if (optimization_parameters.store_iterates)
        times_and_names.emplace_back(
          static_cast<double>(iteration),
          output_iteration(iteration, control, target_state));

      if (pg_norm < optimization_parameters.gradient_tolerance)
        {
          std::cout << std::endl;
          break;
        }

      double         alpha      = optimization_parameters.armijo.alpha0;
      bool           accepted   = false;
      Vector<double> trial      = control;
      Vector<double> direction  = control;
      const double   current_fx = fx;

      for (unsigned int backtrack = 0;
           backtrack < optimization_parameters.armijo.max_backtracks;
           ++backtrack)
        {
          trial = control;
          trial.add(-alpha, grad);
          project_onto_box(trial);

          direction = trial;
          direction -= control;

          const double descent = grad * direction;
          if (direction.l2_norm() > 0 &&
              value(trial, target_state) <=
                current_fx + optimization_parameters.armijo.c1 * descent)
            {
              accepted = true;
              break;
            }

          alpha *= optimization_parameters.armijo.beta;
          if (alpha <= optimization_parameters.armijo.minimum_alpha)
            break;
        }

      if (!accepted)
        {
          trial = control;
          trial.add(-optimization_parameters.armijo.minimum_alpha, grad);
          project_onto_box(trial);
          alpha = optimization_parameters.armijo.minimum_alpha;
        }

      control = trial;
      std::cout << ", alpha=" << alpha << std::endl;
    }

  const Vector<double> final_projected_gradient =
    projected_gradient_residual(control, gradient(control, target_state));
  std::cout << "Final objective value: " << value(control, target_state)
            << std::endl;
  std::cout << "Final projected-gradient norm: "
            << final_projected_gradient.l2_norm() << std::endl;

  if (optimization_parameters.store_iterates)
    {
      std::ofstream pvd_output(output_name + ".pvd");
      DataOutBase::write_pvd_record(pvd_output, times_and_names);
    }

  output_results(control, target_state);
}


int
main(int argc, char **argv)
{
  try
    {
      BoxConstrainedPoissonOptimization<DEAL_II_DIMENSION> optimization_problem;

      const std::string parameter_file =
        (argc > 1 ? argv[1] : "parameters/laplacian_box_constraints.prm");

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

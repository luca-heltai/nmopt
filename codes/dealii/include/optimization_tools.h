// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// Copyright (C) 2026 by Luca Heltai
//
// This file is part of the deal.II Testbench for NMOPT Laboratories.
//
// ------------------------------------------------------------------------

#ifndef dealii_optimization_tools_h
#define dealii_optimization_tools_h


#include <deal.II/base/config.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/linear_operator_tools.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * Algorithms and utilities for optimization methods used in the Numerical
 * Methods for Optimal Control laboratories.
 */
namespace OptimizationTools
{
  /**
   * Choice of nonlinear conjugate-gradient update.
   */
  enum class NLCGBetaType
  {
    fletcher_reeves,
    polak_ribiere_plus
  };


  /**
   * Common parameters for Armijo backtracking.
   */
  template <typename Number>
  struct ArmijoParameters
  {
    Number       alpha0         = Number(1.);
    Number       c1             = Number(1.e-4);
    Number       beta           = Number(0.5);
    unsigned int max_backtracks = 30;
    Number       minimum_alpha  = Number(1.e-16);

    void
    add_parameters(ParameterHandler &prm)
    {
      prm.add_parameter("Initial step length", alpha0);
      prm.add_parameter("Armijo coefficient", c1);
      prm.add_parameter("Backtracking reduction", beta);
      prm.add_parameter("Max backtracks", max_backtracks);
      prm.add_parameter("Minimum step length", minimum_alpha);
    }
  };


  /**
   * Generic stopping and logging parameters for unconstrained optimization
   * algorithms.
   */
  template <typename Number>
  struct OptimizationParameters
  {
    unsigned int             max_iterations      = 200;
    Number                   gradient_tolerance  = Number(1.e-8);
    bool                     store_iterates      = true;
    bool                     log_iterations      = true;
    ArmijoParameters<Number> armijo;

    void
    add_parameters(ParameterHandler &prm)
    {
      prm.add_parameter("Max iterations", max_iterations);
      prm.add_parameter("Gradient tolerance", gradient_tolerance);
      prm.add_parameter("Store iterates", store_iterates);
      prm.add_parameter("Log iterations", log_iterations);

      prm.enter_subsection("Armijo");
      armijo.add_parameters(prm);
      prm.leave_subsection();
    }
  };


  /**
   * Parameters for nonlinear conjugate-gradient.
   */
  template <typename Number>
  struct NLCGParameters : public OptimizationParameters<Number>
  {
    NLCGBetaType beta_type     = NLCGBetaType::polak_ribiere_plus;
    unsigned int restart_every = 50;

    void
    add_parameters(ParameterHandler &prm)
    {
      OptimizationParameters<Number>::add_parameters(prm);
      prm.add_parameter("Beta type", beta_type);
      prm.add_parameter("Restart every", restart_every);
    }
  };


  /**
   * Parameters for limited-memory BFGS.
   */
  template <typename Number>
  struct LBFGSParameters : public OptimizationParameters<Number>
  {
    unsigned int history_size = 10;

    void
    add_parameters(ParameterHandler &prm)
    {
      OptimizationParameters<Number>::add_parameters(prm);
      prm.add_parameter("History size", history_size);
    }
  };


  /**
   * Parameters for the Cauchy-point trust-region method.
   */
  template <typename Number>
  struct TrustRegionParameters : public OptimizationParameters<Number>
  {
    Number delta0    = Number(1.);
    Number delta_max = Number(10.);
    Number eta       = Number(0.1);

    void
    add_parameters(ParameterHandler &prm)
    {
      OptimizationParameters<Number>::add_parameters(prm);
      prm.add_parameter("Initial trust-region radius", delta0);
      prm.add_parameter("Maximum trust-region radius", delta_max);
      prm.add_parameter("Acceptance threshold", eta);
    }
  };


  /**
   * Iteration history shared by the algorithms below.
   */
  template <typename VectorType>
  struct OptimizationResult
  {
    using value_type = typename VectorType::value_type;

    VectorType              x;
    unsigned int            iterations = 0;
    std::vector<value_type> function_values;
    std::vector<value_type> gradient_norms;
    std::vector<value_type> step_lengths;
    std::vector<value_type> trust_region_radii;
    std::vector<VectorType> iterates;
  };


  namespace internal
  {
    template <typename VectorType>
    LinearOperator<VectorType>
    outer_product_operator(const VectorType &left, const VectorType &right)
    {
      const auto left_ptr  = std::make_shared<VectorType>(left);
      const auto right_ptr = std::make_shared<VectorType>(right);

      auto op = dealii::identity_operator<VectorType>(
        [left_ptr](VectorType &v, const bool omit_zeroing_entries) {
          v.reinit(*left_ptr, omit_zeroing_entries);
        });

      op.vmult = [left_ptr, right_ptr](VectorType &dst, const VectorType &src) {
        dst = *left_ptr;
        dst *= ((*right_ptr) * src);
      };

      op.vmult_add = [left_ptr, right_ptr](VectorType &dst,
                                           const VectorType &src) {
        dst.add((*right_ptr) * src, *left_ptr);
      };

      op.Tvmult = [left_ptr, right_ptr](VectorType &dst, const VectorType &src) {
        dst = *right_ptr;
        dst *= ((*left_ptr) * src);
      };

      op.Tvmult_add = [left_ptr, right_ptr](VectorType &dst,
                                            const VectorType &src) {
        dst.add((*left_ptr) * src, *right_ptr);
      };

      return op;
    }


    template <typename VectorType>
    void
    maybe_store_iterate(OptimizationResult<VectorType> &result,
                        const VectorType               &x,
                        const bool                      store_iterates)
    {
      if (store_iterates)
        result.iterates.push_back(x);
    }


    template <typename Number>
    void
    log_iteration(const bool         enabled,
                  const std::string &method,
                  const unsigned int iteration,
                  const Number       fx,
                  const Number       gnorm)
    {
      if (!enabled)
        return;

      std::cout << method << " iteration " << iteration << ": f=" << fx
                << ", ||g||=" << gnorm;
    }


    template <typename Number>
    void
    log_scalar(const bool enabled, const std::string &label, const Number value)
    {
      if (enabled)
        std::cout << ", " << label << "=" << value;
    }


    inline void
    log_endline(const bool enabled)
    {
      if (enabled)
        std::cout << std::endl;
    }
  } // namespace internal


  /**
   * Armijo backtracking line search for a differentiable objective.
   */
  template <typename VectorType, typename ValueFunction, typename GradientFunction>
  typename VectorType::value_type
  armijo_backtracking(
    const ValueFunction                                          &value,
    const GradientFunction                                       &gradient,
    const VectorType                                             &x,
    const VectorType                                             &p,
    const ArmijoParameters<typename VectorType::value_type> &parameters = {})
  {
    using Number = typename VectorType::value_type;

    const VectorType g       = gradient(x);
    const Number     fx      = value(x);
    const Number     g_dot_p = g * p;

    Number alpha = parameters.alpha0;

    for (unsigned int backtrack = 0; backtrack < parameters.max_backtracks;
         ++backtrack)
      {
        const VectorType trial = x + alpha * p;

        if (value(trial) <= fx + parameters.c1 * alpha * g_dot_p)
          return alpha;

        alpha *= parameters.beta;
        if (alpha <= parameters.minimum_alpha)
          break;
      }

    return std::max(alpha, parameters.minimum_alpha);
  }


  /**
   * Gradient descent with Armijo backtracking.
   */
  template <typename VectorType, typename ValueFunction, typename GradientFunction>
  OptimizationResult<VectorType>
  optimize_gd(
    const ValueFunction                                             &value,
    const GradientFunction                                          &gradient,
    const VectorType                                                &x0,
    const OptimizationParameters<typename VectorType::value_type> &parameters =
      {})
  {
    using Number = typename VectorType::value_type;

    OptimizationResult<VectorType> result;
    VectorType                     x = x0;

    internal::maybe_store_iterate(result, x, parameters.store_iterates);

    for (unsigned int k = 0; k < parameters.max_iterations; ++k)
      {
        const VectorType g     = gradient(x);
        const Number     gnorm = g.l2_norm();

        const Number fx = value(x);
        result.function_values.push_back(fx);
        result.gradient_norms.push_back(gnorm);
        result.iterations = k;

        internal::log_iteration(
          parameters.log_iterations, "GD", k, fx, gnorm);

        if (gnorm < parameters.gradient_tolerance)
          {
            internal::log_endline(parameters.log_iterations);
            break;
          }

        const VectorType p     = Number(-1.) * g;
        const Number     alpha =
          armijo_backtracking(value, gradient, x, p, parameters.armijo);

        x = x + alpha * p;

        result.step_lengths.push_back(alpha);
        internal::log_scalar(parameters.log_iterations, "alpha", alpha);
        internal::log_endline(parameters.log_iterations);
        internal::maybe_store_iterate(result, x, parameters.store_iterates);
      }

    result.x = x;
    return result;
  }


  /**
   * Nonlinear conjugate-gradient with Armijo backtracking.
   */
  template <typename VectorType, typename ValueFunction, typename GradientFunction>
  OptimizationResult<VectorType>
  optimize_nlcg(const ValueFunction                                  &value,
                const GradientFunction                               &gradient,
                const VectorType                                     &x0,
                const NLCGParameters<typename VectorType::value_type> &parameters =
                  {})
  {
    using Number = typename VectorType::value_type;

    OptimizationResult<VectorType> result;
    VectorType                     x = x0;
    VectorType                     g = gradient(x);
    VectorType                     p = Number(-1.) * g;

    internal::maybe_store_iterate(result, x, parameters.store_iterates);

    for (unsigned int k = 0; k < parameters.max_iterations; ++k)
      {
        const Number gnorm = g.l2_norm();

        const Number fx = value(x);
        result.function_values.push_back(fx);
        result.gradient_norms.push_back(gnorm);
        result.iterations = k;

        internal::log_iteration(
          parameters.log_iterations, "NLCG", k, fx, gnorm);

        if (gnorm < parameters.gradient_tolerance)
          {
            internal::log_endline(parameters.log_iterations);
            break;
          }

        const Number alpha =
          armijo_backtracking(value, gradient, x, p, parameters.armijo);

        const VectorType x_new = x + alpha * p;
        const VectorType g_new = gradient(x_new);

        Number beta = Number();

        const Number gg = std::max(g * g, std::numeric_limits<Number>::epsilon());

        switch (parameters.beta_type)
          {
            case NLCGBetaType::fletcher_reeves:
              beta = (g_new * g_new) / gg;
              break;

            case NLCGBetaType::polak_ribiere_plus:
              {
                const VectorType y = g_new - g;
                beta               = std::max(Number(), (g_new * y) / gg);
                break;
              }
          }

        if ((k + 1) % parameters.restart_every == 0)
          beta = Number();

        p = Number(-1.) * g_new + beta * p;

        if ((p * g_new) >= Number())
          p = Number(-1.) * g_new;

        x = x_new;
        g = g_new;

        result.step_lengths.push_back(alpha);
        internal::log_scalar(parameters.log_iterations, "alpha", alpha);
        internal::log_scalar(parameters.log_iterations, "beta", beta);
        internal::log_endline(parameters.log_iterations);
        internal::maybe_store_iterate(result, x, parameters.store_iterates);
      }

    result.x = x;
    return result;
  }


  /**
   * Limited-memory BFGS with Armijo backtracking.
   */
  template <typename VectorType, typename ValueFunction, typename GradientFunction>
  OptimizationResult<VectorType>
  optimize_bfgs(
    const ValueFunction                                             &value,
    const GradientFunction                                          &gradient,
    const VectorType                                                &x0,
    const LBFGSParameters<typename VectorType::value_type>         &parameters =
      {})
  {
    using Number = typename VectorType::value_type;

    OptimizationResult<VectorType> result;
    VectorType                     x = x0;
    std::vector<VectorType>        s_history;
    std::vector<VectorType>        y_history;
    std::vector<Number>            rho_history;

    internal::maybe_store_iterate(result, x, parameters.store_iterates);

    for (unsigned int k = 0; k < parameters.max_iterations; ++k)
      {
        const VectorType g     = gradient(x);
        const Number     gnorm = g.l2_norm();

        const Number fx = value(x);
        result.function_values.push_back(fx);
        result.gradient_norms.push_back(gnorm);
        result.iterations = k;

        internal::log_iteration(
          parameters.log_iterations, "LBFGS", k, fx, gnorm);

        if (gnorm < parameters.gradient_tolerance)
          {
            internal::log_endline(parameters.log_iterations);
            break;
          }

        VectorType           q = g;
        const unsigned int   m = s_history.size();
        std::vector<Number>  alphas(m);

        for (int i = static_cast<int>(m) - 1; i >= 0; --i)
          {
            alphas[i] = rho_history[i] * (s_history[i] * q);
            q -= alphas[i] * y_history[i];
          }

        Number gamma = Number(1.);
        if (m > 0)
          {
            const Number yy = y_history.back() * y_history.back();
            if (yy > Number(0.))
              gamma = (s_history.back() * y_history.back()) / yy;
          }

        VectorType r = gamma * q;

        for (unsigned int i = 0; i < m; ++i)
          {
            const Number beta = rho_history[i] * (y_history[i] * r);
            r += (alphas[i] - beta) * s_history[i];
          }

        VectorType p = Number(-1.) * r;

        if ((p * g) >= Number())
          {
            p = Number(-1.) * g;
            s_history.clear();
            y_history.clear();
            rho_history.clear();
          }

        const Number alpha =
          armijo_backtracking(value, gradient, x, p, parameters.armijo);

        const VectorType s     = alpha * p;
        const VectorType x_new = x + s;
        const VectorType y     = gradient(x_new) - g;
        const Number     ys    = y * s;

        if (ys > Number(1.e-12))
          {
            if (s_history.size() == parameters.history_size)
              {
                s_history.erase(s_history.begin());
                y_history.erase(y_history.begin());
                rho_history.erase(rho_history.begin());
              }

            s_history.push_back(s);
            y_history.push_back(y);
            rho_history.push_back(Number(1.) / ys);
          }
        else
          {
            s_history.clear();
            y_history.clear();
            rho_history.clear();
          }

        x = x_new;

        result.step_lengths.push_back(alpha);
        internal::log_scalar(parameters.log_iterations, "alpha", alpha);
        internal::log_scalar(
          parameters.log_iterations, "history", s_history.size());
        internal::log_endline(parameters.log_iterations);
        internal::maybe_store_iterate(result, x, parameters.store_iterates);
      }

    result.x = x;
    return result;
  }


  /**
   * Trust-region method based on the Cauchy step.
   *
   * The Hessian model is provided as a LinearOperator-valued callable.
   */
  template <typename VectorType,
            typename ValueFunction,
            typename GradientFunction,
            typename HessianFunction>
  OptimizationResult<VectorType>
  optimize_trust_region_cauchy(
    const ValueFunction                                        &value,
    const GradientFunction                                     &gradient,
    const HessianFunction                                      &hessian,
    const VectorType                                           &x0,
    const TrustRegionParameters<typename VectorType::value_type> &parameters =
      {})
  {
    using Number = typename VectorType::value_type;

    OptimizationResult<VectorType> result;
    VectorType                     x     = x0;
    Number                         delta = parameters.delta0;

    internal::maybe_store_iterate(result, x, parameters.store_iterates);

    for (unsigned int k = 0; k < parameters.max_iterations; ++k)
      {
        const VectorType g     = gradient(x);
        const Number     gnorm = g.l2_norm();
        const auto       B     = hessian(x);

        const Number fx = value(x);
        result.function_values.push_back(fx);
        result.gradient_norms.push_back(gnorm);
        result.trust_region_radii.push_back(delta);
        result.iterations = k;

        internal::log_iteration(
          parameters.log_iterations, "TR", k, fx, gnorm);
        internal::log_scalar(parameters.log_iterations, "delta", delta);

        if (gnorm < parameters.gradient_tolerance)
          {
            internal::log_endline(parameters.log_iterations);
            break;
          }

        const VectorType Bg  = B * g;
        const Number     gBg = g * Bg;

        Number tau = Number(1.);
        if (gBg > Number())
          tau = std::min((gnorm * gnorm * gnorm) / (delta * gBg), Number(1.));

        const VectorType p = (-tau * delta / gnorm) * g;

        const Number ared = value(x) - value(x + p);
        const Number pred =
          -(g * p + Number(0.5) * (p * (B * p)));
        const Number rho =
          ared / std::max(pred, std::numeric_limits<Number>::epsilon());

        if (rho < Number(0.25))
          delta *= Number(0.25);
        else if (rho > Number(0.75) &&
                 std::abs(p.l2_norm() - delta) < Number(1.e-12))
          delta = std::min(Number(2.) * delta, parameters.delta_max);

        if (rho > parameters.eta)
          {
            x = x + p;
            result.step_lengths.push_back(p.l2_norm());
            internal::maybe_store_iterate(result, x, parameters.store_iterates);
          }

        internal::log_scalar(parameters.log_iterations, "rho", rho);
        internal::log_endline(parameters.log_iterations);
      }

    result.x = x;
    return result;
  }
} // namespace OptimizationTools

DEAL_II_NAMESPACE_CLOSE

#endif

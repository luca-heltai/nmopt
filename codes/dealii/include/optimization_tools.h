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

#include <deal.II/lac/linear_operator_tools.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
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
    ArmijoParameters<Number> armijo;
  };


  /**
   * Parameters for nonlinear conjugate-gradient.
   */
  template <typename Number>
  struct NLCGParameters : public OptimizationParameters<Number>
  {
    NLCGBetaType beta_type     = NLCGBetaType::polak_ribiere_plus;
    unsigned int restart_every = 50;
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

        result.function_values.push_back(value(x));
        result.gradient_norms.push_back(gnorm);
        result.iterations = k;

        if (gnorm < parameters.gradient_tolerance)
          break;

        const VectorType p     = Number(-1.) * g;
        const Number     alpha =
          armijo_backtracking(value, gradient, x, p, parameters.armijo);

        x = x + alpha * p;

        result.step_lengths.push_back(alpha);
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

        result.function_values.push_back(value(x));
        result.gradient_norms.push_back(gnorm);
        result.iterations = k;

        if (gnorm < parameters.gradient_tolerance)
          break;

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
        internal::maybe_store_iterate(result, x, parameters.store_iterates);
      }

    result.x = x;
    return result;
  }


  /**
   * BFGS with Armijo backtracking and inverse-Hessian approximation stored as
   * a LinearOperator.
   */
  template <typename VectorType, typename ValueFunction, typename GradientFunction>
  OptimizationResult<VectorType>
  optimize_bfgs(
    const ValueFunction                                             &value,
    const GradientFunction                                          &gradient,
    const VectorType                                                &x0,
    const OptimizationParameters<typename VectorType::value_type> &parameters =
      {})
  {
    using Number = typename VectorType::value_type;

    OptimizationResult<VectorType> result;
    VectorType                     x = x0;

    auto H = dealii::identity_operator<VectorType>(
      [&x0](VectorType &v, const bool omit_zeroing_entries) {
        v.reinit(x0, omit_zeroing_entries);
      });
    const auto I = dealii::identity_operator<VectorType>(
      [&x0](VectorType &v, const bool omit_zeroing_entries) {
        v.reinit(x0, omit_zeroing_entries);
      });

    internal::maybe_store_iterate(result, x, parameters.store_iterates);

    for (unsigned int k = 0; k < parameters.max_iterations; ++k)
      {
        const VectorType g     = gradient(x);
        const Number     gnorm = g.l2_norm();

        result.function_values.push_back(value(x));
        result.gradient_norms.push_back(gnorm);
        result.iterations = k;

        if (gnorm < parameters.gradient_tolerance)
          break;

        VectorType p = Number(-1.) * (H * g);

        if ((p * g) >= Number())
          {
            p = Number(-1.) * g;
            H = I;
          }

        const Number alpha =
          armijo_backtracking(value, gradient, x, p, parameters.armijo);

        const VectorType s     = alpha * p;
        const VectorType x_new = x + s;
        const VectorType y     = gradient(x_new) - g;
        const Number     ys    = y * s;

        if (ys > Number(1.e-12))
          {
            const Number rho   = Number(1.) / ys;
            const auto   sy    = internal::outer_product_operator(s, y);
            const auto   ys_op = internal::outer_product_operator(y, s);
            const auto   ss    = internal::outer_product_operator(s, s);

            H = (I - rho * sy) * H * (I - rho * ys_op) + rho * ss;
          }
        else
          H = I;

        x = x_new;

        result.step_lengths.push_back(alpha);
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

        result.function_values.push_back(value(x));
        result.gradient_norms.push_back(gnorm);
        result.trust_region_radii.push_back(delta);
        result.iterations = k;

        if (gnorm < parameters.gradient_tolerance)
          break;

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
      }

    result.x = x;
    return result;
  }
} // namespace OptimizationTools

DEAL_II_NAMESPACE_CLOSE

#endif

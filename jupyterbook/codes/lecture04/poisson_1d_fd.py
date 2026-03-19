from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

CODE_ROOT = Path(__file__).resolve().parents[1]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from common.elliptic_control_1d import ReducedPoissonControl1D, reduced_gradient_history
from common.fd1d import UniformGrid1D, l2_inner, l2_norm


def finite_difference_check(problem: ReducedPoissonControl1D, control: np.ndarray) -> None:
    direction = np.sin(np.pi * problem.grid.x)
    gradient = problem.reduced_gradient(control)
    directional_derivative = l2_inner(gradient, direction, problem.grid)

    print("Finite-difference gradient check")
    for step in [1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4]:
        quotient = (problem.reduced_cost(control + step * direction) - problem.reduced_cost(control)) / step
        error = abs(quotient - directional_derivative)
        print(f"  h={step:7.1e}  FD quotient={quotient: .6e}  error={error: .6e}")


def reduced_gradient_descent(
    problem: ReducedPoissonControl1D,
    control0: np.ndarray,
    step_size: float,
    max_iter: int = 25,
) -> np.ndarray:
    print("\nReduced gradient descent")
    history = reduced_gradient_history(problem, control0, step_size=step_size, max_iter=max_iter)
    for item in history[:-1]:
        print(f"  iter={item.iteration:02d}  J(u)={item.value: .6e}  ||grad J||={item.gradient_norm: .6e}")
        if item.gradient_norm < 1.0e-8:
            return item.control
    return history[-1].control


def main() -> None:
    grid = UniformGrid1D(num_interior=80)
    desired_state = np.sin(np.pi * grid.x)
    alpha = 1.0e-2

    problem = ReducedPoissonControl1D(
        grid=grid,
        alpha=alpha,
        desired_state=desired_state,
    )

    control0 = np.zeros_like(grid.x)
    finite_difference_check(problem, control0)

    optimal_control = reduced_gradient_descent(
        problem=problem,
        control0=control0,
        step_size=5.0e-3,
        max_iter=20,
    )

    final_state = problem.state(optimal_control)
    final_adjoint = problem.adjoint(final_state)
    print("\nFinal summary")
    print(f"  ||u||_L2        = {l2_norm(optimal_control, grid):.6e}")
    print(f"  ||y - y_d||_L2  = {l2_norm(final_state - desired_state, grid):.6e}")
    print(f"  ||p + alpha u|| = {l2_norm(final_adjoint + alpha * optimal_control, grid):.6e}")


if __name__ == "__main__":
    main()

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fd1d import UniformGrid1D, l2_inner, l2_norm, solve_poisson_dirichlet


@dataclass
class ReducedPoissonControl1D:
    grid: UniformGrid1D
    alpha: float
    desired_state: np.ndarray
    forcing: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.forcing is None:
            self.forcing = np.zeros_like(self.grid.x)

    def state(self, control: np.ndarray) -> np.ndarray:
        return solve_poisson_dirichlet(control + self.forcing, self.grid)

    def adjoint(self, state: np.ndarray) -> np.ndarray:
        return solve_poisson_dirichlet(state - self.desired_state, self.grid)

    def reduced_cost(self, control: np.ndarray) -> float:
        state = self.state(control)
        tracking = 0.5 * l2_inner(state - self.desired_state, state - self.desired_state, self.grid)
        regularization = 0.5 * self.alpha * l2_inner(control, control, self.grid)
        return tracking + regularization

    def reduced_gradient(self, control: np.ndarray) -> np.ndarray:
        state = self.state(control)
        adjoint = self.adjoint(state)
        return adjoint + self.alpha * control

    def value_and_gradient(self, control: np.ndarray) -> tuple[float, np.ndarray]:
        state = self.state(control)
        adjoint = self.adjoint(state)
        value = (
            0.5 * l2_inner(state - self.desired_state, state - self.desired_state, self.grid)
            + 0.5 * self.alpha * l2_inner(control, control, self.grid)
        )
        gradient = adjoint + self.alpha * control
        return value, gradient


@dataclass
class ReducedIteration:
    iteration: int
    control: np.ndarray
    state: np.ndarray
    adjoint: np.ndarray
    value: float
    gradient_norm: float


def reduced_gradient_history(
    problem: ReducedPoissonControl1D,
    control0: np.ndarray,
    step_size: float,
    max_iter: int = 25,
) -> list[ReducedIteration]:
    control = np.array(control0, dtype=float, copy=True)
    history: list[ReducedIteration] = []

    for iteration in range(max_iter + 1):
        state = problem.state(control)
        adjoint = problem.adjoint(state)
        gradient = adjoint + problem.alpha * control
        history.append(
            ReducedIteration(
                iteration=iteration,
                control=control.copy(),
                state=state.copy(),
                adjoint=adjoint.copy(),
                value=problem.reduced_cost(control),
                gradient_norm=l2_norm(gradient, problem.grid),
            )
        )
        if iteration == max_iter:
            break
        control = control - step_size * gradient

    return history


def reduced_steepest_descent_history(
    problem: ReducedPoissonControl1D,
    control0: np.ndarray,
    max_iter: int = 40,
    tol: float = 1.0e-10,
) -> list[ReducedIteration]:
    """Reduced steepest descent with exact line search for the linear-quadratic model."""
    control = np.array(control0, dtype=float, copy=True)
    history: list[ReducedIteration] = []

    for iteration in range(max_iter + 1):
        state = problem.state(control)
        adjoint = problem.adjoint(state)
        gradient = adjoint + problem.alpha * control
        gradient_norm = l2_norm(gradient, problem.grid)
        history.append(
            ReducedIteration(
                iteration=iteration,
                control=control.copy(),
                state=state.copy(),
                adjoint=adjoint.copy(),
                value=problem.reduced_cost(control),
                gradient_norm=gradient_norm,
            )
        )
        if iteration == max_iter or gradient_norm <= tol:
            break

        descent = -gradient
        state_descent = problem.state(descent)
        numerator = l2_inner(gradient, gradient, problem.grid)
        denominator = l2_inner(state_descent, state_descent, problem.grid) + problem.alpha * numerator
        step_size = numerator / denominator
        control = control + step_size * descent

    return history

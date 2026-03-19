from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class GradientDescentResult:
    iterate: np.ndarray
    values: list[float]
    gradient_norms: list[float]


def armijo_backtracking(
    x: np.ndarray,
    value: float,
    gradient: np.ndarray,
    direction: np.ndarray,
    value_fn,
    step0: float = 1.0,
    beta: float = 0.5,
    c1: float = 1.0e-4,
    max_backtracks: int = 25,
) -> float:
    """Standard Armijo backtracking for unconstrained descent."""
    directional_derivative = float(np.dot(gradient, direction))
    step = step0
    for _ in range(max_backtracks):
        candidate = x + step * direction
        if value_fn(candidate) <= value + c1 * step * directional_derivative:
            return step
        step *= beta
    return step


def gradient_descent(
    x0: np.ndarray,
    value_and_gradient,
    step0: float = 1.0,
    beta: float = 0.5,
    c1: float = 1.0e-4,
    tol: float = 1.0e-8,
    max_iter: int = 200,
) -> GradientDescentResult:
    """Unconstrained gradient descent with Armijo backtracking."""
    x = np.array(x0, dtype=float, copy=True)
    values: list[float] = []
    gradient_norms: list[float] = []

    for _ in range(max_iter):
        value, gradient = value_and_gradient(x)
        gradient_norm = float(np.linalg.norm(gradient))
        values.append(float(value))
        gradient_norms.append(gradient_norm)
        if gradient_norm <= tol:
            break

        direction = -gradient
        step = armijo_backtracking(
            x=x,
            value=float(value),
            gradient=gradient,
            direction=direction,
            value_fn=lambda z: float(value_and_gradient(z)[0]),
            step0=step0,
            beta=beta,
            c1=c1,
        )
        x = x + step * direction

    return GradientDescentResult(iterate=x, values=values, gradient_norms=gradient_norms)

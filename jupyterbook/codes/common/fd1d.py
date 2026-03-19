from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class UniformGrid1D:
    """Uniform grid on [0, 1] with homogeneous Dirichlet boundaries."""

    num_interior: int
    a: float = 0.0
    b: float = 1.0

    @property
    def h(self) -> float:
        return (self.b - self.a) / (self.num_interior + 1)

    @property
    def x(self) -> np.ndarray:
        return np.linspace(self.a + self.h, self.b - self.h, self.num_interior)

    @property
    def x_with_boundary(self) -> np.ndarray:
        return np.linspace(self.a, self.b, self.num_interior + 2)


def poisson_matrix(grid: UniformGrid1D) -> np.ndarray:
    """Finite-difference matrix for -d^2/dx^2 on the interior nodes."""
    main = 2.0 * np.ones(grid.num_interior)
    off = -1.0 * np.ones(grid.num_interior - 1)
    return (np.diag(main) + np.diag(off, 1) + np.diag(off, -1)) / grid.h**2


def solve_poisson_dirichlet(rhs: np.ndarray, grid: UniformGrid1D) -> np.ndarray:
    """Solve -y'' = rhs on (0, 1) with y(0)=y(1)=0."""
    return np.linalg.solve(poisson_matrix(grid), rhs)


def l2_inner(u: np.ndarray, v: np.ndarray, grid: UniformGrid1D) -> float:
    return grid.h * float(np.dot(u, v))


def l2_norm(u: np.ndarray, grid: UniformGrid1D) -> float:
    return np.sqrt(l2_inner(u, u, grid))


def with_boundary_values(values: np.ndarray, left: float = 0.0, right: float = 0.0) -> np.ndarray:
    """Pad an interior vector with boundary values."""
    return np.concatenate(([left], values, [right]))

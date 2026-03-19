from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

CODE_ROOT = Path(__file__).resolve().parents[1]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from common.optim import gradient_descent


def quadratic_value_and_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
    hessian = np.array([[10.0, 0.0], [0.0, 1.0]])
    rhs = np.array([1.0, -2.0])
    value = 0.5 * float(x @ hessian @ x) - float(rhs @ x)
    gradient = hessian @ x - rhs
    return value, gradient


def main() -> None:
    x0 = np.array([2.0, 2.0])
    result = gradient_descent(x0, quadratic_value_and_gradient, step0=1.0, max_iter=100)

    print("Lecture 03 gradient-descent demo")
    print(f"iterations: {len(result.values) - 1}")
    print(f"final iterate: {result.iterate}")
    print(f"final value: {result.values[-1]:.6e}")
    print(f"final gradient norm: {result.gradient_norms[-1]:.6e}")


if __name__ == "__main__":
    main()

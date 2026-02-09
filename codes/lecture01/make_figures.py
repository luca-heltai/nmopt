#!/usr/bin/env python3

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "jupyterbook" / "slides" / "assets" / "lecture01"


def save(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def control_to_state_diagram() -> None:
    # Simple schematic drawn with matplotlib primitives to avoid extra deps.
    fig, ax = plt.subplots(figsize=(8.5, 2.0))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    def box(x: float, y: float, w: float, h: float, text: str) -> None:
        rect = plt.Rectangle((x, y), w, h, fill=False, linewidth=2)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=12)

    box(0.05, 0.35, 0.18, 0.30, "control\n$u$")
    box(0.30, 0.35, 0.18, 0.30, "PDE\nsolve")
    box(0.55, 0.35, 0.18, 0.30, "state\n$y$")
    box(0.80, 0.35, 0.18, 0.30, "obs.\n$\\mathcal{O}(y)$")

    def arrow(x0: float, y0: float, x1: float, y1: float) -> None:
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(arrowstyle="->", linewidth=2),
        )

    arrow(0.23, 0.50, 0.30, 0.50)
    arrow(0.48, 0.50, 0.55, 0.50)
    arrow(0.73, 0.50, 0.80, 0.50)

    ax.text(
        0.50,
        0.13,
        r"Reduced cost: $f(u) = J(S(u), u)$",
        ha="center",
        va="center",
        fontsize=12,
    )

    save(fig, "control_to_state.png")


def reduced_cost_1d() -> None:
    # A tiny LQ example with a 1D control u and linear constraint y = S u.
    alpha = 1e-2
    s = 2.0  # scalar "control-to-state" gain
    y_d = 1.0

    u = np.linspace(-2.0, 2.0, 400)
    y = s * u
    f = 0.5 * (y - y_d) ** 2 + 0.5 * alpha * u**2

    u_star = (s * y_d) / (s**2 + alpha)

    fig, ax = plt.subplots(figsize=(6.0, 3.5))
    ax.plot(u, f, linewidth=2)
    ax.axvline(u_star, linestyle="--", linewidth=1.5)
    ax.set_title("Reduced Cost in a 1D LQ Example")
    ax.set_xlabel("$u$")
    ax.set_ylabel("$f(u)$")
    ax.grid(True, alpha=0.25)
    ax.text(
        u_star + 0.05,
        float(np.min(f)) + 0.05,
        rf"$u^* \approx {u_star:.3f}$",
        fontsize=10,
    )

    save(fig, "reduced_cost_1d.png")


def main() -> None:
    control_to_state_diagram()
    reduced_cost_1d()
    print(f"Wrote figures to: {OUT_DIR}")


if __name__ == "__main__":
    main()


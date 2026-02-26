#!/usr/bin/env python3

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "jupyterbook" / "slides" / "assets"


def save(fig: plt.Figure, name: str, tight: bool = True) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    if tight:
        fig.savefig(path, dpi=200, bbox_inches="tight")
    else:
        fig.savefig(path, dpi=200)
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
        ax.text(x + w / 2, y + h / 2, text,
                ha="center", va="center", fontsize=12)

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

    save(fig, "01_control_to_state.png")


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

    save(fig, "01_reduced_cost_1d.png")


def minimization_geometry() -> None:
    # 2D constrained minimization example:
    # min 0.5 * u^T A u s.t. B u = g.
    A = np.array([[3.0, 1.0], [1.0, 2.0]], dtype=float)
    B = np.array([[2.0, -1.0]], dtype=float)
    g = np.array([-1], dtype=float)

    # Solve KKT system [A, -B^T; -B, 0] [u; lambda] = [0; -g].
    kkt = np.block([[A, -B.T], [-B, np.zeros((1, 1))]])
    rhs = np.concatenate([np.zeros(2), -g])
    sol = np.linalg.solve(kkt, rhs)
    u_star = sol[:2]

    # Build contour plot for f(u) and overlay the linear constraint.
    u1 = np.linspace(-1.2, 1.8, 400)
    u2 = np.linspace(-1.2, 1.8, 400)
    U1, U2 = np.meshgrid(u1, u2)
    F = 0.5 * (
        A[0, 0] * U1**2 + 2.0 * A[0, 1] * U1 * U2 + A[1, 1] * U2**2
    )

    f_star = 0.5 * u_star @ A @ u_star
    min_f = F.min()
    max_f = F.max()
    df = (f_star - min_f)/4.01

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    levels = np.linspace(min_f+0.01, max_f, 21)
    levels[-1] = f_star
    levels.sort()
    cmap = plt.get_cmap("Greens")
    colors = cmap(np.linspace(0.35, 0.95, len(levels)))
    ax.contour(U1, U2, F, levels=levels, colors=colors, linewidths=1.5)

    # Constraint: B u = g => u2 = (g - B[0]*u1) / B[1].
    line_u1 = np.array([u1.min(), u1.max()])
    line_u2 = (g[0] - B[0, 0]*line_u1) / B[0, 1]
    ax.plot(line_u1, line_u2, color="purple",
            linewidth=2.2, label=r"$\varphi(u)=Bu-g=0$")

    # Feasible minimizer and gradients at u*.
    ax.scatter(u_star[0], u_star[1], color="navy", s=50, zorder=5)
    ax.text(u_star[0] + 0.05, u_star[1] + 0.03, r"$u^\star$", color="navy")

    grad_f = A @ u_star
    grad_phi = B.ravel()
    scale = 0.28
    ax.arrow(
        u_star[0],
        u_star[1],
        scale * grad_f[0],
        scale * grad_f[1],
        width=0.004,
        head_width=0.05,
        color="tab:red",
        length_includes_head=True,
    )
    ax.text(
        u_star[0] + scale * grad_f[0] - 0.1,
        u_star[1] + scale * grad_f[1] + 0.05,
        r"$\nabla f(u^\star)$",
        color="tab:red",
    )

    ax.arrow(
        u_star[0],
        u_star[1],
        scale * grad_phi[0],
        scale * grad_phi[1],
        width=0.004,
        head_width=0.05,
        color="tab:blue",
        length_includes_head=True,
    )
    ax.text(
        u_star[0] + scale * grad_phi[0] + 0.03,
        u_star[1] + scale * grad_phi[1] + 0.02,
        r"$\nabla \varphi$",
        color="tab:blue",
    )

    ax.set_title("2D constrained minimization: geometry")
    ax.set_xlabel(r"$u_1$")
    ax.set_ylabel(r"$u_2$")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(u1.min(), u1.max())
    ax.set_ylim(u2.min(), u2.max())
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    # Keep the original canvas size (no tight crop) so paired figures match dimensions.
    save(fig, "01_minimization_geometry.png", tight=False)


def minimization_along_constraint() -> None:
    # Same model, reduced to one scalar variable along the feasible line.
    A = np.array([[3.0, 1.0], [1.0, 2.0]], dtype=float)
    c = np.array([1.0, -1.0], dtype=float)
    g = 0.5

    # Particular feasible point and tangent direction to c.u = g.
    u0 = np.array([g, 0.0], dtype=float)
    d = np.array([1.0, 1.0], dtype=float)  # c.d = 0

    t = np.linspace(-1.5, 1.5, 300)
    U = u0[None, :] + t[:, None] * d[None, :]
    vals = 0.5 * np.einsum("bi,ij,bj->b", U, A, U)

    t_star = t[int(np.argmin(vals))]

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.plot(t, vals, linewidth=2.0)
    ax.axvline(t_star, linestyle="--", linewidth=1.5)
    ax.set_title(r"$f(u(t))$ on the constraint $\varphi(u)=0$")
    ax.set_xlabel(r"$t$ (feasible direction parameter)")
    ax.set_ylabel(r"$f(u(t))$")
    ax.set_aspect("auto")
    ax.grid(True, alpha=0.25)
    ax.text(
        t_star + 0.04,
        float(np.min(vals)) + 0.03,
        rf"$t^\star \approx {t_star:.3f}$",
        fontsize=10,
    )
    # Keep the original canvas size (no tight crop) so paired figures match dimensions.
    save(fig, "01_minimization_on_constraint.png", tight=False)


def polar_single_active_constraint() -> None:
    # One active inequality varphi(u) >= 0 at ubar:
    # T(ubar) = {d : grad(varphi)(ubar) . d >= 0}, T(ubar)^o is a ray.
    grad_varphi = np.array([0.7, 0.35], dtype=float)
    polar_dir = -grad_varphi

    t = np.linspace(-2.0, 2.0, 200)
    perp = np.array([grad_varphi[1], -grad_varphi[0]], dtype=float)
    line = np.outer(t, perp)

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.axhline(0.0, color="0.35", linewidth=1.0)
    ax.axvline(0.0, color="0.35", linewidth=1.0)
    ax.plot(line[:, 0], line[:, 1], color="black", linewidth=1.5)

    ax.quiver(0, 0, grad_varphi[0], grad_varphi[1],
              angles="xy", scale_units="xy", scale=1, color="tab:blue")
    ax.quiver(0, 0, polar_dir[0], polar_dir[1],
              angles="xy", scale_units="xy", scale=1, color="tab:red")

    ax.text(grad_varphi[0] + 0.03, grad_varphi[1] + 0.03,
            r"$\nabla\varphi(\bar u)$", color="tab:blue")
    ax.text(polar_dir[0] - 0.02, polar_dir[1] - 0.06,
            r"$T(\bar u)^\circ$", color="tab:red")
    ax.text(-1.02, -0.68, r"$T(\bar u)=\{d:\nabla\varphi(\bar u)\cdot d\geq 0\}$")

    ax.set_aspect("equal", "box")
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xlabel(r"$d_1$")
    ax.set_ylabel(r"$d_2$")
    ax.set_title("Tangent halfspace and its polar ray")
    ax.grid(True, alpha=0.25)
    save(fig, "02_polar_single_active_constraint.png", tight=False)


def polar_first_quadrant() -> None:
    # K = {d: d1 >= 0, d2 >= 0}, K^o = {v: v1 <= 0, v2 <= 0}.
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.axhline(0.0, color="0.35", linewidth=1.0)
    ax.axvline(0.0, color="0.35", linewidth=1.0)

    # Rays for K
    ax.plot([0.0, 1.1], [0.0, 0.0], color="tab:blue", linewidth=2.0)
    ax.plot([0.0, 0.0], [0.0, 1.1], color="tab:blue", linewidth=2.0)
    ax.text(0.58, 0.25, r"$K$", color="tab:blue")

    # Rays for K^o
    ax.plot([0.0, -1.1], [0.0, 0.0], color="tab:red", linewidth=2.0)
    ax.plot([0.0, 0.0], [0.0, -1.1], color="tab:red", linewidth=2.0)
    ax.text(-0.72, -0.26, r"$K^\circ$", color="tab:red")

    ax.set_aspect("equal", "box")
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xlabel(r"$d_1$")
    ax.set_ylabel(r"$d_2$")
    ax.set_title("First quadrant cone and its polar")
    ax.grid(True, alpha=0.25)
    save(fig, "02_polar_first_quadrant.png", tight=False)


def polar_ray_halfspace() -> None:
    # K = {d = t e1 : t >= 0}; K^o = {v: v1 <= 0}.
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.axhline(0.0, color="0.35", linewidth=1.0)
    ax.axvline(0.0, color="0.35", linewidth=1.0)

    # Ray K along +d1
    ax.arrow(0.0, 0.0, 1.1, 0.0, width=0.01, head_width=0.07,
             color="tab:blue", length_includes_head=True)

    # Boundary of K^o: v1 = 0
    ax.plot([0.0, 0.0], [-1.1, 1.1], color="tab:red", linewidth=2.0)
    ax.text(0.82, -0.12, r"$K$", color="tab:blue")
    ax.text(-1.08, 0.92, r"$K^\circ=\{v_1\leq 0\}$", color="tab:red")

    ax.set_aspect("equal", "box")
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xlabel(r"$d_1$")
    ax.set_ylabel(r"$d_2$")
    ax.set_title("Ray and its polar halfspace")
    ax.grid(True, alpha=0.25)
    save(fig, "02_polar_ray_halfspace.png", tight=False)


def main() -> None:
    control_to_state_diagram()
    reduced_cost_1d()
    minimization_geometry()
    minimization_along_constraint()
    polar_single_active_constraint()
    polar_first_quadrant()
    polar_ray_halfspace()
    print(f"Wrote figures to: {OUT_DIR}")


if __name__ == "__main__":
    main()

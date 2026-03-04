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


def armijo_condition_plot() -> None:
    # 1D model phi(alpha) = f(x + alpha p) for a descent direction p.
    phi0 = 1.0
    dphi0 = -2.0
    curvature = 5.0
    c1 = 1e-4

    alpha = np.linspace(0.0, 1.2, 500)
    phi = phi0 + dphi0 * alpha + 0.5 * curvature * alpha**2
    armijo_rhs = phi0 + c1 * alpha * dphi0

    ok = phi <= armijo_rhs
    alpha_acc = alpha[np.where(ok)[0][-1]] if np.any(ok) else 0.0

    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.plot(alpha, phi, linewidth=2.0, label=r"$\phi(\alpha)=f(u_k+\alpha p_k)$")
    ax.plot(alpha, armijo_rhs, linewidth=2.0, linestyle="--",
            label=r"$\phi(0)+c_1\alpha\phi'(0)$")
    ax.axvline(alpha_acc, color="tab:green", linewidth=1.8,
               linestyle=":", label=rf"max accepted $\alpha\approx {alpha_acc:.2f}$")
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"value")
    ax.set_title("Armijo sufficient decrease condition")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper center", fontsize=9)
    save(fig, "03_armijo_condition.png")


def gd_path_quadratic() -> None:
    # Gradient descent path on an anisotropic quadratic.
    A = np.array([[8.0, 2.5], [2.5, 1.2]], dtype=float)
    A = 0.5 * (A + A.T)

    def f(x: np.ndarray) -> float:
        return 0.5 * x @ A @ x

    def grad(x: np.ndarray) -> np.ndarray:
        return A @ x

    x = np.array([1.6, 1.2], dtype=float)
    path = [x.copy()]
    c1 = 1e-4
    beta = 0.5

    for _ in range(25):
        g = grad(x)
        if np.linalg.norm(g) < 1e-10:
            break
        p = -g
        alpha = 1.0
        fx = f(x)
        while f(x + alpha * p) > fx + c1 * alpha * g @ p:
            alpha *= beta
            if alpha < 1e-12:
                break
        x = x + alpha * p
        path.append(x.copy())

    path = np.array(path)
    grid = np.linspace(-1.7, 1.8, 320)
    X, Y = np.meshgrid(grid, grid)
    Z = 0.5 * (A[0, 0] * X**2 + 2.0 * A[0, 1] * X * Y + A[1, 1] * Y**2)

    fig, ax = plt.subplots(figsize=(5.8, 5.8))
    ax.contour(X, Y, Z, levels=20, cmap="viridis")
    ax.plot(path[:, 0], path[:, 1], marker="o", linewidth=1.8,
            markersize=4, color="tab:red")
    ax.scatter(path[0, 0], path[0, 1], s=55, color="black", label="start")
    ax.scatter(0.0, 0.0, s=55, color="tab:green", label="minimizer")
    ax.set_title("Gradient descent path on anisotropic quadratic")
    ax.set_xlabel(r"$x_1$")
    ax.set_ylabel(r"$x_2$")
    ax.set_aspect("equal", "box")
    ax.grid(True, alpha=0.2)
    ax.legend(loc="upper right")
    save(fig, "03_gd_path_quadratic.png", tight=False)


def gd_vs_cg_convergence() -> None:
    # Convergence on SPD quadratic: compare steepest descent and linear CG.
    n = 80
    D = np.geomspace(1.0, 1e4, n)
    A = np.diag(D)
    b = np.ones(n)
    x0 = np.ones(n)

    def f(x: np.ndarray) -> float:
        return 0.5 * x @ A @ x - b @ x

    x_star = np.linalg.solve(A, b)
    f_star = f(x_star)

    # Steepest descent with exact line search for quadratics.
    x = x0.copy()
    gd = []
    for _ in range(80):
        g = A @ x - b
        gd.append(max(f(x) - f_star, 1e-18))
        if np.linalg.norm(g) < 1e-12:
            break
        p = -g
        alpha = (g @ g) / (p @ A @ p)
        x = x + alpha * p

    # Linear conjugate gradient.
    x = x0.copy()
    r = b - A @ x
    p = r.copy()
    cg = []
    for _ in range(80):
        cg.append(max(f(x) - f_star, 1e-18))
        rr = r @ r
        if np.sqrt(rr) < 1e-12:
            break
        Ap = A @ p
        alpha = rr / (p @ Ap)
        x = x + alpha * p
        r_new = r - alpha * Ap
        beta = (r_new @ r_new) / rr
        p = r_new + beta * p
        r = r_new

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.semilogy(gd, linewidth=2.0, label="GD (exact LS)")
    ax.semilogy(cg, linewidth=2.0, label="Linear CG")
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$f(x_k)-f(x^\star)$")
    ax.set_title("Convergence on SPD quadratic")
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="upper right")
    save(fig, "03_gd_vs_cg_convergence.png")


def bfgs_vs_gd_rosenbrock() -> None:
    # Compare GD and BFGS on Rosenbrock in 2D.
    def f(x: np.ndarray) -> float:
        return 100.0 * (x[1] - x[0]**2)**2 + (1.0 - x[0])**2

    def grad(x: np.ndarray) -> np.ndarray:
        return np.array([
            -400.0 * x[0] * (x[1] - x[0]**2) - 2.0 * (1.0 - x[0]),
            200.0 * (x[1] - x[0]**2),
        ])

    def armijo_step(x: np.ndarray, p: np.ndarray, g: np.ndarray,
                    c1: float = 1e-4, beta: float = 0.5) -> float:
        alpha = 1.0
        fx = f(x)
        while f(x + alpha * p) > fx + c1 * alpha * g @ p:
            alpha *= beta
            if alpha < 1e-14:
                break
        return alpha

    x0 = np.array([-1.2, 1.0], dtype=float)
    nmax = 80

    # GD history.
    x = x0.copy()
    gd = []
    for _ in range(nmax):
        g = grad(x)
        gd.append(np.linalg.norm(g))
        if gd[-1] < 1e-12:
            break
        p = -g
        alpha = armijo_step(x, p, g)
        x = x + alpha * p

    # BFGS history.
    x = x0.copy()
    H = np.eye(2)
    bfgs = []
    for _ in range(nmax):
        g = grad(x)
        bfgs.append(np.linalg.norm(g))
        if bfgs[-1] < 1e-12:
            break
        p = -H @ g
        if g @ p >= 0:
            p = -g
            H = np.eye(2)
        alpha = armijo_step(x, p, g)
        s = alpha * p
        x_new = x + s
        y = grad(x_new) - g
        ys = y @ s
        if ys > 1e-12:
            rho = 1.0 / ys
            I = np.eye(2)
            H = (I - rho * np.outer(s, y)) @ H @ (I - rho * np.outer(y, s)) \
                + rho * np.outer(s, s)
        else:
            H = np.eye(2)
        x = x_new

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.semilogy(gd, linewidth=2.0, label="GD + Armijo")
    ax.semilogy(bfgs, linewidth=2.0, label="BFGS + Armijo")
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\|\nabla f(x_k)\|_2$")
    ax.set_title("Rosenbrock: GD vs BFGS")
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="upper right")
    save(fig, "03_bfgs_vs_gd_rosenbrock.png")


def trust_region_geometry() -> None:
    # Show local quadratic model and trust-region ball around current iterate.
    def f(x1: np.ndarray, x2: np.ndarray) -> np.ndarray:
        return 0.5 * (x1**2 - 1.2 * x2**2) + 0.1 * (x1**4 + x2**4)

    xk = np.array([0.7, 0.8], dtype=float)
    delta = 0.55
    B = np.array([[1.0, 0.0], [0.0, -1.2]], dtype=float)
    g = np.array([xk[0] + 0.4 * xk[0]**3, -1.2 * xk[1] + 0.4 * xk[1]**3])

    # Cauchy-like and Newton-like trial steps for visualization.
    p_sd = -delta * g / max(np.linalg.norm(g), 1e-14)
    try:
        p_newton = -np.linalg.solve(B, g)
    except np.linalg.LinAlgError:
        p_newton = p_sd
    if np.linalg.norm(p_newton) > delta:
        p_newton = delta * p_newton / np.linalg.norm(p_newton)

    grid = np.linspace(-1.5, 1.8, 350)
    X, Y = np.meshgrid(grid, grid)
    Z = f(X, Y)

    fig, ax = plt.subplots(figsize=(5.8, 5.8))
    ax.contour(X, Y, Z, levels=22, cmap="cividis")
    ax.scatter(xk[0], xk[1], s=70, color="black", label=r"$x_k$")

    theta = np.linspace(0.0, 2.0 * np.pi, 400)
    circle = xk[:, None] + delta * np.vstack((np.cos(theta), np.sin(theta)))
    ax.plot(circle[0], circle[1], linestyle="--", linewidth=1.8,
            color="tab:blue", label=r"trust region $\|p\|\leq \Delta_k$")

    x_sd = xk + p_sd
    x_nt = xk + p_newton
    ax.scatter(x_sd[0], x_sd[1], s=55, color="tab:red", label="Cauchy-like step")
    ax.scatter(x_nt[0], x_nt[1], s=55, color="tab:green", label="model minimizer in ball")
    ax.plot([xk[0], x_sd[0]], [xk[1], x_sd[1]], color="tab:red", linewidth=1.5)
    ax.plot([xk[0], x_nt[0]], [xk[1], x_nt[1]], color="tab:green", linewidth=1.5)

    ax.set_title("Trust-region geometry")
    ax.set_xlabel(r"$x_1$")
    ax.set_ylabel(r"$x_2$")
    ax.set_xlim(-1.4, 1.8)
    ax.set_ylim(-1.4, 1.8)
    ax.set_aspect("equal", "box")
    ax.grid(True, alpha=0.2)
    ax.legend(loc="upper left", fontsize=8)
    save(fig, "03_trust_region_geometry.png", tight=False)


def main() -> None:
    control_to_state_diagram()
    reduced_cost_1d()
    minimization_geometry()
    minimization_along_constraint()
    polar_single_active_constraint()
    polar_first_quadrant()
    polar_ray_halfspace()
    armijo_condition_plot()
    gd_path_quadratic()
    gd_vs_cg_convergence()
    bfgs_vs_gd_rosenbrock()
    trust_region_geometry()
    print(f"Wrote figures to: {OUT_DIR}")


if __name__ == "__main__":
    main()

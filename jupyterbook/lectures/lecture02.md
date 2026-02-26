# Karush-Kuhn-Tucker Conditions for Inequality Constraints

## Overview

This lecture introduces first-order and second-order optimality conditions for finite-dimensional optimization problems with inequality constraints, expressed as
$$
\varphi_j(u)\ge 0, \qquad j=1,\dots,m.
$$

---

## 1. Problem Setting

Consider
$$
\min_{u\in\mathbb{R}^n} f(u)
\quad\text{s.t.}\quad
 \varphi_j(u)\ge 0,\; j=1,\dots,m,
$$
where $f,\varphi_j\in C^1(\mathbb{R}^n)$.

Define the feasible set
$$
\mathcal{U}_{\mathrm{ad}}:=\{u\in\mathbb{R}^n: \varphi_j(u)\ge 0\ \forall j\}.
$$

Let $\bar u\in\mathcal{U}_{\mathrm{ad}}$ be a local minimizer and define the active set
$$
\mathcal{A}(\bar u):=\{j\in\{1,\dots,m\}: \varphi_j(\bar u)=0\}.
$$

---

## 2. Geometry and First-Order Necessary Condition

### 2.1 Tangent Cone (Linearized)

At $\bar u$, the linearized feasible directions are
$$
T_{\mathcal{U}_{\mathrm{ad}}}(\bar u)
:=
\{d\in\mathbb{R}^n: \nabla \varphi_j(\bar u)\cdot d\ge 0\ \forall j\in\mathcal{A}(\bar u)\}.
$$

### 2.2 Variational Inequality Form

A first-order necessary condition at a local minimizer is
$$
\nabla f(\bar u)\cdot d\ge 0
\quad\forall d\in T_{\mathcal{U}_{\mathrm{ad}}}(\bar u).
$$

Equivalently,
$$
-\nabla f(\bar u)\in T_{\mathcal{U}_{\mathrm{ad}}}(\bar u)^\circ,
$$
where the *polar cone* of a set $T\subseteq\mathbb{R}^n$ is defined by
$$
T^\circ:=\{v\in\mathbb{R}^n:\; v\cdot d\le 0\ \forall d\in T\}.
$$

### 2.3 Geometric examples of polar cones

#### Example A: one active inequality

For one active constraint $\varphi(\bar u)=0$, the linearized tangent cone is
$$
T(\bar u)=\{d:\nabla\varphi(\bar u)\cdot d\ge 0\},
$$
a halfspace. Its polar is the ray
$$
T(\bar u)^\circ=\{-\lambda\nabla\varphi(\bar u):\lambda\ge 0\}.
$$

<img src="../slides/assets/02_polar_single_active_constraint.png" alt="Tangent halfspace and polar ray" width="50%">

This picture is the key local mechanism behind KKT: the normal direction comes from active constraints.

#### Example B: first quadrant cone

Let
$$
K=\{d:d_1\ge 0,\ d_2\ge 0\}.
$$
Then
$$
K^\circ=\{v:v_1\le 0,\ v_2\le 0\}.
$$

<img src="../slides/assets/02_polar_first_quadrant.png" alt="First quadrant cone and its polar" width="50%">

So the polar collects vectors making non-positive scalar product with every direction in $K$.

#### Example C: a ray and its polar halfspace

Let
$$
K=\{d=t\,e_1:\ t\ge 0\}.
$$
Then
$$
K^\circ=\{v:v_1\le 0\},
$$
while $v_2$ is free.

<img src="../slides/assets/02_polar_ray_halfspace.png" alt="Ray and its polar halfspace" width="50%">

This is the dual picture of Example A: halfspace $\leftrightarrow$ ray.

---

## 3. KKT Conditions

Under suitable regularity (also known as **constraint qualification** or **CQ**, see below), there exist multipliers
$\lambda_j\ge 0$ such that
$$
\nabla f(\bar u)-\sum_{j=1}^m \lambda_j\nabla \varphi_j(\bar u)=0.
$$

Indeed, under CQ one can represent the polar as
$$ T_{\mathcal{U}_{\mathrm{ad}}}(\bar u)^\circ = \left\{-\sum_{j\in\mathcal A(\bar u)}\lambda_j\nabla\varphi_j(\bar u):\lambda_j\ge 0\right\},$$
and therefore
$$
-\nabla f(\bar u)\in T_{\mathcal{U}_{\mathrm{ad}}}(\bar u)^\circ
\;\Longleftrightarrow\;
\nabla f(\bar u)-\sum_{j\in\mathcal A(\bar u)}\lambda_j\nabla\varphi_j(\bar u)=0.
$$

The multiplier $\lambda_j$ is associated with the $j$-th constraint $\varphi_j(u)\ge 0$, and the full KKT system is:
$$
\begin{aligned}
\nabla f(\bar u)-\sum_{j=1}^m \lambda_j\nabla \varphi_j(\bar u) &= 0,\\
\varphi_j(\bar u) &\ge 0,\quad j=1,\dots,m,\\
\lambda_j &\ge 0,\quad j=1,\dots,m,\\
\lambda_j \varphi_j(\bar u) &= 0,\quad j=1,\dots,m.
\end{aligned}
$$

Complementarity can be interpreted as:

- if $\varphi_j(\bar u)>0$ (the $j$-th constraint is inactive), then $\lambda_j=0$;
- if $\lambda_j>0$, then $\varphi_j(\bar u)=0$ (the $j$-th constraint is active).

---

## 4. Constraint Qualifications

KKT multipliers do not automatically exist at every local minimizer.

### 4.1 LICQ (**Linear Independence CQ**)

LICQ holds at $\bar u$ if
$$
\{\nabla \varphi_j(\bar u): j\in\mathcal{A}(\bar u)\}
$$
are linearly independent, where $\mathcal{A}(\bar u)$ is the active set of constraints at $\bar u$, i.e., the set of indices $j$ such that $\varphi_j(\bar u)=0$.

Consequences:

- existence of KKT multipliers;
- uniqueness of multipliers.

### 4.2 MFCQ (**Mangasarian-Fromovitz CQ**)

MFCQ holds at $\bar u$ if there exists $d\in\mathbb{R}^n$ such that
$$
\nabla \varphi_j(\bar u)\cdot d>0
\quad\forall j\in\mathcal{A}(\bar u).
$$

Consequences:

- existence of KKT multipliers;
- multipliers may be non-unique.

LICQ implies MFCQ.

### 4.3 Slater Condition (Convex Problems)

If $f$ is convex, each $\varphi_j$ is concave (so $\varphi_j(u)\ge 0$ defines a convex feasible set), and there exists $u_0$ such that
$$
\varphi_j(u_0)>0\quad\forall j,
$$
then KKT conditions are necessary and sufficient, and strong duality holds.

---

## 5. Second-Order Conditions

Assume KKT holds at $(\bar u,\lambda)$ and define the critical cone
$$
\mathcal{C}(\bar u,\lambda)
:=
\left\{d\in\mathbb{R}^n:
\begin{array}{l}
\nabla \varphi_j(\bar u)\cdot d=0\ \text{for } j\in\mathcal{A}(\bar u)\ \text{with } \lambda_j>0,\\
\nabla \varphi_j(\bar u)\cdot d\ge 0\ \text{for } j\in\mathcal{A}(\bar u)\ \text{with } \lambda_j=0
\end{array}
\right\}.
$$

Use the Lagrangian
$$
\mathcal{L}(u,\lambda):=f(u)-\sum_{j=1}^m \lambda_j \varphi_j(u).
$$

Second-order necessary condition:
$$
d^T\nabla^2_{uu}\mathcal{L}(\bar u,\lambda)d\ge 0
\quad\forall d\in\mathcal{C}(\bar u,\lambda).
$$

Second-order sufficient condition:
$$
d^T\nabla^2_{uu}\mathcal{L}(\bar u,\lambda)d>0
\quad\forall d\in\mathcal{C}(\bar u,\lambda)\setminus\{0\}
$$
which implies that $\bar u$ is a strict local minimizer.

---

## 6. Example: Box Constraints

Consider
$$
\min_{u\in\mathbb{R}^n} f(u)
\quad\text{s.t.}\quad
u_i-a_i\ge 0,\quad b_i-u_i\ge 0,
\quad i=1,\dots,n.
$$

Componentwise KKT conditions become
$$
\begin{cases}
\partial_i f(\bar u)\ge 0 & \text{if } \bar u_i=a_i,\\
\partial_i f(\bar u)\le 0 & \text{if } \bar u_i=b_i,\\
\partial_i f(\bar u)=0 & \text{if } a_i<\bar u_i<b_i.
\end{cases}
$$

This is the finite-dimensional prototype of bound-constrained control variables in PDE-constrained optimization.

---

## 7. Summary

At a local minimizer:

1. first-order optimality is a variational inequality on feasible directions;
2. with suitable CQ (LICQ/MFCQ), this becomes KKT;
3. complementarity links active constraints and positive multipliers;
4. second-order conditions refine local optimality.

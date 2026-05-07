# General Optimality Systems for PDE-Constrained Optimization

## Overview

The previous lectures developed a complete theory for linear-quadratic elliptic optimal control:

- reduced formulation through the control-to-state map;
- adjoint-based gradient representation;
- variational inequality for control constraints;
- discrete KKT system and its numerical solution.

This lecture places these results in an abstract framework that covers
nonlinear state equations and general objective functionals.

The logical path of the lecture:

1. fix the functional analytic setting — spaces, duality pairings, and adjoint operators;
2. define Fréchet differentiability of maps between Banach spaces;
3. derive the differentiability of the control-to-state map via the Implicit Function Theorem;
4. compute the Fréchet derivative of the reduced functional via the chain rule;
5. introduce the adjoint equation and rewrite the gradient without nested PDE solves;
6. formulate the Lagrangian, derive the full KKT system, and state the second-order conditions;
7. define the normal cone precisely and recover the variational inequality and projection formula;
8. compare the reduced and all-at-once approaches in terms of structure and computational cost;
9. verify the framework on a semilinear elliptic example.

---

## Functional Analytic Setting

Let $X$, $Y$, $Z$, $U$ be real Banach spaces.
A **Banach space** is a complete normed vector space.
A **Hilbert space** is a Banach space whose norm is induced by an inner product:
$\|v\|^2 = (v,v)$.
Hilbert spaces are reflexive, meaning $(Y')' \cong Y$ isometrically.

**Dual spaces and duality pairings.**
The **dual space** $X'$ consists of all bounded linear functionals $\ell: X \to \mathbb{R}$.
For $\ell \in X'$ and $x \in X$ we write the **duality pairing**
$$
\langle \ell, x \rangle_{X', X} := \ell(x).
$$
The norm of $X'$ is $\|\ell\|_{X'} := \sup_{\|x\|_X \le 1} |\langle \ell, x \rangle|$.

**Spaces used in this lecture.**

- $Y$: the **state space**.
  For elliptic problems: $Y = H_0^1(\Omega)$ or $Y = H^1(\Omega)$.
  In the abstract framework: a Hilbert space.

- $U$: the **control space**.
  Typically: $U = L^2(\Omega)$ or $U = L^2(\Gamma)$.
  In the abstract framework: a Hilbert space.

- $Z'$: the **constraint residual space**, i.e., the codomain of the state equation.
  For elliptic problems: $Z' = H^{-1}(\Omega)$, the dual of $H_0^1(\Omega)$.
  In the abstract framework: the dual of some Banach space $Z$.

- $p \in Z$: the **adjoint state** or **Lagrange multiplier** for the state equation.

**Adjoint operators.**
Let $A: X \to Z'$ be a bounded linear operator.
Its **adjoint** $A^*: Z \to X'$ is defined by
$$
\langle Ax, z \rangle_{Z', Z}
= \langle x, A^* z \rangle_{X, X'}
\qquad \forall x \in X, \; z \in Z.
$$
If $X$ and $Z$ are Hilbert spaces and we identify $Z \cong Z'$ via the Riesz isomorphism,
this becomes $(Ax, z)_Z = (x, A^* z)_X$.

**Riesz isomorphism.**
In a Hilbert space $X$, the **Riesz map** $\mathcal{R}_X: X \to X'$ is defined by
$$
\langle \mathcal{R}_X u, v \rangle_{X', X} = (u, v)_X \qquad \forall v \in X.
$$
It is an isometric isomorphism.
Its inverse $\mathcal{R}_X^{-1}: X' \to X$ identifies the **gradient**:
given $\ell \in X'$, the unique element $g := \mathcal{R}_X^{-1} \ell \in X$ satisfies
$$
\langle \ell, v \rangle_{X', X} = (g, v)_X \qquad \forall v \in X,
$$
and we write $g = \nabla f(u)$ when $\ell = f'(u)$.

---

## Fréchet Differentiability

Let $X$ and $W$ be Banach spaces, and let $F: X \to W$.
$F$ is **Fréchet differentiable** at $x \in X$ if there exists a bounded linear operator
$F'(x): X \to W$ such that
$$
\lim_{\|h\|_X \to 0}
\frac{\|F(x+h) - F(x) - F'(x)h\|_W}{\|h\|_X} = 0.
$$
The operator $F'(x)$ is called the **Fréchet derivative** of $F$ at $x$.

For real-valued functionals $F: X \to \mathbb{R}$, the Fréchet derivative is an element of $X'$:
$$
F'(x) \in X'.
$$

**Partial derivatives.**
For a map $F: X_1 \times X_2 \to W$, the **partial Fréchet derivative** with respect to $X_1$
at $(x_1, x_2)$ in direction $h_1 \in X_1$ is
$$
F_{x_1}(x_1, x_2) h_1
:= \lim_{t \to 0} \frac{F(x_1 + t h_1, x_2) - F(x_1, x_2)}{t}.
$$
For fixed $(x_1, x_2)$, the map $h_1 \mapsto F_{x_1}(x_1, x_2) h_1$ is a bounded linear operator
from $X_1$ to $W$.

**Chain rule.**
If $G: X \to W$ and $H: W \to V$ are Fréchet differentiable, then
$(H \circ G): X \to V$ is Fréchet differentiable and
$$
(H \circ G)'(x) = H'(G(x)) \circ G'(x).
$$
For our purposes: given $S: U \to Y$ and $J: Y \times U \to \mathbb{R}$, with $f(u) := J(S(u), u)$,

$$
f'(u)h
= \langle J_y(S(u), u), S'(u) h \rangle_{Y', Y}
+ \langle J_u(S(u), u), h \rangle_{U', U}
\qquad \forall h \in U.
$$

---

## Abstract Problem

We now state the abstract optimal control problem and fix notation for all derivatives.

**Objective functional.**
$$
J: Y \times U \to \mathbb{R}
$$
is assumed Fréchet differentiable, with partial derivatives
$$
J_y(y, u) \in Y', \qquad J_u(y, u) \in U'.
$$

**Constraint mapping.**
$$
e: Y \times U \to Z'
$$
is assumed Fréchet differentiable, with partial derivatives
$$
e_y(y, u): Y \to Z', \qquad e_u(y, u): U \to Z'.
$$
These are bounded linear operators for each fixed $(y, u)$.

**Admissible controls.**
$$
U_{\mathrm{ad}} \subset U
$$
is a nonempty, closed, and convex subset.

**Optimal control problem.**
Find
$$
(\bar y, \bar u) \in Y \times U_{\mathrm{ad}}
$$
such that
$$
J(\bar y, \bar u)
= \inf\bigl\{ J(y, u):\; e(y, u) = 0 \text{ in } Z',\; u \in U_{\mathrm{ad}} \bigr\}.
$$

**Standing assumptions.**

1. For every $u \in U_{\mathrm{ad}}$, the equation $e(y, u) = 0$ in $Z'$ admits a unique solution $y \in Y$.
2. $J$ and $e$ are of class $C^1$ (all Fréchet derivatives exist and are continuous).
3. At every feasible pair $(y, u)$, the operator $e_y(y, u): Y \to Z'$ is an isomorphism.

Assumption 3 is the regularity (constraint qualification) condition.
For the Poisson state equation, it holds because $-\Delta: H_0^1(\Omega) \to H^{-1}(\Omega)$ is a homeomorphism.

---

## Control-to-State Map

Under Assumption 1, define
$$
S: U_{\mathrm{ad}} \to Y, \qquad u \mapsto y = S(u),
$$
where $S(u)$ is the unique solution of $e(y, u) = 0$ in $Z'$.
The map $S$ is the **control-to-state operator** (or **solution operator**) of the state equation.

**Differentiability via the Implicit Function Theorem.**
The Implicit Function Theorem in Banach spaces states:
if $e: Y \times U \to Z'$ is $C^1$ and $e_y(y_0, u_0): Y \to Z'$ is an isomorphism at some
feasible point $(y_0, u_0)$, then locally around $u_0$ the map $u \mapsto S(u)$ is $C^1$ and
$$
S'(u)h = -e_y(y, u)^{-1} e_u(y, u) h \qquad \forall h \in U,
$$
where $y = S(u)$.

**Derivation.**
Differentiate the identity $e(S(u), u) = 0$ with respect to $u$ in direction $h$:
$$
e_y(y, u)\bigl[S'(u) h\bigr] + e_u(y, u)\bigl[h\bigr] = 0 \qquad \text{in } Z'.
$$
Solving for $S'(u)h$ using invertibility of $e_y$:
$$
S'(u) h = -e_y(y, u)^{-1} e_u(y, u) h.
$$

This formula shows that:

- the state variation $S'(u)h$ due to a control perturbation $h$ is obtained by one linearized state solve;
- the linearized state equation is $e_y(y,u)[v] = -e_u(y,u)[h]$, driven by the perturbation;
- the map $h \mapsto S'(u)h$ is a bounded linear operator from $U$ to $Y$.

**Example (linear Poisson model).**
For $e(y, u) = -\Delta y - u$:
$$
e_y(y, u) = -\Delta : H_0^1(\Omega) \to H^{-1}(\Omega), \qquad
e_u(y, u) = -I : L^2(\Omega) \to H^{-1}(\Omega).
$$
Hence $S'(u)h = (-\Delta)^{-1}h = S(h)$,
consistent with the linearity of $S$ in that case.

---

## Reduced Functional and Chain Rule

Define the **reduced functional**
$$
f: U_{\mathrm{ad}} \to \mathbb{R}, \qquad f(u) := J(S(u), u).
$$
The PDE-constrained problem becomes an optimization problem in $u$ alone:
$$
\min_{u \in U_{\mathrm{ad}}} f(u).
$$

**Fréchet derivative of $f$.**
By the chain rule,

$$
f'(u) h
= \langle J_y(y, u),\, S'(u) h \rangle_{Y', Y}
+ \langle J_u(y, u),\, h \rangle_{U', U}
\qquad \forall h \in U,
$$

where $y = S(u)$.
Substituting the formula for $S'(u)h$:

$$
f'(u) h
= -\langle J_y(y, u),\, e_y(y, u)^{-1} e_u(y, u) h \rangle_{Y', Y}
+ \langle J_u(y, u),\, h \rangle_{U', U}.
$$

This expression is not yet computationally convenient: evaluating the first term for every
direction $h$ requires one linearized state solve per direction tested.
The **adjoint equation** removes this bottleneck by moving the operator to the other argument.

---

## Adjoint Equation

**Setup.**
Let $(y, u)$ be a feasible pair.
We want to rewrite the directional derivative $f'(u)h$ so that $h$ appears only in one
bounded linear functional, without intermediate state solves.

**Definition of the adjoint state.**
Introduce $p \in Z$ as the solution of the **adjoint equation**
$$
e_y(y, u)^* p = J_y(y, u) \qquad \text{in } Y'.
$$
The adjoint state $p$ is uniquely determined because $e_y(y, u)^*: Z \to Y'$ is
an isomorphism (dual of the isomorphism $e_y(y,u): Y \to Z'$).

Note: $p$ lives in $Z$, the pre-dual of $Z'$.
For the Poisson model, $Z = Z' = H^{-1}(\Omega)$ and the Riesz identification sends $p$ to $H_0^1(\Omega)$.

**Gradient formula via the adjoint.**
With $p$ as above, use the duality identity: for any invertible $A: Y \to Z'$,
$$
\langle q, A^{-1} w \rangle_{Y', Y}
= \langle (A^*)^{-1} q, w \rangle_{Z, Z'} \qquad \forall q \in Y', \; w \in Z'.
$$
Apply this with $A = e_y(y,u)$, $q = J_y(y, u)$, $w = e_u(y, u) h$, and $p = (e_y^*)^{-1} J_y$:
$$
\langle J_y(y, u),\, e_y(y, u)^{-1} e_u(y, u) h \rangle_{Y', Y}
= \langle p,\, e_u(y, u) h \rangle_{Z, Z'}
= \langle e_u(y, u)^* p,\, h \rangle_{U', U}.
$$
Hence
$$
f'(u) h
= \langle J_u(y, u) - e_u(y, u)^* p,\; h \rangle_{U', U}
\qquad \forall h \in U.
$$

**Gradient of the reduced functional.**
Identifying $f'(u) \in U'$ with its Riesz representative $\nabla f(u) \in U$:
$$
\nabla f(u)
= \mathcal{R}_U^{-1}\bigl(J_u(y, u) - e_u(y, u)^* p\bigr)
\in U.
$$

**Cost count.** Regardless of how many directions $h$ must be tested:

- one state solve to get $y = S(u)$;
- one adjoint solve to get $p$;
- one evaluation of $J_u - e_u^* p$.

**Example (Poisson model).**
$J_y = y - y_d \in L^2(\Omega)$, $e_y^* = -\Delta$, so the adjoint equation is
$-\Delta p = y - y_d$ with $p = 0$ on $\partial\Omega$.
$J_u = \alpha u \in L^2(\Omega)$, $e_u^* = -I$, so $\nabla f(u) = p + \alpha u$,
exactly as derived in Lecture 4.

---

## Lagrangian Formulation

The Lagrangian naturally encodes both the objective and the constraint.

**Definition.**
$$
\mathcal{L}: Y \times U \times Z \to \mathbb{R}, \qquad
\mathcal{L}(y, u, p) := J(y, u) - \langle e(y, u), p \rangle_{Z', Z}.
$$
Here $p \in Z$ is the **Lagrange multiplier** for the equality constraint $e(y, u) = 0$.
The duality pairing is well-defined because $e(y, u) \in Z'$ and $p \in Z$.

**Partial derivatives of $\mathcal{L}$.**

With respect to $y$ in direction $v \in Y$:

$$
\mathcal{L}_y(y, u, p)[v]
= \langle J_y(y, u),\, v \rangle_{Y', Y}
- \langle e_y(y, u) v,\, p \rangle_{Z', Z}
= \langle J_y(y, u) - e_y(y, u)^* p,\; v \rangle_{Y', Y}.
$$

With respect to $u$ in direction $h \in U$:

$$
\mathcal{L}_u(y, u, p)[h]
= \langle J_u(y, u),\, h \rangle_{U', U}
- \langle e_u(y, u) h,\, p \rangle_{Z', Z}
= \langle J_u(y, u) - e_u(y, u)^* p,\; h \rangle_{U', U}.
$$

With respect to $p$ in direction $q \in Z$:

$$
\mathcal{L}_p(y, u, p)[q] = -\langle e(y, u),\, q \rangle_{Z', Z}.
$$

**Interpretation of the conditions $\mathcal{L}_p = 0$, $\mathcal{L}_y = 0$.**

- $\mathcal{L}_p(y, u, p) = 0$ means $e(y, u) = 0$ in $Z'$: the state equation.
- $\mathcal{L}_y(y, u, p) = 0$ means $e_y(y, u)^* p = J_y(y, u)$ in $Y'$: the adjoint equation.
- $\mathcal{L}_u(y, u, p) = 0$ (unconstrained case) means $J_u - e_u^* p = 0$ in $U'$: optimality.

---

## First-Order Optimality System

We now state the full KKT conditions.

**Unconstrained case ($U_{\mathrm{ad}} = U$).**
At a local minimizer $\bar u$ with corresponding state $\bar y = S(\bar u)$,
there exists a unique adjoint state $\bar p \in Z$ such that:
$$
\begin{cases}
e(\bar y, \bar u) = 0 & \text{in } Z', \\[4pt]
e_y(\bar y, \bar u)^* \bar p = J_y(\bar y, \bar u) & \text{in } Y', \\[4pt]
J_u(\bar y, \bar u) - e_u(\bar y, \bar u)^* \bar p = 0 & \text{in } U'.
\end{cases}
$$

**Constrained case ($U_{\mathrm{ad}} \subsetneq U$).**
At a local minimizer $\bar u \in U_{\mathrm{ad}}$, the optimality condition for $u$ becomes
the variational inequality
$$
f'(\bar u)(u - \bar u) \ge 0 \qquad \forall u \in U_{\mathrm{ad}},
$$
which in terms of the Lagrangian is
$$
0 \in \mathcal{L}_u(\bar y, \bar u, \bar p) + N_{U_{\mathrm{ad}}}(\bar u) \qquad \text{in } U'.
$$

The full constrained KKT system is therefore:
$$
\begin{cases}
e(\bar y, \bar u) = 0 & \text{in } Z', \\[4pt]
e_y(\bar y, \bar u)^* \bar p = J_y(\bar y, \bar u) & \text{in } Y', \\[4pt]
0 \in J_u(\bar y, \bar u) - e_u(\bar y, \bar u)^* \bar p + N_{U_{\mathrm{ad}}}(\bar u) & \text{in } U'.
\end{cases}
$$

The three equations live, respectively, in $Z'$, $Y'$, and $U'$.
The three unknowns are $\bar y \in Y$, $\bar p \in Z$, and $\bar u \in U_{\mathrm{ad}}$.

---

## Normal Cone and Variational Inequality

**Definition.**
Let $K \subset U$ be a nonempty closed convex set and $u \in K$.
The **normal cone** to $K$ at $u$ is the closed convex cone
$$
N_K(u) :=
\bigl\{\xi \in U' :\; \langle \xi, v - u \rangle_{U', U} \le 0 \quad \forall v \in K \bigr\}.
$$
Geometrically: $\xi \in N_K(u)$ if and only if $u$ maximizes the linear functional
$v \mapsto \langle \xi, v \rangle$ over $K$ (i.e., $\xi$ points outward from $K$ at $u$).

**Equivalence with the variational inequality.**
The normal-cone inclusion $-f'(\bar u) \in N_{U_{\mathrm{ad}}}(\bar u)$ is equivalent to
$$
f'(\bar u)(u - \bar u) \ge 0 \qquad \forall u \in U_{\mathrm{ad}}.
$$

**Projection.** When $U$ is a Hilbert space and the Riesz identification $U \cong U'$ is used,
for any $w \in U$ the projection $\Pi_K(w)$ is the unique element of $K$ satisfying
$$
\Pi_K(w) \in K, \qquad
(w - \Pi_K(w),\, v - \Pi_K(w))_U \le 0 \qquad \forall v \in K.
$$
The normal-cone inclusion reads: $w - u \in N_K(u)$ if and only if $u = \Pi_K(w)$.

**Projection formula for the control.**
The KKT condition $0 \in \nabla f(\bar u) + N_{U_{\mathrm{ad}}}(\bar u)$ is equivalent to
$$
\bar u = \Pi_{U_{\mathrm{ad}}}(\bar u - \rho\, \nabla f(\bar u))
\qquad \text{for any } \rho > 0.
$$
This is both a characterization of optimality and a natural fixed-point iteration
(projected gradient step).

**Box constraints in $L^2(\Omega)$.**
For
$$
U_{\mathrm{ad}} = \{u \in U :\; u_{\min}(x) \le u(x) \le u_{\max}(x) \text{ a.e.}\}
$$
with $u_{\min}, u_{\max} \in L^\infty(\Omega)$, the normal cone at $\bar u \in U_{\mathrm{ad}}$ is
characterized pointwise: $\xi \in N_{U_{\mathrm{ad}}}(\bar u)$ if and only if for a.e. $x \in \Omega$,
$$
\begin{cases}
\xi(x) \ge 0 & \text{if } \bar u(x) = u_{\min}(x), \\
\xi(x) = 0 & \text{if } u_{\min}(x) < \bar u(x) < u_{\max}(x), \\
\xi(x) \le 0 & \text{if } \bar u(x) = u_{\max}(x).
\end{cases}
$$
For the linear-quadratic case, the KKT condition $0 \in \alpha\bar u + \bar p + N_{U_{\mathrm{ad}}}(\bar u)$
is equivalent to the pointwise projection
$$
\bar u(x) = \Pi_{[u_{\min}(x),\, u_{\max}(x)]}
\!\left(-\tfrac{1}{\alpha} \bar p(x)\right)
\quad \text{a.e. in } \Omega.
$$

---

## Semismooth Newton for Control Constraints

The projection characterization of the KKT condition suggests a nonsmooth root equation.
Fix $\rho > 0$ and define

$$
F(u)
:=
u - \Pi_{U_{\mathrm{ad}}}\bigl(u - \rho\, \nabla f(u)\bigr)
\in U.
$$

Then

$$
F(\bar u)=0
\iff
\bar u = \Pi_{U_{\mathrm{ad}}}(\bar u - \rho\, \nabla f(\bar u))
\iff
0 \in \nabla f(\bar u) + N_{U_{\mathrm{ad}}}(\bar u).
$$

Hence solving the constrained first-order condition is equivalent to solving
the nonsmooth equation $F(u)=0$ in $U$.

### Proof of the three equivalences for $F(\bar u)=0$

We prove the equivalence chain in detail.

Fix $K := U_{\mathrm{ad}}$, and define

$$
w := \bar u - \rho\, \nabla f(\bar u).
$$

1. **First equivalence**

  $$
  F(\bar u)=0
  \iff
  \bar u - \Pi_K(w)=0
  \iff
  \bar u = \Pi_K(w).
  $$

  This follows directly from the definition
  $F(u)=u-\Pi_K(u-\rho\nabla f(u))$.

1. **Second equivalence**
  By the projection-normal-cone characterization in a Hilbert space,
  for any $u \in K$ and $w \in U$:

  $$
  u = \Pi_K(w)
  \iff
  w-u \in N_K(u).
  $$
  
  Applying this with $u=\bar u$ and $w=\bar u-\rho\nabla f(\bar u)$ gives
  
  $$
  \bar u = \Pi_K(\bar u-\rho\nabla f(\bar u))
  \iff
  -\rho\nabla f(\bar u) \in N_K(\bar u).
  $$
  
  Since $N_K(\bar u)$ is a cone, $\eta \in N_K(\bar u)$ implies
  $(1/\rho)\eta \in N_K(\bar u)$ for every $\rho>0$, hence
  
  $$
  -\rho\nabla f(\bar u) \in N_K(\bar u)
  \iff
  -\nabla f(\bar u) \in N_K(\bar u)
  \iff
  0 \in \nabla f(\bar u)+N_K(\bar u).
  $$

Combining 1 and 2 yields

$$
F(\bar u)=0
\iff
\bar u = \Pi_{U_{\mathrm{ad}}}(\bar u - \rho\nabla f(\bar u))
\iff
0 \in \nabla f(\bar u)+N_{U_{\mathrm{ad}}}(\bar u).
$$

### Semismoothness and generalized derivatives

In finite dimensions, a locally Lipschitz map $G: \mathbb{R}^n \to \mathbb{R}^m$
is **semismooth** at $x$ if for every sequence $h_k \to 0$ and every choice
$V_k \in \partial_C G(x+h_k)$ (Clarke generalized Jacobian),
$$
\|G(x+h_k)-G(x)-V_k h_k\| = o(\|h_k\|).
$$
It is **strongly semismooth** if the remainder is $O(\|h_k\|^2)$.

In Hilbert spaces, the same idea is used with generalized derivatives of
nonsmooth operators (Clarke/Bouligand derivatives or Newton derivatives).
The key property is that each linearization captures the first-order behavior
with a small remainder, so a Newton correction can be defined.

For box constraints, the projection is piecewise affine pointwise:
$$
\Pi_{[a,b]}(s)=\min\{\max\{s,a\},b\}.
$$
Therefore $\Pi_{[a,b]}$ is globally Lipschitz and strongly semismooth in
finite dimensions, and so is the induced superposition operator in
$L^q(\Omega)$ spaces ($1<q<\infty$).

### Semismooth Newton step

Given an iterate $u^k$, choose a generalized derivative
$$
M_k \in \partial F(u^k),
$$
and compute $\delta u^k \in U$ from
$$
M_k\, \delta u^k = -F(u^k),
$$
then update
$$
u^{k+1} = u^k + \delta u^k.
$$

For constrained PDE-control, each Newton step requires coupled evaluation
of state and adjoint variables because $\nabla f(u)$ depends on
$(y(u),p(u))$. In practice, one linearizes the reduced mapping
$u \mapsto F(u)$ and solves the resulting linear system with PDE blocks
using the same solver technology as in the all-at-once KKT setting.

### Primal-dual construction with an additional multiplier for $F(u)=0$

The semismooth Newton step can be obtained from a KKT system with an explicit
multiplier for the equation constraint $F(u)=0$.

At iteration $k$, consider the proximal feasibility problem

$$
\min_{u\in U}\; \frac12\|u-u^k\|_U^2
\qquad \text{subject to } F(u)=0.
$$

Its Lagrangian is

$$
\mathscr L_k(u,\lambda)
= \frac12\|u-u^k\|_U^2 + \langle \lambda, F(u)\rangle_{U,U},
$$

where $\lambda \in U$ is an additional Lagrange multiplier.

The first-order system is
$$
\begin{cases}
u-u^k + F'(u)^*\lambda = 0,\\
F(u)=0.
\end{cases}
$$

Replace $F'(u^k)$ with a generalized derivative
$M_k\in\partial F(u^k)$ and linearize at $(u^k,\lambda^k)$:

$$
\begin{pmatrix}
I & M_k^* \\
M_k & 0
\end{pmatrix}
\begin{pmatrix}
\delta u^k \\
\delta\lambda^k
\end{pmatrix}
=-
\begin{pmatrix}
u^k-u^k+M_k^*\lambda^k \\
F(u^k)
\end{pmatrix}
=-
\begin{pmatrix}
M_k^*\lambda^k \\
F(u^k)
\end{pmatrix}.
$$

If we initialize and keep $\lambda^k=0$, the first block row gives
$\delta u^k = -M_k^*\delta\lambda^k$, and substituting in the second row yields
the Schur-complement equation
$$
M_k M_k^*\,\delta\lambda^k = F(u^k),
$$
with
$$
\delta u^k = -M_k^*\delta\lambda^k.
$$
This primal-dual step is equivalent to solving a regularized Newton equation
for $F(u)=0$ and leads to the same local model as the direct semismooth
Newton correction when $M_k$ is nonsingular.

Hence semismooth Newton may be viewed in two equivalent ways:

- **reduced form**: solve $M_k\delta u^k=-F(u^k)$;
- **primal-dual form**: solve a saddle-point system in $(\delta u^k,\delta\lambda^k)$
  with the extra multiplier enforcing the linearized equality constraint.

The primal-dual perspective is useful for block preconditioning and for
embedding the step in all-at-once PDE-constrained solvers.

### Explicit form for box-constrained linear-quadratic control

For
$$
J(y,u)=\tfrac12\|y-y_d\|_{L^2}^2 + \tfrac\alpha2\|u\|_{L^2}^2,
\qquad
U_{\mathrm{ad}} = \{u: u_{\min}\le u \le u_{\max}\},
$$
we have
$$
\nabla f(u)=\alpha u + p(u),
$$
and the root equation is
$$
F(u)=u-\Pi_{[u_{\min},u_{\max}]}
\bigl(u-\rho(\alpha u + p(u))\bigr)=0.
$$

Define the active and inactive sets at iteration $k$ from
$$
w^k := u^k - \rho(\alpha u^k + p^k):
$$
$$
\mathcal A_-^k := \{x:\, w^k(x) \le u_{\min}(x)\},
\qquad
\mathcal A_+^k := \{x:\, w^k(x) \ge u_{\max}(x)\},
\qquad
\mathcal I^k := \Omega \setminus (\mathcal A_-^k \cup \mathcal A_+^k).
$$
Then the generalized derivative of the projection is multiplication by
the characteristic function $\chi_{\mathcal I^k}$, and the Newton step is
equivalent to:

- enforce control bounds directly on active sets,
  $$
  u^{k+1}=u_{\min} \text{ on } \mathcal A_-^k,
  \qquad
  u^{k+1}=u_{\max} \text{ on } \mathcal A_+^k;
  $$
- solve the unconstrained optimality equation on $\mathcal I^k$.

This is exactly the primal-dual active-set (PDAS) method.
Thus, for box constraints, PDAS can be interpreted as a semismooth Newton
method applied to the projection equation.

### Local convergence statement

Under standard assumptions:

1. $\nabla f$ is (locally) Lipschitz differentiable in a neighborhood of $\bar u$;
2. $F$ is semismooth at $\bar u$;
3. every generalized derivative $M \in \partial F(\bar u)$ is uniformly invertible,
   i.e., $\|M^{-1}\| \le c$.

Then semismooth Newton is locally superlinearly convergent:
$$
\|u^{k+1}-\bar u\|_U = o(\|u^k-\bar u\|_U).
$$
If $F$ is strongly semismooth and the derivative approximation is exact,
one obtains local quadratic convergence:
$$
\|u^{k+1}-\bar u\|_U = O(\|u^k-\bar u\|_U^2).
$$

These rates explain why active-set/Newton-type methods are often much faster
than projected gradient iterations once the active set is close to the final one.

---

## Reduced vs All-at-Once

The first-order optimality conditions can be approached in two complementary ways.

**Reduced approach.**
The control $\bar u \in U_{\mathrm{ad}}$ is the only optimization unknown.
The state $\bar y = S(\bar u)$ and adjoint $\bar p$ are functions of $\bar u$.
The problem reduces to
$$
\min_{u \in U_{\mathrm{ad}}} f(u),
$$
and algorithms iterate on $u$ only, solving state and adjoint equations at each step.

At each optimization iterate $u^k$:

1. state solve: find $y^k \in Y$ such that $e(y^k, u^k) = 0$;
2. adjoint solve: find $p^k \in Z$ such that $e_y(y^k, u^k)^* p^k = J_y(y^k, u^k)$;
3. gradient: $g^k = \mathcal{R}_U^{-1}(J_u(y^k, u^k) - e_u(y^k, u^k)^* p^k) \in U$;
4. control update: $u^{k+1} = \Pi_{U_{\mathrm{ad}}}(u^k - \tau_k g^k)$.

Advantages:

- the optimization variable lives only in $U$;
- off-the-shelf iterative solvers can be reused for state/adjoint;
- memory footprint is one control vector.

Disadvantages:

- state and adjoint must be solved to sufficient accuracy at each step;
- second-order information (Hessian of $f$) requires a second adjoint solve (second-order adjoint);
- slow convergence when inner solves are expensive.

**All-at-once approach.**
The triple $(\bar y, \bar u, \bar p)$ is the unknown.
One applies Newton's method directly to the full KKT system.
The Newton step $(\delta y, \delta u, \delta p)$ solves the block system

$$
\begin{pmatrix}
\mathcal{L}_{yy} & \mathcal{L}_{yu} & e_y^* \\[4pt]
\mathcal{L}_{uy} & \mathcal{L}_{uu} & e_u^* \\[4pt]
e_y & e_u & 0
\end{pmatrix}
\begin{pmatrix}
\delta y \\[4pt] \delta u \\[4pt] \delta p
\end{pmatrix}
= -
\begin{pmatrix}
\mathcal{L}_y \\[4pt] \mathcal{L}_u + N_{U_{\mathrm{ad}}}^{-1}(\bar u) \\[4pt] e
\end{pmatrix},
$$

where $\mathcal{L}_{yy}$, $\mathcal{L}_{yu}$, $\mathcal{L}_{uu}$ are second partial derivatives
of the Lagrangian and $N_{U_{\mathrm{ad}}}^{-1}$ is the linearization of the normal-cone
inclusion (active-set method or semismooth Newton regularization).

For the linear-quadratic case, the block system is the saddle-point system seen in Lecture 7,
which is linear and can be solved directly.
For nonlinear problems, Newton linearization is needed and the block structure is the same,
but with variable-dependent operators.

Advantages:

- quadratic convergence near the solution (Newton's method);
- one global solve per Newton iteration rather than many inner iterations;
- second-order information built in.

Disadvantages:

- the system lives in $Y \times U \times Z$, which is much larger than $U$ alone;
- saddle-point preconditioning is needed for efficiency.

---

## A Worked Example: Semilinear Elliptic Control

We verify that the abstract framework applies to a nonlinear state equation.

**Model problem.**
Let $\Omega \subset \mathbb{R}^d$ be a bounded Lipschitz domain and $d: \mathbb{R} \to \mathbb{R}$
a smooth function with $d(0) = 0$ (the **semilinearity**).
Consider

$$
\min_{u \in U_{\mathrm{ad}}} J(y, u)
:= \frac{1}{2}\|y - y_d\|_{L^2(\Omega)}^2
+ \frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2
$$

subject to

$$
\begin{cases}
-\Delta y + d(y) = u & \text{in } \Omega, \\[4pt]
y = 0 & \text{on } \partial\Omega.
\end{cases}
$$

**Spaces.**

$$
Y = H_0^1(\Omega), \quad U = L^2(\Omega), \quad Z' = H^{-1}(\Omega), \quad Z = H_0^1(\Omega).
$$

**Constraint mapping.**

$$
e: H_0^1(\Omega) \times L^2(\Omega) \to H^{-1}(\Omega),
$$

defined by

$$
\langle e(y, u), v \rangle_{H^{-1}, H_0^1}
:= \int_\Omega \nabla y \cdot \nabla v \, dx
+ \int_\Omega d(y) v \, dx
- \int_\Omega u v \, dx
\qquad \forall v \in H_0^1(\Omega).
$$

**Partial derivatives of $e$.**

Linear operator $e_y(y, u): H_0^1(\Omega) \to H^{-1}(\Omega)$:

$$
\langle e_y(y, u) w, v \rangle_{H^{-1}, H_0^1}
= \int_\Omega \nabla w \cdot \nabla v \, dx
+ \int_\Omega d'(y) w\, v \, dx
\qquad \forall w, v \in H_0^1(\Omega).
$$

Under the assumption $d'(y) \ge -\sigma_0$ for some $\sigma_0 < \lambda_1(\Omega)$
(first Dirichlet eigenvalue of $-\Delta$), coercivity of the bilinear form is preserved,
so $e_y(y, u)$ is an isomorphism (Assumption 3 is satisfied).

Linear operator $e_u(y, u): L^2(\Omega) \to H^{-1}(\Omega)$:

$$
\langle e_u(y, u) h, v \rangle_{H^{-1}, H_0^1}
= -\int_\Omega h\, v \, dx
= -(h, v)_{L^2(\Omega)}
\qquad \forall h \in L^2(\Omega), \; v \in H_0^1(\Omega).
$$

Thus $e_u(y, u) h = -h$ (embedding of $L^2$ into $H^{-1}$ via the $L^2$ inner product).

**Partial derivatives of $J$.**

$$
J_y(y, u) = y - y_d \in L^2(\Omega) \subset H^{-1}(\Omega) = Y',
\qquad
J_u(y, u) = \alpha u \in L^2(\Omega) = U'.
$$

**Adjoint operators.**

$e_y(y, u)^*: H_0^1(\Omega) \to H^{-1}(\Omega)$ (identifying $Z \cong H_0^1$ via Riesz):

$$
\langle e_y(y, u)^* p, v \rangle_{H^{-1}, H_0^1}
= \int_\Omega \nabla v \cdot \nabla p \, dx
+ \int_\Omega d'(y)\, v\, p \, dx
\qquad \forall p, v \in H_0^1(\Omega).
$$

(The symmetry of the bilinear form and the symmetry of multiplication by $d'(y)$ mean
$e_y^* = e_y$ in this case; this would differ for non-symmetric operators such as transport.)

$e_u(y, u)^*: H_0^1(\Omega) \to L^2(\Omega) = U'$:

$$
\langle e_u(y, u)^* p, h \rangle_{U', U}
= \langle p, e_u(y, u) h \rangle_{Z, Z'}
= -(p, h)_{L^2(\Omega)}
\qquad \forall h \in L^2(\Omega),
$$

so $e_u(y, u)^* p = -p \in L^2(\Omega)$ (restriction of the $H_0^1$ function $p$ to $L^2$).

**Adjoint equation.**
Find $\bar p \in H_0^1(\Omega)$ such that

$$
e_y(\bar y, \bar u)^* \bar p = J_y(\bar y, \bar u) = \bar y - y_d
\quad \text{in } H^{-1}(\Omega),
$$

i.e., in weak form:

$$
\int_\Omega \nabla v \cdot \nabla \bar p \, dx
+ \int_\Omega d'(\bar y)\, v\, \bar p \, dx
= \int_\Omega (\bar y - y_d) v \, dx
\qquad \forall v \in H_0^1(\Omega).
$$

The corresponding strong form is

$$
\begin{cases}
-\Delta \bar p + d'(\bar y)\, \bar p = \bar y - y_d & \text{in } \Omega, \\[4pt]
\bar p = 0 & \text{on } \partial\Omega.
\end{cases}
$$

```{admonition} Key observation
:class: important

The adjoint equation is **always linear** in the adjoint state $\bar p$,
even when the state equation is nonlinear in $y$.
This is because the adjoint equation involves only the linearization $e_y$,
not the full nonlinear mapping $e$.
```

**Gradient of the reduced functional.**

$$
\nabla f(\bar u)
= J_u(\bar y, \bar u) - e_u(\bar y, \bar u)^* \bar p
= \alpha \bar u - (-\bar p)
= \alpha \bar u + \bar p \in L^2(\Omega).
$$

**Full optimality system for box constraints.**

$$
\begin{cases}
-\Delta \bar y + d(\bar y) = \bar u & \text{in } \Omega, \\[4pt]
-\Delta \bar p + d'(\bar y)\, \bar p = \bar y - y_d & \text{in } \Omega, \\[4pt]
\bar u = \Pi_{[u_{\min},\, u_{\max}]}\!\left(-\tfrac{1}{\alpha} \bar p\right)
  & \text{a.e. in } \Omega, \\[4pt]
\bar y = \bar p = 0 & \text{on } \partial\Omega.
\end{cases}
$$

The structure is the same as the linear-quadratic case from Lectures 4 and 5.
The nonlinearity only affects the state equation (via $d(\bar y)$) and
the adjoint equation (via $d'(\bar y)\bar p$).
The projection formula for the control is unchanged.

---

## Existence for the General Problem

For completeness we state when a minimizer exists in the abstract setting.

**Theorem.**
Suppose:

1. $U_{\mathrm{ad}}$ is nonempty, bounded, closed, and convex in the reflexive Banach space $U$;
2. the reduced functional $f: U_{\mathrm{ad}} \to \mathbb{R}$ is weakly lower semicontinuous;
3. $f$ is bounded below on $U_{\mathrm{ad}}$.

Then there exists at least one global minimizer $\bar u \in U_{\mathrm{ad}}$.

**Remarks.**

- If $U_{\mathrm{ad}}$ is not bounded, coercivity of $f$ (level sets bounded) replaces Condition 1.
- Condition 2 holds if $f$ is convex and continuous (e.g., linear-quadratic case),
  or more generally if $J$ is convex and the state map $S$ is weakly continuous.
- Strict convexity of $f$ implies uniqueness of the minimizer.

---

## Second-Order Conditions

First-order conditions are necessary but not sufficient for a strict local minimizer.

**Lagrangian Hessian.**
At a KKT point $(\bar y, \bar u, \bar p)$, the second-order derivative of $\mathcal{L}$
with respect to $(y, u)$ in the direction $(v, h) \in Y \times U$ is

$$
\mathcal{L}''(\bar y, \bar u, \bar p)[(v, h), (v, h)]
= \langle J_{yy}(\bar y, \bar u) v, v \rangle_{Y', Y}
+ 2\langle J_{yu}(\bar y, \bar u) h, v \rangle_{Y', Y}
+ \langle J_{uu}(\bar y, \bar u) h, h \rangle_{U', U}
- \langle e_{yy}(\bar y, \bar u)[v, v],\, \bar p \rangle_{Z', Z}.
$$

The term involving $e_{yy}$ vanishes for linear state equations,
so for the Poisson model the Hessian of the Lagrangian coincides with the Hessian of $J$.

**Critical cone.**

$$
C(\bar u) :=
\bigl\{h \in U :\; h \text{ satisfies the first-order admissibility conditions at } \bar u,\;
f'(\bar u) h = 0 \bigr\}.
$$

**Second-order necessary condition.**
If $\bar u$ is a local minimizer, then for every $h \in C(\bar u)$ and $v = S'(\bar u) h$,

$$
\mathcal{L}''(\bar y, \bar u, \bar p)[(v, h), (v, h)] \ge 0.
$$

**Second-order sufficient condition.**
If $(\bar y, \bar u, \bar p)$ is a KKT point and

$$
\mathcal{L}''(\bar y, \bar u, \bar p)[(v, h), (v, h)] > 0
\qquad \forall (v, h) \in Y \times C(\bar u),\; (v, h) \ne (0, 0),
$$

then $\bar u$ is a strict local minimizer.
For the linear-quadratic case with $\alpha > 0$ this condition holds globally:

$$
\mathcal{L}''[(v,h),(v,h)] = \|v\|_{L^2}^2 + \alpha\|h\|_{L^2}^2 > 0,
$$

confirming that the linear-quadratic problem has a unique global minimizer.

---

## Summary

The abstract framework identifies three functional spaces and two adjoint operators
that are the backbone of every PDE-constrained optimality system.

- State : $\bar y \in Y$ satisfies $e(\bar y, \bar u) = 0$ in $Z'$.
- Adjoint : $\bar p \in Z$ satisfies $e_y(\bar y, \bar u)^* \bar p = J_y(\bar y, \bar u)$ in $Y'$.
- Control : $\bar u \in U_{\mathrm{ad}}$ satisfies
  $0 \in J_u(\bar y, \bar u) - e_u(\bar y, \bar u)^* \bar p + N_{U_{\mathrm{ad}}}(\bar u)$ in $U'$.

The reduced gradient is
$$
\nabla f(\bar u) = \mathcal{R}_U^{-1}\bigl(J_u(\bar y, \bar u) - e_u(\bar y, \bar u)^* \bar p\bigr) \in U.
$$

For the linear-quadratic elliptic model of Lectures 4 and 5, all operators are linear,
the Lagrangian Hessian is positive definite, and the system reduces exactly to the
adjoint-based formulation derived there.
For the semilinear model of this lecture, the structure is identical but the
state equation and adjoint equation carry additional nonlinear terms.

## References

- F. Tröltzsch, *Optimal Control of Partial Differential Equations*, AMS, 2010 —
  Chapters 4 and 5 (semilinear elliptic control, second-order conditions)
- A. Manzoni, A. Quarteroni, S. Salsa, *Optimal Control of Partial Differential Equations*,
  Springer, 2021 — Chapters 3 and 4
- J. C. De los Reyes, *Numerical PDE-Constrained Optimization*, Springer, 2015 —
  Chapters 2 and 3
- M. Hinze, R. Pinnau, M. Ulbrich, S. Ulbrich, *Optimization with PDE Constraints*,
  Springer, 2009 — Chapter 1 (abstract framework)

# Finite Element Discretization of Elliptic Optimal Control

## Overview

The previous lectures developed the continuous theory for linear-quadratic elliptic
optimal control, derived the adjoint-based gradient formula, and established the
first-order optimality system in both unconstrained and box-constrained settings.

Before introducing finite elements, we place the elliptic control problem
in a general analytical framework:

1. a general variational formulation for elliptic state equations;
2. well-posedness of the control-to-state map;
3. passage from infinite-dimensional analysis to finite-dimensional approximation.

---

## A General Elliptic Boundary Value Problem

Let $\Omega\subset\mathbb R^d$ be a bounded Lipschitz domain with boundary $\Gamma$.
We consider the second-order linear elliptic operator

$$
A y
:=
-\nabla\cdot(\mathbf A(x)\nabla y)
-\mathbf b(x)\cdot\nabla y
+c(x)y,
$$
with coefficients satisfying standard assumptions:

- $\mathbf A\in L^\infty(\Omega;\mathbb R^{d\times d})$ is uniformly elliptic, i.e. there exist
  constants $0<\alpha\le M<\infty$ such that

  $$
  \alpha |\xi|^2 \le \xi^T \mathbf A(x)\xi \le M |\xi|^2
  \qquad \forall \xi\in\mathbb R^d,\ \text{for q.o. }x\in\Omega;
  $$
- $\mathbf b\in L^\infty(\Omega)^d$ and

  $$
  \|\mathbf b\|_{L^\infty(\Omega)}\le \beta;
  $$
- $c\in L^\infty(\Omega)$ and

  $$
  \|c\|_{L^\infty(\Omega)}\le \gamma;
  $$
- $\operatorname{div}\mathbf b\in L^\infty(\Omega)$;
- the lower-order part is not too negative, in the sense that there exists
  $\alpha_0\ge 0$ such that

  $$
  \frac12 \operatorname{div}\mathbf b(x) + c(x)\ge -\alpha_0
  \qquad \text{for q.o. } x\in\Omega,
  $$

  and $\alpha_0$ is small enough compared with the ellipticity constant $\alpha$, namely

  $$
  \alpha_0 < \frac{\alpha}{C_P^2},
  $$

  where $C_P$ is the Poincare constant of $\Omega$.

We impose homogeneous Dirichlet boundary conditions

$$
y=0 \qquad \text{on }\Gamma.
$$

We work in the standard Sobolev space

$$
V:=H_0^1(\Omega).
$$

The weak formulation of the state equation is:
find $y\in V$ such that

$$
a(y,v)=(u,v)_{L^2(\Omega)}+\langle F,v\rangle_{V',V}
\qquad \forall v\in V,
$$
where the bilinear form is

$$
\langle Ay, v\rangle_{V',V} := a(y,v)
=
\int_\Omega \mathbf A\nabla y\cdot \nabla v\,dx
- \int_\Omega \mathbf b\cdot \nabla y\,v\,dx
+ \int_\Omega c\,y\,v\,dx.
$$

Special case:

$$
\mathbf A(x)=I, \qquad \mathbf b(x)=0, \qquad c(x)=0.
$$

The model studied in the previous lectures is the simplest representative of this class.

---

## Trace Theorem and Boundary Data

Since $\Omega$ is Lipschitz, the trace operator

$$
\gamma:H^s(\Omega)\to H^{s-\frac12}(\Gamma),
\qquad s>\frac12,
$$
is continuous.

Consequences:

- boundary conditions are meaningful in the Sobolev setting;
- boundary control problems can also be formulated rigorously, even though in this lecture
  we focus on distributed controls.

For $s=1$:

$$
\gamma:H^1(\Omega)\to H^{1/2}(\Gamma),
$$

The trace operator admits a continuous right inverse:
there exists a bounded lifting operator

$$
E:H^{s-\frac12}(\Gamma)\to H^s(\Omega),
\qquad s>\frac12,
$$

such that

$$
\gamma(Eg)=g.
$$

This yields the standard lifting of nonhomogeneous Dirichlet data.
For homogeneous Dirichlet data, the variational space is $H_0^1(\Omega)$.

---

## Continuity and Coercivity of the Bilinear Form

### Continuity

Since $\mathbf A$, $\mathbf b$, and $c$ are bounded, each term is continuous.

$$
\left|\int_\Omega \mathbf A\nabla y\cdot \nabla v\,dx\right|
\le M \|\nabla y\|_{L^2(\Omega)} \|\nabla v\|_{L^2(\Omega)}.
$$

Transport term:

$$
\left|\int_\Omega \mathbf b\cdot \nabla y\,v\,dx\right|
\le \beta \|\nabla y\|_{L^2(\Omega)} \|v\|_{L^2(\Omega)}.
$$

Zero-order term:

$$
\left|\int_\Omega c\,y\,v\,dx\right|
\le \gamma \|y\|_{L^2(\Omega)} \|v\|_{L^2(\Omega)}.
$$

For $y,v\in H_0^1(\Omega)$, Poincare's inequality yields

$$
\|y\|_{L^2(\Omega)}\le C_P \|\nabla y\|_{L^2(\Omega)},
\qquad
\|v\|_{L^2(\Omega)}\le C_P \|\nabla v\|_{L^2(\Omega)}.
$$

Hence

$$
|a(y,v)|
\le
\left(M+\beta C_P+\gamma C_P^2\right)
\|\nabla y\|_{L^2(\Omega)}\|\nabla v\|_{L^2(\Omega)}.
$$

On $H_0^1(\Omega)$, the norm $\|w\|_{H^1(\Omega)}$ is equivalent to
$\|\nabla w\|_{L^2(\Omega)}$. Therefore

$$
|a(y,v)|\le C \|y\|_{H^1(\Omega)}\|v\|_{H^1(\Omega)}
\qquad \forall y,v\in V.
$$

with the explicit choice

$$
C=M+\beta C_P+\gamma C_P^2.
$$

### Coercivity

Let $v\in V=H_0^1(\Omega)$. Then

$$
a(v,v)
=
\int_\Omega \mathbf A\nabla v\cdot \nabla v\,dx
- \int_\Omega \mathbf b\cdot \nabla v\,v\,dx
+ \int_\Omega c\,v^2\,dx.
$$

Uniform ellipticity gives

$$
\int_\Omega \mathbf A\nabla v\cdot \nabla v\,dx
\ge \alpha \|\nabla v\|_{L^2(\Omega)}^2.
$$

Integrating the first-order term by parts:

$$
-\int_\Omega \mathbf b\cdot \nabla v\,v\,dx
=
-\frac12 \int_\Omega \mathbf b\cdot \nabla(v^2)\,dx
=
\frac12 \int_\Omega (\operatorname{div}\mathbf b)\,v^2\,dx,
$$

The boundary term vanishes since $v=0$ on $\Gamma$.
Hence

$$
a(v,v)
=
\int_\Omega \mathbf A\nabla v\cdot \nabla v\,dx
+ \int_\Omega \left(c+\frac12\operatorname{div}\mathbf b\right) v^2\,dx.
$$

Under the assumption

$$
c(x)+\frac12\operatorname{div}\mathbf b(x)\ge -\alpha_0
\qquad \text{for q.o. }x\in\Omega,
$$

$$
a(v,v)\ge \alpha \|\nabla v\|_{L^2(\Omega)}^2-\alpha_0 \|v\|_{L^2(\Omega)}^2.
$$

Poincare's inequality gives

$$
\|v\|_{L^2(\Omega)}\le C_P \|\nabla v\|_{L^2(\Omega)},
$$

Hence

$$
a(v,v)\ge (\alpha-\alpha_0 C_P^2)\|\nabla v\|_{L^2(\Omega)}^2.
$$

If

$$
\alpha_0 < \frac{\alpha}{C_P^2},
$$

the right-hand side is strictly positive. Since on $H_0^1(\Omega)$ the
seminorm $\|\nabla v\|_{L^2(\Omega)}$ is equivalent to the full norm,

$$
a(v,v)\ge \alpha_0 \|v\|_{H^1(\Omega)}^2
\qquad \forall v\in V,
$$

for some $\alpha_0>0$.

---

## Nonhomogeneous Dirichlet and Neumann Conditions

Consider

$$
\begin{cases}
-\nabla\cdot(\mathbf A(x)\nabla y)-\mathbf b(x)\cdot\nabla y+c(x)y=f
& \text{in }\Omega,\\
y=g_D & \text{on }\Gamma_D,\\
(\mathbf A\nabla y)\cdot n=g_N & \text{on }\Gamma_N,
\end{cases}
$$

with

$$
\Gamma=\Gamma_D\cup\Gamma_N,
\qquad
\Gamma_D\cap\Gamma_N=\varnothing.
$$

If $g_D\neq 0$, one works in an affine space with prescribed trace, or reduces to the homogeneous case by lifting. If Neumann data are present, the boundary term appears on the right-hand side of the weak formulation.

For pure Dirichlet data on $\Gamma$, let

$$
\gamma y=g_D,
\qquad
g_D\in H^{1/2}(\Gamma).
$$

Choose a lifting

$$
w:=Eg_D\in H^1(\Omega),
\qquad
\gamma w=g_D,
\qquad
\|w\|_{H^1(\Omega)}\le C_{\mathrm{tr}}\|g_D\|_{H^{1/2}(\Gamma)}.
$$

Set

$$
z:=y-w.
$$

Then $z\in H_0^1(\Omega)$ and

$$
a(z,v)=(u,v)_{L^2(\Omega)}+\langle F,v\rangle_{V',V}-a(w,v)
\qquad \forall v\in H_0^1(\Omega).
$$

Define

$$
\ell(v):=(u,v)_{L^2(\Omega)}+\langle F,v\rangle_{V',V}-a(w,v).
$$

By continuity of $a(\cdot,\cdot)$,

$$
|a(w,v)|
\le C \|w\|_{H^1(\Omega)}\|v\|_{H^1(\Omega)}
\le C\,C_{\mathrm{tr}}\|g_D\|_{H^{1/2}(\Gamma)}\|v\|_{H^1(\Omega)}.
$$

Hence

$$
|\ell(v)|
\le
\Bigl(
C_u\|u\|_{L^2(\Omega)}
+ \|F\|_{V'}
+ C\,C_{\mathrm{tr}}\|g_D\|_{H^{1/2}(\Gamma)}
\Bigr)\|v\|_{H^1(\Omega)}.
$$

Lax-Milgram gives

$$
\|z\|_{H^1(\Omega)}
\le
C\Bigl(
\|u\|_{L^2(\Omega)}
+ \|F\|_{V'}
+ \|g_D\|_{H^{1/2}(\Gamma)}
\Bigr),
$$

and therefore

$$
\|y\|_{H^1(\Omega)}
\le
C\Bigl(
\|u\|_{L^2(\Omega)}
+ \|F\|_{V'}
+ \|g_D\|_{H^{1/2}(\Gamma)}
\Bigr).
$$

---

## Well-Posedness of the State Equation

By Lax-Milgram, for every $u\in L^2(\Omega)$ and every $F\in V'$,
there exists a unique state $y\in V$ such that

$$
a(y,v)=(u,v)_{L^2(\Omega)}+\langle F,v\rangle_{V',V}
\qquad \forall v\in V.
$$

The solution depends continuously on the data.
This defines the control-to-state operator

$$
S:L^2(\Omega)\to H_0^1(\Omega),
\qquad
u\mapsto y.
$$

---

## Cea and Bramble-Hilbert Lemmas

Let $V_h\subset V$ be a conforming finite-dimensional space and let $y_h\in V_h$ be the
Galerkin approximation of the state:

$$
a(y_h,v_h)=(u,v_h)_{L^2(\Omega)}+\langle F,v_h\rangle_{V',V}
\qquad \forall v_h\in V_h.
$$

Cea's lemma gives the quasi-optimality estimate

$$
\|y-y_h\|_{H^1(\Omega)}
\le
\frac{C}{\alpha_0}\inf_{v_h\in V_h}\|y-v_h\|_{H^1(\Omega)}.
$$

The finite element error is controlled, up to the factor $C/\alpha_0$, by the best
approximation error in the chosen discrete space.

The approximation term

$$
\inf_{v_h\in V_h}\|y-v_h\|_{H^1(\Omega)}
$$

is estimated by the Bramble-Hilbert lemma.
For standard finite element spaces of polynomial degree $k$ on a shape-regular mesh,
if the exact solution is sufficiently smooth, for example $y\in H^{k+1}(\Omega)$, then
the interpolation error satisfies

$$
\|y-I_h y\|_{H^m(\Omega)}
\le
C h^{k+1-m}|y|_{H^{k+1}(\Omega)},
\qquad 0\le m\le k+1.
$$

Combining this with Cea's lemma gives

$$
\|y-y_h\|_{H^1(\Omega)}
\le
C h^k |y|_{H^{k+1}(\Omega)}.
$$

- the weak formulation is posed in Sobolev spaces;
- Galerkin approximations inherit stability from coercivity;
- approximation theory converts best approximation into concrete mesh-dependent error bounds.

---

## Linear-Quadratic Optimal Control Problem

Let

$$
Q:=L^2(\Omega).
$$

Consider

$$
\min_{(y,u)} J(y,u)
:=
\frac12\|y-y_d\|_{L^2(\Omega)}^2
+\frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2,
$$

subject to

$$
a(y,v)=(u,v)_{L^2(\Omega)}+\langle F,v\rangle_{V',V}
\qquad \forall v\in V.
$$

The state equation defines the control-to-state map

$$
S:Q\to V,
\qquad
u\mapsto y.
$$

The reduced functional is

$$
f(u):=J(Su,u).
$$

The reduced problem is

$$
\min_{u\in Q} f(u).
$$

---

## Continuous Optimality System

State equation:

$$
a(y,v)=(u,v)_{L^2(\Omega)}+\langle F,v\rangle_{V',V}
\qquad \forall v\in V.
$$

Adjoint equation:

$$
a(v,p)=(y-y_d,v)_{L^2(\Omega)}
\qquad \forall v\in V.
$$

Gradient condition in the unconstrained case:

$$
\alpha u+p=0
\qquad \text{in }L^2(\Omega).
$$

Equivalently, in dual notation,

$$
\alpha R_Q u+B^*p=0
\qquad \text{in }Q'.
$$

If control constraints are present, this last equation is replaced by the corresponding
variational inequality or projection formula.

---

## Optimize-Then-Discretize and Discretize-Then-Optimize

Two routes are available.

### Optimize-Then-Discretize

1. derive the continuous state equation;
2. derive the continuous adjoint equation;
3. derive the continuous optimality system;
4. discretize the resulting system.

### Discretize-Then-Optimize

1. discretize the state equation;
2. define the discrete cost functional;
3. formulate the finite-dimensional optimization problem;
4. derive its KKT conditions.

In linear-quadratic elliptic settings these routes often lead to the same algebraic system. In more general settings they may differ.

---

## Choice of Discrete Spaces

Let $\mathcal T_h$ be a finite element mesh of $\Omega$.
We choose a finite-dimensional subspace

$$
V_h\subset V=H_0^1(\Omega)
$$
for the state and adjoint variables.

Typical choice:

- $V_h$ = continuous piecewise polinomial finite elements.

For the control, there are several possible choices.
The two most common are:

1. choose a discrete control space

   $$
   Q_h\subset Q=L^2(\Omega),
   $$

   for example piecewise polynomials or continuous piecewise polynomials;
2. do not discretize the control explicitly at first, and only discretize the state and adjoint.

In this lecture we focus first on the standard fully discrete setting

$$
V_h\subset V,\qquad Q_h\subset Q.
$$

Why are these choices important?

- the state equation must be well posed in the discrete space;
- the control space determines the number and meaning of the optimization variables;
- the discrete gradient and the discrete KKT system depend on the pairing between $V_h$ and $Q_h$.

The passage from the continuous problem to a numerical method is not only
the replacement of functions by vectors:
it also requires a choice of finite-dimensional spaces.

---

## Discrete State Equation

Given $u_h\in Q_h$, the finite element state equation is:
find $y_h\in V_h$ such that

$$
a(y_h,v_h)=(u_h,v_h)_Q+\langle F,v_h\rangle_{V',V}
\qquad \forall v_h\in V_h.
$$

This is the Galerkin approximation of the continuous weak problem.

Choose bases

$$
V_h=\operatorname{span}\{\varphi_1,\dots,\varphi_{n_y}\},
\qquad
Q_h=\operatorname{span}\{\psi_1,\dots,\psi_{n_u}\}.
$$

Expand

$$
y_h=\sum_{i=1}^{n_y} y_i\varphi_i,
\qquad
u_h=\sum_{j=1}^{n_u} u_j\psi_j.
$$

Then the discrete state equation becomes the linear system

$$
A y - B u = f,
$$
where

$$
A_{ij}=a(\varphi_j,\varphi_i),
\qquad
B_{ij}=(\psi_j,\varphi_i)_Q,
\qquad
f_i=\langle F,\varphi_i\rangle_{V',V}.
$$

Interpretation of the matrices:

- $A$ is the stiffness matrix of the PDE;
- $B$ transfers the control unknowns into the state equation;
- $f$ is the load vector.

If $Q_h=V_h$ and the same basis is used, then $B$ is simply the mass matrix.
If $Q_h$ is piecewise constant, then $B$ is a rectangular coupling matrix.

The important structural point is:
for each discrete control vector $u$, there is a unique discrete state vector $y$,
provided the stiffness matrix $A$ is invertible.

The discrete control-to-state map is

$$
S_h:Q_h\to V_h,
\qquad
u_h\mapsto y_h.
$$

---

## Discrete Cost Functional

Once the state has been discretized, the cost functional becomes

$$
J_h(y_h,u_h)
=
\frac12\|y_h-y_{d,h}\|_{L^2(\Omega)}^2
+\frac{\alpha}{2}\|u_h\|_{L^2(\Omega)}^2.
$$

In matrix form, this reads

$$
J_h(y,u)
=
\frac12 (y-y_d)^T M_y (y-y_d)
+\frac{\alpha}{2} u^T M_u u,
$$
where

$$
(M_y)_{ij}=(\varphi_j,\varphi_i)_Q,
\qquad
(M_u)_{ij}=(\psi_j,\psi_i)_Q.
$$

The fully discrete constrained problem is

$$
\min_{(y,u)\in \mathbb R^{n_y}\times\mathbb R^{n_u}} J_h(y,u)
\qquad \text{subject to} \qquad
A y-Bu=f.
$$

The conceptual path of Lecture 1 reappears in finite element form:

1. choose a control vector $u$;
2. solve the discrete state equation to get $y=S_h(u)$;
3. evaluate the discrete objective;
4. apply finite-dimensional optimization ideas.

---

## Reduced Discrete Problem

Since the discrete state is uniquely determined by the discrete control,
we may eliminate $y$ and define the reduced discrete functional

$$
f_h(u):=J_h(S_hu,u).
$$

Then the discrete optimal control problem becomes

$$
\min_{u\in \mathbb R^{n_u}} f_h(u),
$$
or, with bounds,

$$
\min_{u\in U_{\mathrm{ad}}^h} f_h(u),
$$
where $U_{\mathrm{ad}}^h\subset \mathbb R^{n_u}$ is the discrete admissible set.

Remarks:

- the unknown is now a vector of control coefficients;
- each evaluation of $f_h(u)$ still requires the solution of one PDE-like linear system;
- all numerical optimization methods from Lecture 3 now apply directly to $f_h$.

From the algorithmic point of view, the reduced finite element problem is a standard
optimization problem in which function and gradient evaluations require PDE solves.

---

## Discrete Adjoint Equation

To compute the gradient of $f_h$, we follow the continuous theory.
Given $u_h\in Q_h$, first solve the discrete state equation for $y_h$, and then define
the discrete adjoint $p_h\in V_h$ by

$$
a(v_h,p_h)=(y_h-y_{d,h},v_h)_Q
\qquad \forall v_h\in V_h.
$$

In matrix form, this becomes

$$
A^T p = M_y(y-y_d).
$$

For symmetric elliptic operators such as Poisson with standard Dirichlet conditions,
the stiffness matrix is symmetric, so $A^T=A$.
The transpose is kept in the general derivation:
the adjoint equation is governed by the adjoint operator.

Once $p_h$ has been computed, the reduced gradient is

$$
\nabla f_h(u)=\alpha M_u u + B^T p.
$$

Finite element counterpart of the continuous identity:

$$
\nabla f(u)=\alpha u+p
$$
when the control space is identified with $L^2(\Omega)$ through the Riesz map.

Computational pattern:

1. state solve;
2. adjoint solve;
3. gradient assembly.

The adjoint principle is unchanged by discretization.

---

## Discrete KKT System

Instead of reducing the problem to the control alone, we may keep state, control,
and adjoint unknowns together.
In the unconstrained case, the discrete first-order system is

$$
\begin{aligned}
A y - B u &= f,\\
A^T p - M_y y &= -M_y y_d,\\
\alpha M_u u + B^T p &= 0.
\end{aligned}
$$

Reordering the equations as before, this becomes the block system

$$
\begin{pmatrix}
M_y & 0 & -A^T\\
0 & \alpha M_u & B^T\\
A & -B & 0
\end{pmatrix}
\begin{pmatrix}
y\\
u\\
p
\end{pmatrix}
=
\begin{pmatrix}
M_y y_d\\
0\\
f
\end{pmatrix}.
$$

This is the discrete saddle-point system associated with the optimality conditions.

Structure:

- the $(1,1)$ block comes from state tracking;
- the $(2,2)$ block comes from Tikhonov regularization;
- the off-diagonal blocks encode state and adjoint couplings;
- the zero block in the lower-right corner is the saddle-point signature.

Discretization preserves the continuous KKT structure and yields a sparse block linear system.

This algebraic form motivates:

- sparse direct solvers;
- Krylov methods;
- block preconditioners;
- all-at-once implementations.

---

## Box Constraints in the Discrete Problem

Suppose now that the control is constrained componentwise:

$$
u_{\min}\le u \le u_{\max}.
$$

Then the reduced discrete problem is

$$
\min_{u\in U_{\mathrm{ad}}^h} f_h(u),
$$
where $U_{\mathrm{ad}}^h$ is a box in $\mathbb R^{n_u}$.

The first-order condition becomes the finite-dimensional variational inequality

$$
\nabla f_h(\bar u)^T(v-\bar u)\ge 0
\qquad \forall v\in U_{\mathrm{ad}}^h.
$$

Using the discrete gradient formula,

$$
(\alpha M_u \bar u+B^T\bar p)^T(v-\bar u)\ge 0
\qquad \forall v\in U_{\mathrm{ad}}^h.
$$

If the control basis is chosen so that the discrete control inner product is diagonal
or easily invertible, this leads to the discrete projection formula

$$
\bar u=\Pi_{U_{\mathrm{ad}}^h}\!\left(\bar u-\rho(\alpha M_u\bar u+B^T\bar p)\right).
$$

In the simplest lumped-mass or nodal setting, this reduces componentwise to

$$
\bar u_i=\min\!\left\{u_{\max,i},\max\{u_{\min,i},-\alpha^{-1}\bar p_i\}\right\}.
$$

The discrete constrained problem has the same geometry as the finite-dimensional
KKT picture from Lecture 2 and the continuous projection picture from Lecture 5.

---

## A Minimal 1D Finite Element Thought Experiment

To make the previous formulas concrete, imagine $\Omega=(0,1)$ and let $V_h$ be the space
of continuous piecewise linear functions on a uniform mesh.

Then:

- the stiffness matrix $A$ is tridiagonal;
- the mass matrix $M_y$ is sparse and symmetric positive definite;
- the discrete state solve means solving one sparse linear system;
- the discrete adjoint solve has the same matrix and a different right-hand side.

One reduced-gradient iteration has the form:

1. choose the current control vector $u_k$;
2. solve

   $$
   A y_k = B u_k + f;
   $$

3. solve

   $$
   A^T p_k = M_y(y_k-y_d);
   $$

4. compute

   $$
   g_k=\alpha M_u u_k+B^T p_k;
   $$

5. update the control with a finite-dimensional optimization step.

This is the finite element version of the state-adjoint-gradient loop of the
continuous reduced formulation.

---

## References

- Fredi Troeltzsch, *Optimal Control of Partial Differential Equations*, Chapters 3 and 4.
- Juan Carlos De Los Reyes, *Numerical PDE-Constrained Optimization*.
- Alfio Borzi and Volker Schulz, *Computational Optimization of Systems Governed by Partial Differential Equations*.

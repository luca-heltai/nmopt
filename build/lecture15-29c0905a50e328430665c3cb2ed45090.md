# One-Shot KKT and PDAS for Inverse Poisson Coefficient Identification

## Overview

In the previous lecture we discussed parameter estimation in elliptic
problems at the continuous level.  We now turn that theory into a concrete
all-at-once finite element algorithm for the identification of the diffusion
coefficient in Poisson's equation.

The goal of this lecture is twofold:

- derive the optimality system for a coefficient-identification problem with
  box constraints;
- explain in detail the deal.II implementation in
  [`codes/dealii/source/inverse_poisson_kkt.cc`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc),
  together with the class declaration
  [`codes/dealii/include/inverse_poisson_kkt.h`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/include/inverse_poisson_kkt.h),
  the driver
  [`codes/dealii/execs/inverse_poisson_kkt.cc`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/execs/inverse_poisson_kkt.cc),
  and the default parameter file
  [`codes/dealii/parameters/inverse_poisson_kkt.prm`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/parameters/inverse_poisson_kkt.prm).

The computational strategy is not reduced optimization.  Instead, we assemble
and solve the Karush-Kuhn-Tucker system for the state, adjoint, and control
together.  Because the control enters the operator, the KKT system is
nonlinear.  The code therefore combines:

- a one-shot KKT formulation;
- a Newton linearization of the nonlinear optimality system;
- a primal-dual active set strategy for the box constraints;
- `AffineConstraints<double>` to freeze active control unknowns directly in
  the linearized solve.

Throughout the lecture, let $\Omega\subset\mathbb R^d$ be a bounded polygonal
or polyhedral Lipschitz domain.

---

## The Continuous Inverse Problem

We want to reconstruct a diffusion coefficient $a$ from measurements of the
corresponding state $y$.  The forward model is

$$
\begin{cases}
-\nabla\cdot(a\nabla y)=f & \text{in }\Omega,\\
y=0 & \text{on }\partial\Omega.
\end{cases}
$$

The coefficient is the control variable.  We assume pointwise bounds

$$
a_a(x)\le a(x)\le a_b(x)
\qquad\text{a.e. in }\Omega,
$$

with

$$
0< a_a(x)\le \underline a(x) \le a_b(x)
$$

almost everywhere.  The positivity assumption is essential because it keeps
the state equation uniformly elliptic.

Given a desired state $y_d\in L^2(\Omega)$ and a reference coefficient
$a_{\mathrm{ref}}\in L^2(\Omega)$, we consider

$$
\min_{(y,a)}
J(y,a)
:=
\frac12\|y-y_d\|_{L^2(\Omega)}^2
+
\frac{\alpha}{2}\|a-a_{\mathrm{ref}}\|_{L^2(\Omega)}^2
$$

subject to the PDE and the box constraints, where $\alpha>0$.

The admissible set is

$$
U_{\mathrm{ad}}
:=
\{a\in L^\infty(\Omega): a_a\le a\le a_b \text{ a.e. in }\Omega\}.
$$

The natural weak state equation is: find $y\in H_0^1(\Omega)$ such that

$$
\int_\Omega a\nabla y\cdot \nabla v\,dx
=
\int_\Omega f v\,dx
\qquad \forall v\in H_0^1(\Omega).
$$

If $a\in U_{\mathrm{ad}}$, then

$$
b_a(y,v):=\int_\Omega a\nabla y\cdot\nabla v\,dx
$$

is bounded on $H_0^1(\Omega)\times H_0^1(\Omega)$ and coercive:

$$
b_a(v,v)\ge \underline a \|\nabla v\|_{L^2(\Omega)}^2.
$$

Hence Lax-Milgram gives a unique state $y(a)\in H_0^1(\Omega)$.

---

## Sensitivity with Respect to the Coefficient

Unlike distributed control, the control does not appear additively in the
right-hand side.  It multiplies the highest-order term.  This is the key
structural difference.

Let $a\in U_{\mathrm{ad}}$ and perturb it by $h$.  Formally differentiating
the weak state equation gives the sensitivity $z=S'(a)h\in H_0^1(\Omega)$:

$$
\int_\Omega a\nabla z\cdot \nabla v\,dx
=
-\int_\Omega h \nabla y(a)\cdot \nabla v\,dx
\qquad \forall v\in H_0^1(\Omega).
$$

Thus the derivative is obtained by solving another coercive elliptic problem
with right-hand side depending on the current state gradient.

This immediately shows why coefficient identification is more nonlinear than
source control:

- the state depends nonlinearly on the control;
- the sensitivity equation itself depends on the current state;
- the optimality system will contain products of state and adjoint gradients.

---

## Lagrangian and First-Order Conditions

Introduce the Lagrangian

$$
\mathcal L(y,p,a)
:=
\frac12\|y-y_d\|_{L^2(\Omega)}^2
+
\frac{\alpha}{2}\|a-a_{\mathrm{ref}}\|_{L^2(\Omega)}^2
+
\int_\Omega a\nabla y\cdot \nabla p\,dx
-
\int_\Omega f p\,dx,
$$

where $p\in H_0^1(\Omega)$ is the adjoint variable.

The derivative with respect to $p$ gives back the state equation:

$$
\int_\Omega a\nabla y\cdot \nabla v\,dx
=
\int_\Omega f v\,dx
\qquad \forall v\in H_0^1(\Omega).
$$

The derivative with respect to $y$ yields the adjoint equation:

$$
\int_\Omega a\nabla p\cdot \nabla v\,dx
=
-\int_\Omega (y-y_d)v\,dx
\qquad \forall v\in H_0^1(\Omega).
$$

The sign depends on the sign convention chosen in the Lagrangian.  The code
uses the equivalent convention

$$
\int_\Omega a\nabla p\cdot \nabla v\,dx
=
\int_\Omega (y_d-y)v\,dx,
$$

which is the same equation written with the right-hand side moved to the
other side.

For the control variable, the Gâteaux derivative is

$$
D_a\mathcal L(y,p,a)[h]
=
\alpha\int_\Omega (a-a_{\mathrm{ref}})h\,dx
+
\int_\Omega h \nabla y\cdot \nabla p\,dx.
$$

Therefore the gradient density is

$$
g(a)
=
\alpha(a-a_{\mathrm{ref}})
+
\nabla y\cdot \nabla p.
$$

The constrained first-order condition is the variational inequality

$$
\int_\Omega g(a)(\widetilde a-a)\,dx \ge 0
\qquad
\forall \widetilde a\in U_{\mathrm{ad}}.
$$

Equivalently, if the control space is interpreted with the $L^2$ inner
product, one has the projection formula

$$
a
=
P_{[a_a,a_b]}
\bigl(a-\gamma g(a)\bigr)
\qquad \forall \gamma>0.
$$

This is the starting point for active-set methods.

---

## Primal-Dual Active Set Interpretation

Define the multiplier-like quantity

$$
\lambda := g(a).
$$

Then the box-constrained optimality conditions can be written pointwise as

$$
\begin{cases}
\lambda = 0 & \text{where } a_a < a < a_b,\\
\lambda \ge 0 & \text{where } a=a_a,\\
\lambda \le 0 & \text{where } a=a_b.
\end{cases}
$$

In a primal-dual active set method, one predicts the active sets through a
shifted complementarity test.  Given a parameter $c>0$, define

$$
\mathcal A_-(a,\lambda)
:=
\{x\in\Omega:\lambda(x)+c(a(x)-a_a(x))<0\},
$$

$$
\mathcal A_+(a,\lambda)
:=
\{x\in\Omega:\lambda(x)+c(a(x)-a_b(x))>0\}.
$$

Then the inactive set is

$$
\mathcal I := \Omega\setminus (\mathcal A_-\cup \mathcal A_+).
$$

The active-set update is:

- on $\mathcal A_-$, impose $a=a_a$;
- on $\mathcal A_+$, impose $a=a_b$;
- on $\mathcal I$, solve the stationarity equation $\lambda=0$ together with
  state and adjoint equations.

For linear-quadratic box-constrained control this is equivalent to a
semismooth Newton method.  In our coefficient-identification problem, the
state equation is already nonlinear in the control, so PDAS is combined with
an additional Newton linearization of the KKT system.

---

## One-Shot KKT Formulation

Let the unknown be the triple

$$
X := (y,p,a).
$$

The optimality system can be written as a nonlinear operator equation

$$
F(X)=0,
$$

where the three blocks are:

1. state residual,
2. adjoint residual,
3. control stationarity residual.

The one-shot philosophy means that we do not eliminate the state with a
reduced map.  We solve directly for all variables at once.

This has two computational consequences:

- the assembled linear system is a coupled saddle-point-type Jacobian;
- each Newton step simultaneously updates state, adjoint, and coefficient.

The price is a larger linear algebra problem.  The benefit is that the code
reflects the structure of the continuous optimality system almost literally.

---

## Newton Linearization of the KKT System

Suppose $X^k=(y^k,p^k,a^k)$ is the current iterate.  Newton's method solves

$$
DF(X^k)\,\delta X^k = -F(X^k),
$$

and updates

$$
X^{k+1}=X^k+\delta X^k.
$$

Because the diffusion coefficient multiplies both $\nabla y$ and $\nabla p$,
the Jacobian contains the couplings

$$
\delta a\, \nabla y^k\cdot \nabla v,
\qquad
a^k \nabla \delta y\cdot \nabla v,
$$

for the state block, and similarly

$$
\delta a\, \nabla p^k\cdot \nabla v,
\qquad
a^k \nabla \delta p\cdot \nabla v
$$

for the adjoint block.

The derivative of the stationarity equation produces

$$
\alpha\,\delta a
+
\nabla \delta y\cdot \nabla p^k
+
\nabla y^k\cdot \nabla \delta p.
$$

This explains why the matrix assembled in the code is not the same as the
linear-quadratic KKT matrix from distributed Poisson control.

---

## Discretization Used in the Code

The code uses a single `FESystem` with three scalar components:

- state: continuous `FE_Q`;
- adjoint: continuous `FE_Q`;
- control: discontinuous `FE_DGQ`.

The choice

$$
Y_h \subset H_0^1(\Omega),
\qquad
P_h \subset H_0^1(\Omega),
\qquad
U_h \subset L^2(\Omega)
$$

is natural:

- $y$ and $p$ need $H^1$ conformity;
- the coefficient only appears as an $L^\infty$ or $L^2$ quantity in the
  variational equations;
- a discontinuous coefficient space makes pointwise box projection simple and
  local.

The implementation creates exactly this structure in
[`create_fe_system()` is implicit in the constructor and parameter parsing](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L45)
through
`FE_Q(state_degree)`, `FE_Q(state_degree)`, and `FE_DGQ(control_degree)`.

Homogeneous Dirichlet conditions for state and adjoint are imposed with
`VectorTools::interpolate_boundary_values(...)` on the corresponding
components in
[setup_system()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L186).

---

## Algorithm Structure

The overall solver is easier to understand from the following conceptual
diagram than from the exact file layout.

```{mermaid}
flowchart TD
    A[Choose mesh, FE spaces, data, bounds] --> B[Initialize coefficient a^0]
    B --> C[Assemble nonlinear KKT residual]
    C --> D[Compute stationarity lambda]
    D --> E[Predict lower and upper active sets]
    E --> F[Assemble KKT Jacobian]
    F --> G[Impose active controls with AffineConstraints]
    G --> H[Solve linearized KKT system for Newton correction]
    H --> I[Line search / damping on the one-shot update]
    I --> J[Project inactive controls back into the box]
    J --> K[Write VTU and append to PVD history]
    K --> L{Residual improved?}
    L -- yes --> M{Converged?}
    M -- no --> C
    M -- yes --> N[Write final state, adjoint, diffusion]
    L -- no --> O[Stop: no further residual decrease]
```

This is not a literal copy of the source code.  It is the mathematical logic
implemented by the solver.

---

## Mapping the Theory to the deal.II Class

The main class is
[`InversePoissonKKT`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/include/inverse_poisson_kkt.h).
Its public interface is intentionally small: the driver constructs the class,
reads the parameter file, and calls `run()`.

The most important internal methods are:

- [setup_system()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L186):
  distributes DoFs, renumbers by block, and builds the global constraint
  object for homogeneous Dirichlet data and hanging nodes.
- [initialize_control_mass_matrix()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L276):
  assembles the mass matrix on the control block and factorizes it.  This is
  used to transform the control residual into an $L^2$-type stationarity
  indicator.
- [assemble_system()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L366):
  assembles both the nonlinear KKT residual and its Jacobian at the current
  iterate.
- [compute_stationarity()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L579):
  extracts the control residual and applies the inverse control mass matrix to
  obtain the stationarity quantity used by PDAS.
- [compute_lower_active_set()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L595)
  and
  [compute_upper_active_set()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L618):
  implement the shifted active-set tests.
- [make_control_constraints()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L643):
  converts the active-set prediction into `AffineConstraints<double>` by
  freezing active control degrees of freedom to lower or upper bounds.
- [solve_linearized_system()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L676):
  condenses the matrix with the active constraints and solves the Newton
  system with `SparseDirectUMFPACK`.
- [run()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L833):
  orchestrates the full PDAS/Newton loop, damping, convergence logic, and
  output.

---

## The Residual Assembled by the Code

In
[assemble_system()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L366)
the code evaluates the current state, adjoint, and control at quadrature
points and assembles three residual contributions on each cell.

For a test basis function in the state component, it inserts

$$
\int_K a_h \nabla y_h\cdot \nabla v_h\,dx
-
\int_K f v_h\,dx.
$$

For a test basis function in the adjoint component, it inserts

$$
\int_K a_h \nabla p_h\cdot \nabla v_h\,dx
+
\int_K (y_h-y_d)v_h\,dx.
$$

For a test basis function in the control component, it inserts

$$
\int_K
\Bigl(
\alpha(a_h-a_{\mathrm{ref},h})
-
\nabla y_h\cdot \nabla p_h
\Bigr)
w_h\,dx.
$$

The sign in the last term reflects the sign convention chosen in the discrete
Lagrangian.  This is why the stationarity residual appears in the source as

```text
regularization * (coefficient - reference) - state_gradient * adjoint_gradient
```

rather than with a plus sign.

---

## Why `AffineConstraints<double>` Are a Good Fit for PDAS

The key idea of the implementation is that an active set is not handled by a
special projection after solving the linear system.  Instead, it is inserted
into the linear algebra structure itself.

If a control DoF is predicted active at the lower bound, the code adds the
constraint

$$
\delta a_i = a_{a,i}-a_i^k,
$$

and similarly at the upper bound

$$
\delta a_i = a_{b,i}-a_i^k.
$$

This means the Newton correction already respects the active-set decision.
The constrained linear system is built by:

1. copying the current Jacobian;
2. merging the PDE constraints with the active-set constraints;
3. condensing matrix and right-hand side;
4. solving for the correction;
5. redistributing the constrained values.

This is exactly the role of
[make_control_constraints()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L643)
and
[solve_linearized_system()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L676).

For teaching purposes, this is an excellent example of how a nonsmooth
optimization idea can be translated into standard finite element tools.

---

## Damping and Stopping

Because the KKT system is nonlinear in the coefficient, a full Newton step may
be too aggressive.  The code therefore performs a simple residual-based
damping strategy inside
[run()](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/inverse_poisson_kkt.cc#L833):

- compute the Newton update;
- try step lengths $1,1/2,1/4,\ldots$;
- accept the first step that decreases the KKT residual norm;
- otherwise keep the best trial found during backtracking.

This is not yet a polished globalization strategy, but it is mathematically
honest and easy to inspect.

The solver stops in one of two situations:

- the active set is unchanged and both the update norm and residual norm are
  below tolerance;
- the active set is unchanged and no damped Newton step improves the residual.

The second condition is important in a teaching code: it prevents endless
iterations when the current globalization strategy stalls.

---

## Build and Run

From the deal.II directory:

```bash
cd codes/dealii
cmake -S . -B build
cmake --build build --target inverse_poisson_kkt_2d.g
./build/inverse_poisson_kkt_2d.g parameters/inverse_poisson_kkt.prm
```

The default parameter file is written for a two-dimensional manufactured
solution.  In particular:

- the desired state is `sin(pi*x)*sin(pi*y)`;
- the exact diffusion coefficient is `1+x+y`;
- the forcing term is chosen consistently with that exact pair;
- the coefficient bounds are positive, ensuring ellipticity.

The output consists of:

- one `.vtu` file for each PDAS/Newton iteration;
- one `.pvd` file collecting the full nonlinear history;
- one final `.vtu` file with the last iterate.

The VTU fields include:

- `state`;
- `adjoint`;
- `diffusion`;
- projected desired state;
- reference and exact coefficient;
- lower and upper bounds;
- lower-active, upper-active, and inactive indicators.

---

## What to Learn from This Example

This laboratory shows a useful transition point in PDE-constrained
optimization.

For distributed linear-quadratic control, the KKT system is already linear.
For coefficient identification, three new ingredients appear simultaneously:

- nonlinearity through the control-dependent operator;
- nonsmoothness through box constraints;
- mixed regularity, since the coefficient naturally lives in a discontinuous
  finite element space.

The deal.II implementation shows that these difficulties can still be handled
within the same finite element framework, provided one is careful about:

- the variational structure of the residual;
- the Jacobian couplings;
- the choice of discrete spaces;
- the enforcement of active constraints.

In this sense, the code is not just an application of the previous lectures.
It is a prototype for much more advanced inverse problems and PDE-constrained
Newton methods.

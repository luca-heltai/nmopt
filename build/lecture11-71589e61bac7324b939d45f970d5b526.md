# Active Sets and Affine Constraints for Box-Constrained KKT Systems

## Overview

The previous lectures introduced:

- variational inequalities for pointwise control constraints;
- reduced and all-at-once formulations for elliptic optimal control;
- the discrete KKT system for the unconstrained linear-quadratic problem.

This lecture focuses directly on the next numerical question:

> how do we solve the discrete KKT system when the control is subject to box
> constraints?

The answer is a **primal-dual active set** strategy implemented through
`AffineConstraints<double>`.

The technical inspiration is the deal.II tutorial
[`step_41`](https://www.dealii.org/current/doxygen/deal.II/step_41.html),
but we will not discuss the obstacle problem itself.
We only reuse the algorithmic idea:

> once the active set is known, the inequality-constrained problem becomes a
> linear problem with additional affine equality constraints.

The logical path of the lecture is:

1. write the KKT system for the box-constrained control problem;
2. identify the complementarity structure on the control variable;
3. define the lower, upper, and inactive sets;
4. explain the PDAS iteration;
5. show why active control constraints are naturally represented by `AffineConstraints`;
6. connect this with the deal.II implementation in `kkt_box_constraints.cc`.

---

## The PDE-Constrained Problem

Consider the standard distributed control problem:

$$
\min_{(y,u)} \;
\frac12\|y-y_d\|_{L^2(\Omega)}^2
+\frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2
$$

subject to

$$
\mathcal A y = u + f
\qquad \text{in }\Omega,
$$

with suitable boundary conditions on the state, and with box constraints

$$
u_a(x)\le u(x)\le u_b(x)
\qquad \text{a.e. in }\Omega.
$$

Here:

- $y$ is the state;
- $u$ is the control;
- $y_d$ is the desired state;
- $\alpha>0$ is the regularization parameter.

Without the box constraints, the optimality system is linear.
With the box constraints, the control equation becomes a variational
inequality, or equivalently a complementarity system.

---

## Continuous Optimality Structure

Introduce the adjoint state $p$.
The unconstrained first-order condition would read

$$
\alpha u - p = 0.
$$

With box constraints, this is replaced by

$$
u \in U_{ad}
\qquad\text{and}\qquad
(\alpha u-p, v-u)_{L^2(\Omega)} \ge 0
\quad \forall v\in U_{ad},
$$

where

$$
U_{ad} = \{v\in L^2(\Omega): u_a\le v\le u_b \text{ a.e.}\}.
$$

Pointwise, this means:

- if $u_a(x) < u(x) < u_b(x)$, then $\alpha u(x)-p(x)=0$;
- if $u(x)=u_a(x)$, then $\alpha u(x)-p(x)\ge 0$;
- if $u(x)=u_b(x)$, then $\alpha u(x)-p(x)\le 0$.

So the bound-constrained control problem already contains an active-set
decomposition at the continuous level.

---

## Discrete KKT Viewpoint

After finite element discretization, we obtain a KKT system in the unknowns

$$
(Y,P,U),
$$

where:

- $Y$ is the vector of state coefficients;
- $P$ is the vector of adjoint coefficients;
- $U$ is the vector of control coefficients.

Ignoring the bounds for a moment, the assembled monolithic system has the form

$$
\mathcal K
\begin{pmatrix}
Y\\ P\\ U
\end{pmatrix}
=
\mathcal F.
$$

The box constraints act only on the control block.
Therefore, the nonlinearity is not in the PDE operator itself, but in the
fact that for each control degree of freedom we must decide whether:

- the lower bound is active;
- the upper bound is active;
- the degree of freedom is inactive.

This is exactly what an active-set method does.

---

## Lower, Upper, and Inactive Sets

Let $r_u(U)$ denote the discrete stationarity residual on the control block,
obtained by plugging the current iterate into the unconstrained KKT system.
This residual is a **dual** quantity.

To compare it with the primal coefficient gaps $U-U_a$ and $U-U_b$, we should
first map it back to the control space through the control mass matrix
$M_u$:

$$
g_u(U) := M_u^{-1} r_u(U).
$$

The vector $g_u(U)$ is the discrete $L^2$-gradient, i.e. the Riesz
representative of the control residual in the primal control space.

For every control degree of freedom $i$, we distinguish three cases:

- lower active:
  $$
  U_i = (U_a)_i;
  $$
- upper active:
  $$
  U_i = (U_b)_i;
  $$
- inactive:
  the unconstrained stationarity equation is enforced.

This leads to three sets:

$$
\mathcal A_-,
\qquad
\mathcal A_+,
\qquad
\mathcal I.
$$

They form a partition of the control indices.

---

## PDAS Update Rule

A primal-dual active set strategy updates these sets using both primal and dual
information.
One convenient discrete criterion is:

$$
\mathcal A_-^{k+1}
=
\left\{
i:\;
-g_{u,i}(U^k) + c\bigl(U_i^k-(U_a)_i\bigr) < 0
\right\},
$$

$$
\mathcal A_+^{k+1}
=
\left\{
i:\;
-g_{u,i}(U^k) + c\bigl(U_i^k-(U_b)_i\bigr) > 0
\right\},
$$

with a parameter $c>0$.
This scaling is important: if we used the raw algebraic residual $r_u$ in the
test, then the parameter $c$ would inherit mesh-dependent scaling from the
dual basis.
Using $M_u^{-1}r_u$ instead makes the comparison with the primal bound gaps
consistent.
The remaining indices define the inactive set:

$$
\mathcal I^{k+1}
=
\{1,\dots,N_u\}\setminus(\mathcal A_-^{k+1}\cup\mathcal A_+^{k+1}).
$$

Given these sets, the next iterate is computed by solving the KKT system with

$$
U_i = (U_a)_i \qquad \forall i\in\mathcal A_-^{k+1},
$$

$$
U_i = (U_b)_i \qquad \forall i\in\mathcal A_+^{k+1},
$$

and with the original control optimality equation retained only on
$\mathcal I^{k+1}$.

This is the core idea:

1. predict which bounds are active;
2. freeze those DoFs to the corresponding bounds;
3. solve the resulting linear system;
4. repeat until the active sets stop changing.

---

## Why This Becomes a Linear Problem

The original box-constrained KKT system is nonlinear because we do not know in
advance which indices satisfy:

- $U_i=(U_a)_i$;
- $U_i=(U_b)_i$;
- the unconstrained stationarity equation.

But once $\mathcal A_-$ and $\mathcal A_+$ are fixed, that ambiguity
disappears.
The active inequalities are replaced by explicit equalities:

$$
U_i = (U_a)_i
\quad\text{or}\quad
U_i = (U_b)_i.
$$

At that point, the problem is again a linear KKT system with additional affine
constraints.

This is the precise sense in which PDAS converts an inequality-constrained
problem into a sequence of linear solves.

---

## Why `AffineConstraints<double>` Is the Right Tool

In deal.II, a condition of the form

$$
U_i = \text{given value}
$$

is exactly an affine constraint.
So an active control index should be handled in the same language as:

- Dirichlet boundary values;
- hanging-node constraints;
- any other prescribed linear relation among DoFs.

For a lower-active control DoF, we add the line

$$
U_i = (U_a)_i.
$$

For an upper-active control DoF, we add the line

$$
U_i = (U_b)_i.
$$

The resulting `AffineConstraints<double>` object is then merged with the
existing constraints of the KKT problem.

This gives a very clean implementation pattern:

1. start from the usual constraints already present in the discretization;
2. add the active control lines;
3. condense matrix and right-hand side with the combined constraint object;
4. solve;
5. distribute the constrained values back to the solution vector.

This is mathematically cleaner than editing rows of the matrix by hand, and it
fits the normal deal.II workflow.

---

## Translation into the Course Code

The implementation in
[`codes/dealii/execs/kkt_box_constraints.cc`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/execs/kkt_box_constraints.cc)
follows exactly this structure.

### Step 1: Build the unconstrained KKT system

The class [`KKT`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/include/kkt.h)
assembles the monolithic state-adjoint-control matrix and the corresponding
right-hand side.

### Step 2: Compute the control residual and scale it correctly

Given a current iterate, the code computes the residual of the full KKT system
and extracts the control block.
Then it applies the inverse control mass matrix to obtain the discrete
$L^2$-gradient on the control space.
This scaled quantity is what enters the PDAS test.

### Step 3: Detect the active sets

Using the lower and upper bound functions, the code marks:

- the indices where the lower bound is active;
- the indices where the upper bound is active.

### Step 4: Build control constraints

For every active control DoF, the code adds an inhomogeneous line to an
`AffineConstraints<double>` object:

- lower active $\Rightarrow$ constrain to the lower bound value;
- upper active $\Rightarrow$ constrain to the upper bound value.

Then it merges these lines with the constraints already provided by the KKT
discretization.

### Step 5: Condense and solve

The combined constraint object is used to:

- condense the block matrix and right-hand side;
- solve the resulting linear system;
- distribute the constrained values to the full solution vector.

This is the direct transcription of the PDAS idea into deal.II.

---

## Why This Matters

The main conceptual gain is that the box-constrained control problem is no
longer treated as a special exception.
Instead, it is handled as:

- a KKT system;
- plus a changing set of affine constraints on the control block.

This viewpoint is useful because it scales naturally to:

- semismooth Newton interpretations of PDAS;
- more complicated complementarity systems;
- situations in which the monolithic structure matters more than the reduced one.

It also clarifies the difference with projected-gradient methods:

- in the reduced formulation, the constraint appears through a projection step;
- in the KKT formulation, the constraint appears through active equalities on a subset of control DoFs.

---

## Relation with Reduced Methods

The reduced and monolithic approaches solve the same constrained problem, but
they organize the numerics differently.

### Reduced formulation

- eliminate the state equation;
- minimize only in the control variable;
- project each iterate onto the admissible interval.

This is often the cleanest entry point if the emphasis is on optimization.

### Monolithic KKT formulation

- keep state, adjoint, and control together;
- read the bounds as complementarity conditions on the control block;
- use PDAS to decide which control DoFs are frozen to the bounds.

This is often the natural entry point if the emphasis is on:

- saddle-point systems;
- primal-dual variables;
- semismooth Newton and active-set methods.

---

## Files for This Lecture

- Technical reference for the active-set strategy:
  [`step_41`](https://www.dealii.org/current/doxygen/deal.II/step_41.html)
- Existing unconstrained KKT code:
  [`codes/dealii/include/kkt.h`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/include/kkt.h)
  and
  [`codes/dealii/source/kkt.cc`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/source/kkt.cc)
- Box-constrained executable:
  [`codes/dealii/execs/kkt_box_constraints.cc`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/execs/kkt_box_constraints.cc)
- Parameter file:
  [`codes/dealii/parameters/kkt_box_constraints.prm`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/parameters/kkt_box_constraints.prm)

---

## Suggested In-Class Discussion

1. Why does fixing the active sets turn the box-constrained KKT problem into a linear one?
2. Why are active control bounds naturally represented as affine constraints?
3. Which rows of the KKT system remain active on the inactive set?
4. What changes if the control is continuous (`FE_Q`) instead of discontinuous (`FE_DGQ`)?
5. Why must the PDAS test use the mass-scaled control residual instead of the raw algebraic residual?
6. What is gained by using `AffineConstraints<double>` instead of manual row editing?

The practical message of the lecture is therefore:

> for box constraints on the control, PDAS is implemented by turning the
> active bounds into inhomogeneous affine constraints on the control block of
> the KKT system.

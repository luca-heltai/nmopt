# Lecture 1: Introduction & Motivation (PDE-Constrained Optimal Control)

Numerical Methods for Optimal Control (NMOPT)

----

## What Is Optimal Control?

- Choose a control $u$ to minimize a cost $J(y,u)$.
- The state $y$ is constrained by a model (often a PDE).
- Optimization + differential constraints = PDE-constrained optimization.

----

## Abstract Formulation

Minimize
$$
\min_{y,u} \; J(y,u)
$$
subject to the state equation
$$
\mathcal{E}(y,u) = 0.
$$

Typical pieces:
- state equation (PDE)
- admissible controls $u \in \mathcal{U}_{\mathrm{ad}}$
- constraints (box constraints, etc.)

----

## Forward vs Optimal Control

- Forward (direct) problem:
  - input data given
  - solve PDE once to compute $y$
- Optimal control:
  - control $u$ unknown
  - PDE is solved repeatedly while optimizing $J$

Key point: the PDE is a *constraint*, not the objective.

----

## Reduced Formulation (The Big Idea)

If the state equation defines a control-to-state map $S$:
$$
y = S(u),
$$
then define the reduced cost
$$
f(u) := J(S(u),u),
$$
and solve
$$
\min_{u} f(u).
$$

Why it matters:
- separates "PDE solve" from "optimization step"
- makes gradients/adjoins systematic

---

## Tiny Example Notebook

- Finite-dimensional LQ reduction + a 1D slice plot.
- Generates figures used in these slides.

File: `codes/lecture01/finite_dim_analogy.ipynb`

----

## Example: Reduced Cost Shape

![Reduced cost in 1D](assets/lecture01/reduced_cost_1d.png)

----

## Control-to-State Pipeline

![Control-to-state diagram](assets/lecture01/control_to_state.png)

----

## A Finite-Dimensional Analogy

Given $A y = B u$ with invertible $A$:
$$
y = A^{-1} B u =: S u
$$
Reduced problem:
$$
\min_u f(u) := J(Su,u).
$$

This is the template we will reuse in infinite dimensions.

----

## Why PDEs Show Up

States: temperature, displacement, velocity, concentration, ...

Controls:
- distributed sources
- boundary inputs
- initial conditions
- parameters / coefficients

----

## Typical Examples

- Elliptic: stationary heating / diffusion
- Parabolic: time-dependent heating
- Nonlinear: Navier-Stokes flow control
- Inverse problems: parameter estimation / data assimilation

----

## Constraints

Common: box constraints on the control
$$
u_{\min} \le u \le u_{\max}.
$$

Leads to:
- variational inequalities
- KKT conditions
- active-set / projection-type algorithms

----

## Course Trajectory

Continuous theory -> discretization -> algorithms -> implementation -> applications

Next up:
- adjoints and gradients
- elliptic linear-quadratic control as a core model problem

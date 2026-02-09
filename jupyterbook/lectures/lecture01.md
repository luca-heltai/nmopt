# Lecture 1 — Introduction and Motivation for PDE-Constrained Optimal Control

## Overview

This lecture introduces the general framework of **optimal control problems (OCPs)** with a particular focus on problems governed by **partial differential equations (PDEs)**.  
The goal is to clarify *what* an optimal control problem is, *why* PDEs arise naturally in this context, and *how* control and optimization interact.

## Structure of the Lecture

- Motivation and real-world examples
- General formulation of optimal control problems
- Finite-dimensional analogy and reduced formulation
- From ODEs to PDEs: infinite-dimensional issues
- Types of controls, observations, and constraints
- Discussion, perspective, and roadmap of the course

---

## 1. What Is an Optimal Control Problem?

An **optimal control problem** consists of determining a control variable  
$$
u
$$
such that a given **cost functional**
$$
J(y,u)
$$
is minimized (or maximized), subject to the constraint that the **state variable**
$$
y
$$
solves a governing equation.

### Abstract structure

- **State equation**
  $$
  \mathcal{E}(y,u) = 0
  $$
  (ODE, PDE, or algebraic system)

- **Control variable**
  $$
  u \in \mathcal{U}_{\mathrm{ad}}
  $$

- **Cost functional**
  $$
  J(y,u)
  $$

- **Constraints**
  - control constraints: $u \in \mathcal{U}_{\mathrm{ad}}$
  - possibly state constraints: $y \in \mathcal{Y}_{\mathrm{ad}}$

An optimal control problem can be seen as an **optimization problem with constraints**, where the constraints are given by a differential equation.

---

## 2. Forward Problems vs Optimal Control Problems

### Forward (direct) problem

- Input data (coefficients, sources, boundary conditions) are given.
- The PDE is solved once to compute the state $y$.

### Optimal control problem

- Some input data (the control $u$) are *unknown*.
- The PDE solution depends on $u$: $y = y(u)$.
- We seek the control $u$ that optimizes a performance criterion.

**Key difference**:  
> In optimal control, the PDE is part of the constraint, not the objective.

---

## 2.1 A Finite-Dimensional Analogy

To understand the structure of optimal control problems without technical overhead,  
we begin with a finite-dimensional constrained optimization problem.

Let
$$
J(y,u) \in \mathbb{R}, \quad y \in \mathbb{R}^n, \quad u \in \mathbb{R}^m,
$$
and consider
$$
\min_{y,u} J(y,u) \quad \text{subject to} \quad A y = B u,
$$
where $A \in \mathbb{R}^{n \times n}$ is invertible.

Solving the constraint yields
$$
y = A^{-1} B u =: S u,
$$
where $S$ is the **control-to-state operator**.

Substituting into the cost functional leads to the **reduced problem**
$$
\min_{u} f(u) := J(Su,u).
$$

This reduction is the conceptual blueprint for PDE-constrained optimization.

---

## 3. Why PDEs?

Many physical, biological, and engineering systems are naturally modeled by PDEs:

- heat conduction and diffusion;
- fluid dynamics;
- elasticity and structural mechanics;
- electromagnetism;
- reaction–diffusion systems;
- transport and advection.

In these contexts:

- the **state** represents a physical quantity (temperature, velocity, concentration);
- the **control** represents an external action (source term, boundary input, coefficient).

---

## 3.1 From Finite to Infinite Dimensions

When passing from finite-dimensional systems to PDEs, several new mathematical  
issues arise:

- the state and control live in infinite-dimensional function spaces;
- bounded sets are no longer compact;
- existence of minimizers is nontrivial;
- derivatives must be understood in Banach or Hilbert spaces.

As a result, tools from functional analysis become essential.  
Throughout the course, we will primarily work in Hilbert spaces,  
where inner products allow a natural generalization of gradients and adjoints.

---

## 4. Typical Examples

### 4.1 Optimal Heating (Elliptic PDE)

- **State**: temperature $y(x)$
- **PDE**: stationary heat equation
- **Control**: heat source or boundary temperature
- **Cost**: tracking a desired temperature distribution

This is the prototypical **linear–quadratic elliptic control problem**.

---

### 4.2 Time-Dependent Heating (Parabolic PDE)

- **State**: temperature $y(x,t)$
- **PDE**: heat equation
- **Control**: time-dependent source or boundary input
- **Cost**: tracking at final time or over time

Introduces:

- adjoint equations backward in time;
- space–time discretization issues.

---

### 4.3 Flow Control (Navier–Stokes)

- **State**: velocity and pressure fields
- **PDE**: Navier–Stokes equations
- **Control**: body force or boundary velocity
- **Cost**: drag reduction, vorticity minimization

Leads to **nonlinear PDE-constrained optimization**.

---

### 4.4 Parameter Estimation and Data Assimilation

- **Control**: unknown parameters or initial conditions
- **State**: PDE solution
- **Cost**: mismatch between model output and observations

These problems are common in:

- meteorology;
- inverse problems;
- uncertainty quantification.

---

## 4.5 Classes of Control and Observation

Depending on how the control enters the PDE, we distinguish between:

- **distributed controls** (source terms in the domain);
- **boundary controls** (Dirichlet, Neumann, Robin);
- **initial controls** (for time-dependent problems);
- **parameter controls** (coefficients in the PDE).

Similarly, observations may be:

- distributed in the domain;
- restricted to subdomains;
- located on the boundary;
- taken at final time or over time intervals.

---

## 5. Control vs Controllability

- **Controllability**: can we drive the system exactly to a desired state?
- **Optimal control**: can we *approximately* reach a target in an optimal way?

For many PDEs:

- exact controllability may fail;
- optimal control remains well-posed and meaningful.

---

## 6. Components of a PDE-Constrained OCP

A PDE-constrained optimal control problem is characterized by:

1. **State equation**
2. **Control space**
3. **Observation operator**
4. **Cost functional**
5. **Constraints**

Graphically:

```
u  ──►  PDE  ──►  y  ──►  observation
 \__________________________/
              |
           cost
```

---

## 6.1 Control Constraints and Modeling Considerations

In realistic applications, controls are subject to constraints, such as:
$$
u_{\min} \le u \le u_{\max}.
$$

These constraints model physical, technological, or safety limitations and  
lead to **variational inequalities** and **Karush–Kuhn–Tucker (KKT) conditions**  
in the optimality system.

Control constraints play a crucial role both in the analysis and in the  
design of numerical algorithms.

---

## 7. Outlook: What Comes Next

In the next lectures we will:

- introduce the **control-to-state map** $u \mapsto y(u)$;
- derive **first-order optimality conditions** using adjoint equations;
- study **linear–quadratic elliptic problems** in detail;
- move progressively toward numerical approximation and implementation.

---

## References for This Lecture

Suggested reading:

- F. Tröltzsch, *Optimal Control of Partial Differential Equations*, Chapter 1  
- A. Manzoni, A. Quarteroni, S. Salsa, *Optimal Control of PDEs*, Chapter 1  
- J. C. De los Reyes, *Numerical PDE-Constrained Optimization*, Section 1

---

## Exercises and Discussion

1. Write down a forward PDE problem you are familiar with and identify which  
   quantities could realistically act as controls.  
2. Explain why the reduced formulation is computationally preferable to a  
   naive optimization over both state and control.  
3. Discuss the role of control constraints in modeling real systems.  
4. Which difficulties do you expect when moving from elliptic to parabolic  
   optimal control problems?

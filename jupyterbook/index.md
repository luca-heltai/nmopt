```{include} ../README.md
:start-after: Master’s Degree Course
```

Use the table of contents to navigate the 20 sessions (each ~2h).

```{admonition} Course at a Glance
:class: tip

- Format: 20 sessions (~2h each), lectures + laboratories.
- Structure: 3 blocks (Foundations -> Numerics -> Nonlinear/advanced topics).
- Quick links: [Lecture 01](lectures/lecture01.md), [Slides 01](/slideshow/slides01.html).
```

Below is a concise, structured summary of the course topics, following the three-block / 20-lecture split. This is suitable for the Jupyter Book landing page, syllabus summary, or course webpage.

## Course Topic Summary

The course is organized into three main blocks, progressing from foundational concepts to advanced numerical methods and applications. The emphasis is on PDE-constrained optimal control, with systematic integration of theory and numerical implementation.

### Block I: Foundations and Linear Elliptic Optimal Control (Lectures 1-6)

This block introduces the mathematical framework of optimal control problems governed by PDEs, with a focus on linear-quadratic elliptic problems.

#### Topics

- Introduction to optimal control problems: motivation, examples, and applications.
- Abstract formulation of PDE-constrained optimization problems.
- Control-to-state mapping and reduced cost functional.
- Adjoint equations and reduced gradients.
- Linear-quadratic optimal control of elliptic PDEs.
- Distributed and boundary controls.
- Existence and uniqueness of optimal controls.
- Control constraints and Karush-Kuhn-Tucker (KKT) conditions.
- Conceptual introduction to numerical discretization strategies: optimize-then-discretize vs discretize-then-optimize.

By the end of this block, students understand the continuous optimality system and its analytical properties.

### Block II: Numerical Approximation and Parabolic Optimal Control (Lectures 7-12)

This block focuses on the numerical solution of optimal control problems and extends the theory to time-dependent (parabolic) PDEs.

#### Topics

- Finite element discretization of elliptic optimal control problems.
- Reduced vs all-at-once formulations.
- Saddle-point systems and block structure.
- Iterative solvers and preconditioning strategies.
- Numerical treatment of control constraints.
- Linear-quadratic optimal control of parabolic PDEs.
- Adjoint equations backward in time.
- Space-time discretization strategies.
- Numerical algorithms for time-dependent optimal control problems.
- Computational aspects and scalability considerations.

Laboratories in this block introduce deal.II-based implementations for elliptic and parabolic problems.

### Block III: Nonlinear Problems and Advanced Topics (Lectures 13-20)

The final block addresses nonlinear and nonsmooth optimal control problems, together with advanced numerical algorithms and applications.

#### Topics

- Optimal control problems governed by semilinear PDEs.
- Optimization in Banach spaces.
- First-order optimality conditions for nonlinear problems.
- Second-order optimality conditions and stability.
- Newton and Sequential Quadratic Programming (SQP) methods.
- Primal-Dual Active Set (PDAS) methods and semismooth Newton methods.
- Nonsmooth problems: box constraints, sparsity-promoting costs.
- Selected advanced applications (e.g., flow control, parameter identification).
- Discussion of current research directions in PDE-constrained optimization.

This block prepares students to read current research literature and develop research-level numerical codes.

### Overall Learning Trajectory

The course progresses along the following axis:

continuous theory -> discrete formulations -> algorithms -> implementation -> applications

By the end of the course, students have:

- A solid theoretical understanding of PDE-constrained optimal control.
- Practical experience with finite element implementations in deal.II.
- The ability to critically assess numerical methods and modeling choices.

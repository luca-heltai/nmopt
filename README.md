# Numerical Methods for Optimal Control (NMOPT)

| Item | Link |
| --- | --- |
| Jupyter Book Status | [![Jupyter Book Status](https://github.com/luca-heltai/nmopt/actions/workflows/build-book.yml/badge.svg)](https://github.com/luca-heltai/nmopt/actions/workflows/build-book.yml) |
| Rendered book | <https://luca-heltai.github.io/nmopt/> |
| Course repository | <https://github.com/luca-heltai/nmopt> |
| Author | [Luca Heltai](https://github.com/luca-heltai) |

Master’s Degree Course  
**Numerical Methods for Optimal Control**  
Teacher: Luca Heltai (University of Pisa)
Total workload: 40 hours (lectures + laboratories)  
Main numerical library: **deal.II**

This repository contains the teaching material for the course, organized as a **Jupyter Book**.  
It includes lecture notes, theoretical material, numerical laboratories, and references.

---

## Course Objectives

The course introduces the mathematical theory and numerical methods for **PDE-constrained optimal control problems**.  
The focus is on:

- rigorous formulation of optimal control problems governed by PDEs;
- derivation of first-order (and selected second-order) optimality conditions;
- finite element discretization of optimal control problems;
- implementation of numerical algorithms using **deal.II**.

At the end of the course, students will be able to develop and analyze research-level prototype codes for PDE-constrained optimization.

---

## Prerequisites

- Advanced calculus
- Partial differential equations (elliptic and parabolic)
- Numerical methods for PDEs (finite element methods)
- Programming in **C++** and **Python** is helpful for the laboratories, but not required for the lectures
- Basic familiarity with Linux and scientific software tools.

Previous experience with **deal.II** is helpful but not mandatory.

---

## Numerical Laboratories (deal.II, Python)

The laboratories are an integral part of the course and include:

- distributed control of the Poisson equation;
- state–adjoint formulation and reduced gradient methods;
- box constraints on the control (projection methods, PDAS);
- time-dependent optimal control problems;
- analysis of computational performance and scalability.

The codes:

- are written in **C++** and **Python**;
- are based on **deal.II**
- are designed as starting points for extensions and final projects
- will be provided with detailed documentation and comments for you to understand and modify.

---

## References

- F. Tröltzsch, *Optimal Control of Partial Differential Equations*, AMS, 2010  
- A. Manzoni, A. Quarteroni, S. Salsa, *Optimal Control of Partial Differential Equations*, Springer, 2021  
- J. C. De los Reyes, *Numerical PDE-Constrained Optimization*, Springer, 2015  

Additional references are provided within the individual chapters.

---

## License

The content of this repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

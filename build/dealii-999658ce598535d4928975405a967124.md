---
title: deal.II laboratories
---

# deal.II laboratories

This section documents the C++/**deal.II** laboratory infrastructure used in the course for finite element experiments on PDE-constrained optimal control problems.

The source tree lives in [`codes/dealii`](https://github.com/luca-heltai/nmopt/tree/main/codes/dealii). It is organized as a course testbench: a small, inspectable code base that students can compile, modify, test, and extend during the laboratories.

## Purpose

The deal.II laboratory track complements the Python notebooks in [`jupyterbook/codes`](https://github.com/luca-heltai/nmopt/tree/main/jupyterbook/codes).

- The notebooks are used for fast experimentation and visualization.
- The deal.II code is used for realistic finite element implementations.
- The goal is to let students move from reduced models to research-style PDE control prototypes.

## Testbench structure

The directory is meant to contain:

- reference finite element applications for state, adjoint, and coupled optimality systems;
- a `CMake` build system with dimension-dependent executables;
- regression tests in deal.II style under [`codes/dealii/tests`](https://github.com/luca-heltai/nmopt/tree/main/codes/dealii/tests);
- unit tests under [`codes/dealii/gtests`](https://github.com/luca-heltai/nmopt/tree/main/codes/dealii/gtests);
- a Doxygen configuration under [`codes/dealii/doc`](https://github.com/luca-heltai/nmopt/tree/main/codes/dealii/doc).

## Planned evolution

The testbench hosts laboratories on:

- distributed control of elliptic equations;
- state, adjoint, and gradient assembly;
- reduced cost functionals and optimization loops;
- box constraints and active-set ideas;
- solver design, preconditioning, and performance comparisons;
- verification with tests that are simple enough to be read and modified by students.

## Working model for the laboratories

Students should be able to:

1. configure and compile the code;
2. run a reference state or optimality-system problem;
3. inspect the finite element assembly;
4. extend the code to include controls, constraints, and optimization algorithms;
5. validate changes through tests and documented experiments.

## Repository guide

For a directory-level overview and build instructions, see [`codes/dealii/README.md`](https://github.com/luca-heltai/nmopt/blob/main/codes/dealii/README.md).

The generated API documentation is available in the companion page [API reference](dealii-api.md).

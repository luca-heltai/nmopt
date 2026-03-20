deal.II Testbench for NMOPT Laboratories
========================================

| Doxygen | Indentation | Tests |
|----------|-------------|-------|
| [![Jupyter Book Status](https://github.com/luca-heltai/nmopt/actions/workflows/build-book.yml/badge.svg "Doxygen")](https://github.com/luca-heltai/nmopt/actions/workflows/build-book.yml) | [![Indentation Check](https://github.com/luca-heltai/nmopt/actions/workflows/dealii-indentation.yml/badge.svg "Indentation Check")](https://github.com/luca-heltai/nmopt/actions/workflows/dealii-indentation.yml) | [![deal.II CI](https://github.com/luca-heltai/nmopt/actions/workflows/dealii-tests.yml/badge.svg "Tests")](https://github.com/luca-heltai/nmopt/actions/workflows/dealii-tests.yml) |

This directory contains the C++/**deal.II** laboratory testbench for the course
*Numerical Methods for Optimal Control*.

It is organized as a compact, teaching-oriented code base for experiments on
PDE-constrained optimization. The goal is not to provide a large framework, but
a small and readable project that students can compile, inspect, modify, and
test during the laboratories.

Planned use in the course
-------------------------

This directory hosts laboratories on:

- elliptic state equations;
- adjoint equations and optimality systems;
- distributed and boundary control problems;
- reduced formulations and all-at-once formulations;
- constrained optimization algorithms;
- testing and reproducible numerical experiments.

Students should be able to use this directory as a testbench to:

- start from a working finite element code;
- introduce controls, adjoints, and cost functionals;
- compare alternative formulations and solvers;
- validate implementations with automated tests.

Directory structure
-------------------

The structure of the directory is the following:

 ./source
 ./execs
 ./include
 ./tests
 ./gtests
 ./doc

The directories are used to organize:

- application code for state, adjoint, and optimization modules;
- executable drivers in `./execs`;
- a test directory that uses deal.II-style regression tests;
- a test directory that uses GoogleTest;
- a documentation directory with a `Doxyfile`.

Build and test
--------------

The `CMakeLists.txt` generates executables and shared libraries from the files
in `./source`, `./include`, and `./execs`. The libraries are linked to the
tests, so that the same application code can be exercised through automated
checks.

After configuring and compiling the application, you can run

```sh
make test
```

or

```sh
ctest
```

to start the testsuite.

Documentation and course integration
------------------------------------

The course-level documentation for this directory is mirrored in the Jupyter
Book under the section **deal.II laboratories**. That section explains the role
of this testbench within the course and will track the progressive addition of
optimal control examples.

For general references on testing:

- deal.II testsuite:
  <https://www.dealii.org/developer/developers/testsuite.html>
- GoogleTest primer:
  <https://github.com/google/googletest/blob/master/googletest/docs/primer.md>

License
=======

See the file [LICENSE.md](./LICENSE.md) for details

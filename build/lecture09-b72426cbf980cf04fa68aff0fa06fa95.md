# deal.II Box-Constrained Optimal Control

## Laboratory Goal

In this laboratory we compare two numerical viewpoints for the elliptic
optimal control problem with pointwise box constraints

$$
u_a(x) \le u(x) \le u_b(x),
$$

cost functional

$$
J(y,u)=\frac12\|y-y_d\|_{L^2(\Omega)}^2
+\frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2,
$$

and state equation

$$
-\Delta y = u + f
\qquad \text{in }\Omega,
\qquad
y = 0
\qquad \text{on }\partial\Omega.
$$

The two formulations we want to experiment with are:

- **reduced formulation**: eliminate the state and optimize only with respect to the control by a projected-gradient step;
- **all-at-once formulation**: solve the discrete KKT system with a primal-dual active set (PDAS) strategy.

The point of the laboratory is not only to compute a solution, but to see how
the same constrained problem changes when viewed through:

- a reduced functional and a projection step;
- a monolithic KKT system and an active-set update;
- two different code organizations in deal.II.

This lecture should now be read together with
[Lecture 11](/Users/heltai/latex/courses/2025/02_nmopt/jupyterbook/lectures/lecture11.md),
where the active-set viewpoint is reformulated directly in terms of box
constraints on the control and `AffineConstraints<double>`.

---

## Relevant Files

The two executables used in this laboratory are:

- [`codes/dealii/execs/laplacian_box_constraints.cc`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/execs/laplacian_box_constraints.cc)
  for the reduced formulation;
- [`codes/dealii/execs/kkt_box_constraints.cc`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/execs/kkt_box_constraints.cc)
  for the all-at-once KKT formulation.

The associated parameter files are:

- [`codes/dealii/parameters/laplacian_box_constraints.prm`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/parameters/laplacian_box_constraints.prm)
- [`codes/dealii/parameters/kkt_box_constraints.prm`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/parameters/kkt_box_constraints.prm)

The underlying unconstrained KKT assembly is provided by:

- [`codes/dealii/include/kkt.h`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/include/kkt.h)
- [`codes/dealii/source/kkt.cc`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/source/kkt.cc)

---

## Part 1: Reduced Formulation

In the reduced viewpoint, the state is eliminated through the control-to-state
map.
The algorithmic loop is:

1. fix a control $u^k$;
2. solve the state equation for $y(u^k)$;
3. solve the adjoint equation for $p(u^k)$;
4. build the reduced gradient $\nabla \hat J(u^k)$;
5. take a descent step and project back onto the admissible interval.

The update has the form

$$
u^{k+1} = P_{[u_a,u_b]}\bigl(u^k-\tau_k\nabla\hat J(u^k)\bigr).
$$

This is the natural discretization of the variational inequality for the
reduced problem:

$$
(\nabla \hat J(\bar u), v-\bar u) \ge 0
\qquad \forall v\in U_{ad}.
$$

### What to observe

- the correct residual is not the unconstrained gradient, but the **projected-gradient residual**;
- once a control degree of freedom hits the bounds, an unconstrained descent direction is no longer admissible;
- the active/inactive information appears as a consequence of the projection, not as an explicit constraint object;
- the `.vtu` output stores the fields `lower_active`, `upper_active`, `inactive`.

### Suggested experiments

1. Reduce `Regularization` from `1e-4` to `1e-6` and observe how the control saturates.
2. Enlarge the bounds from `[-0.2,0.2]` to `[-1,1]` and compare with the nearly unconstrained case.
3. Change `Target state` from a step to `sin(pi*x)` and see whether the active region shrinks.
4. Increase `Global refinements` and check whether the active sets stabilize.
5. Modify the `Initial control` and compare the number of projected-gradient iterations.

---

## Part 2: All-at-Once KKT

In the monolithic formulation we keep state, adjoint, and control together.
Without box constraints, the discrete KKT system is linear.
With box constraints, the control optimality equation becomes a
complementarity system.

The PDAS strategy handles this by predicting:

- a lower active set $A_-$, where $u=u_a$;
- an upper active set $A_+$, where $u=u_b$;
- an inactive set, where the unconstrained stationarity condition is enforced.

At each PDAS iteration we therefore:

1. compute the stationarity indicator on the control block;
2. update $A_-$ and $A_+$;
3. freeze the active control DoFs;
4. solve the KKT system with those control equalities imposed;
5. repeat until the active sets no longer change.

This is the same structural idea explained in Lecture 11:
the nonlinear inequality problem is replaced by a sequence of linear solves in
which the active inequalities are turned into affine equalities.

In the present code, these equalities are imposed through
`AffineConstraints<double>` on the **control block** of the KKT system.
That is the correct deal.II interpretation of

$$
u_i=(u_a)_i
\qquad\text{or}\qquad
u_i=(u_b)_i.
$$

### What to observe

- the method often converges in a small number of active-set iterations;
- active control DoFs behave like Dirichlet conditions on the control variable;
- the code is closer to complementarity systems and semismooth Newton ideas than the reduced projected-gradient method;
- the use of `AffineConstraints` makes the implementation conceptually aligned with other deal.II constraints.

### Suggested experiments

1. Compare the number of PDAS iterations with the number of projected-gradient iterations.
2. Change `Active-set parameter` among `0.1`, `1.0`, `10.0` and observe the effect on the active-set evolution.
3. Switch `Continuous control = false` to `true` and compare the resulting control profile.
4. Increase `Control degree` and check how the localization of the active bounds changes.
5. Introduce a nonzero `Forcing term` and verify whether the active region moves.

---

## Reduced vs Monolithic View

This laboratory is useful precisely because the constrained problem is the
same, but the numerical interpretation is different.

### Reduced formulation

- the constraint enters through the feasible set for the control;
- the algorithm is built from state solve, adjoint solve, gradient step, projection;
- the active set is implicit in the projection;
- this viewpoint is natural if one wants to think primarily in optimization terms.

### KKT formulation

- the constraint enters through complementarity conditions on the control block;
- the algorithm is built from active-set prediction and constrained linear solves;
- the active set is explicit and updated from primal-dual information;
- this viewpoint is natural if one wants to think in terms of multipliers, saddle-point systems, and PDAS methods.

The reduced method is often easier to explain first.
The monolithic method is often easier to extend once we move toward more
general complementarity systems.

---

## Laboratory Outline

A useful in-class workflow is:

1. Compile the testbench and run the reduced formulation.
2. Inspect where the control saturates on the bounds.
3. Run the KKT formulation with the same parameters.
4. Compare the final state, adjoint, and control.
5. Compare the iteration histories.
6. Discuss which formulation is easier to extend to:
   variable bounds, discontinuous controls, complementarity systems, richer KKT structures.

---

## Useful Commands

Assuming you are working inside [`codes/dealii`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii):

```bash
cmake -S . -B build
cmake --build build --target laplacian_box_constraints_2d.g kkt_box_constraints_2d.g
./build/laplacian_box_constraints_2d.g parameters/laplacian_box_constraints.prm
./build/kkt_box_constraints_2d.g parameters/kkt_box_constraints.prm
```

If you use ParaView, it is convenient to open:

- `laplacian_box_constraints.pvd` to inspect the projected-gradient evolution;
- `laplacian_box_constraints.vtu` and `kkt_box_constraints.vtu` for the final comparison.

---

## Guiding Questions

- Why is the reduced optimality condition a variational inequality rather than $\nabla \hat J(u)=0$?
- In what sense does projected gradient follow the tangent cone to the admissible set?
- Why does the PDAS method turn box constraints into equality constraints on subsets of control DoFs?
- Why is `AffineConstraints<double>` the natural deal.II mechanism for those active control constraints?
- Which formulation would you choose if the next step were semismooth Newton or more general complementarity constraints?

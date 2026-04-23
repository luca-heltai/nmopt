# deal.II Box-Constrained Optimal Control

## Laboratory Goal

In this laboratory we compare two strategies for an elliptic optimal control problem with **box constraints**

$$
u_a(x) \le u(x) \le u_b(x),
$$

and cost functional

$$
J(y,u) = \frac12 \|y-y_d\|_{L^2(\Omega)}^2 + \frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2,
$$

subject to

$$
-\Delta y = u + f
\qquad \text{in }\Omega,
\qquad
y = 0
\qquad \text{on }\partial\Omega.
$$

The two formulations we want to experiment with are:

- **reduced formulation**: eliminate the state and optimize only with respect to the control by means of a **projected gradient method**;
- **all-at-once formulation**: solve the KKT system with box constraints through a **primal-dual active set** strategy.

The main didactic goal is not only to compute a solution, but also to see how the following change:

- the code structure;
- the meaning of discrete optimality;
- the numerical behavior as $\alpha$, the mesh, and the bounds vary.

---

## Relevant Files

Two new executables have been added to the deal.II testbench:

- [`codes/dealii/execs/laplacian_box_constraints.cc`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/execs/laplacian_box_constraints.cc)
  reduced formulation with projected gradient;
- [`codes/dealii/execs/kkt_box_constraints.cc`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/execs/kkt_box_constraints.cc)
  all-at-once formulation with primal-dual active set.

The associated parameter files are:

- [`codes/dealii/parameters/laplacian_box_constraints.prm`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/parameters/laplacian_box_constraints.prm)
- [`codes/dealii/parameters/kkt_box_constraints.prm`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/parameters/kkt_box_constraints.prm)

In addition, [`codes/dealii/include/kkt.h`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/include/kkt.h)
now exposes a minimal hook to save an externally constructed solution, which is useful for the active-set loop.

---

## Part 1: Reduced Formulation

The idea is:

1. fix a control $u$;
2. solve the state equation for $y(u)$;
3. solve the adjoint equation for $p(u)$;
4. build the reduced gradient;
5. take a descent step and then **project** the control back onto the admissible interval.

In the code the projection is performed coefficient-wise:

$$
u^{k+1} = P_{[u_a,u_b]}(u^k - \alpha_k \nabla \hat J(u^k)).
$$

### What to observe

- the optimality residual is not the plain gradient, but the **projected-gradient residual**;
- once the control hits the bounds, unconstrained descent is no longer the correct dynamics;
- the `.vtu` output also stores the fields `lower_active`, `upper_active`, `inactive`.

### Suggested experiments

1. Reduce `Regularization` from `1e-4` to `1e-6` and observe how the control tends to saturate.
2. Enlarge the bounds from `[-0.2,0.2]` to `[-1,1]` and compare with the almost-unconstrained case.
3. Change `Target state` from a step to `sin(pi*x)` and check whether the active set becomes smaller.
4. Increase `Global refinements` and see whether the active regions stabilize.
5. Modify the `Initial control` and compare the number of projected-gradient iterations.

---

## Part 2: All-at-Once KKT

Here we do not eliminate the control and we do not work only with the reduced functional.
Instead, we solve the discrete optimality system and impose the bounds through a **primal-dual active set** strategy:

- estimate the lower active set $A_-$;
- estimate the upper active set $A_+$;
- impose $u=u_a$ on $A_-$ and $u=u_b$ on $A_+$;
- solve the KKT system on the remaining degrees of freedom;
- update the active sets until convergence.

This procedure is often interpreted as a semismooth Newton method for the complementarity conditions.

### What to observe

- the method often converges in a small number of active-set iterations;
- active degrees of freedom become Dirichlet-type conditions on the control block;
- the monolithic formulation is more natural when adding multipliers, penalties, or richer constraints.

### Suggested experiments

1. Compare the number of PDAS iterations with the number of projected-gradient iterations.
2. Change `Active-set parameter` among `0.1`, `1.0`, `10.0` and observe the effect on the active-set evolution.
3. Switch `Continuous control = false` to `true` and compare the resulting control profile.
4. Increase `Control degree` and check how the localization of active bounds changes.
5. Introduce a nonzero `Forcing term` and verify whether the active region moves.

---

## Laboratory Outline

A possible in-class workflow is:

1. Compile the testbench and run the reduced formulation.
2. Open the `.vtu` output and identify where the control saturates on the bounds.
3. Run the KKT formulation with the same parameters.
4. Compare control, state, and final cost.
5. Change `Regularization` and repeat.
6. Discuss which formulation is easier to extend to:
   pointwise constraints, spatially varying bounds, discontinuous controls, iterative solvers.

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

- In what sense does projected gradient "see" the tangent cone?
- Why is the optimality condition no longer simply $\nabla \hat J(u)=0$?
- When is it preferable to eliminate the state, and when is it preferable to keep the full KKT system?
- What happens if the bounds become incompatible with the target?
- How would you extend one of the two codes to time-dependent bounds or parametric bounds?

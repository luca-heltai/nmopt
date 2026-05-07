# deal.II Laboratory: Forward-Backward Projected Gradient for Parabolic Control

## Laboratory Goal

This laboratory turns the fully discrete optimality system of the previous
lecture into a reduced algorithm.

We do not rederive the parabolic theory, the active-set machinery, or the
discrete adjoint.  We use them.

The computational pattern is:

1. choose a space-time control sequence;
2. solve the heat equation forward in time;
3. solve the discrete adjoint equation backward in time;
4. compute the reduced gradient at every time level;
5. project the gradient step onto the box constraints.

The algorithm is implemented in the new deal.II class
[`ParabolicProjectedGradient`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/include/parabolic_projected_gradient.h)
and in the source file
[`codes/dealii/source/parabolic_projected_gradient.cc`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/source/parabolic_projected_gradient.cc).

The executable driver is
[`codes/dealii/execs/parabolic_projected_gradient.cc`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/execs/parabolic_projected_gradient.cc),
with default parameters in
[`codes/dealii/parameters/parabolic_projected_gradient.prm`](/Users/heltai/latex/courses/2025/02_nmopt/codes/dealii/parameters/parabolic_projected_gradient.prm).

The default experiment is one-dimensional and intentionally small, so that
the whole optimization history can be written and inspected quickly.

---

## Problem Solved in the Code

The state equation is the linear parabolic problem

$$
y_t-\nabla\cdot(\mu\nabla y)+\sigma y = u+f,
\qquad y(0)=y_0,
$$

with homogeneous Dirichlet boundary conditions in the default parameter file.

The cost functional is

$$
J(y,u)
=
\frac12\int_0^T \|y(t)-y_d(t)\|_{L^2(\Omega)}^2\,dt
+\frac{\alpha}{2}\int_0^T \|u(t)\|_{L^2(\Omega)}^2\,dt
+\frac{\beta}{2}\|y(T)-y_T\|_{L^2(\Omega)}^2.
$$

The control is constrained by pointwise bounds

$$
u_a(x,t)\le u(x,t)\le u_b(x,t).
$$

All data are parsed functions, so the same code can be used for stationary
targets, moving targets, disappearing targets, different bounds, forcing
terms, and terminal targets.

The default desired state is a spatial step that appears after an initial
waiting time and disappears before the final time:

```text
x > 0.35 ? (x < 0.65 ? (t > 0.25 ? (t < 0.75 ? 1 : 0) : 0) : 0) : 0
```

The default terminal target is zero.

---

## Discrete Loop

The implementation uses the implicit Euler equations from the previous
lecture.

For each optimization iterate, the state solve is the forward loop

$$
(M+\tau K)Y^k
=
MY^{k-1}+\tau B U^k+\tau F^k,
\qquad k=1,\ldots,N.
$$

The adjoint solve is the backward loop

$$
(M+\tau K)^T P^k
=
M^T P^{k+1}
+\tau M(Y^k-Y_d^k),
\qquad k=N,\ldots,1,
$$

with the terminal tracking term added at the final time level.

The reduced gradient is computed as the Riesz representative in the control
space:

$$
G^k
=
\alpha U^k + M_u^{-1}B^T P^k.
$$

The projected update is level-wise:

$$
U^{k,+}
=
P_{[U_a^k,U_b^k]}
\left(U^k-sG^k\right).
$$

The code measures the projected-gradient residual

$$
\left\|
\frac{U-P_{[U_a,U_b]}(U-sG)}{s}
\right\|_{L^2(0,T;L^2(\Omega))}.
$$

This is the stopping quantity for the constrained problem.

---

## Build and Run

From the deal.II laboratory directory:

```bash
cd codes/dealii
cmake -S . -B build
cmake --build build --target parabolic_projected_gradient_1d.g
./build/parabolic_projected_gradient_1d.g parameters/parabolic_projected_gradient.prm
```

The code also configures dimension-dependent targets in the same style as the
other laboratories:

```bash
cmake --build build --target parabolic_projected_gradient_2d.g
cmake --build build --target parabolic_projected_gradient_3d.g
```

The supplied parameter file is written for the default 1D laboratory.  In 2D,
change the parsed-function variable names from

```text
x,t
```

to

```text
x,y,t
```

and adapt the target expression if it should depend on the second coordinate.

---

## Output

Every optimization iteration writes a complete time sequence:

```text
parabolic_projected_gradient-opt0000.pvd
parabolic_projected_gradient-opt0001.pvd
...
```

Each `.pvd` file points to all time levels for that optimization iterate.
Open one of these files in ParaView and inspect:

- `state`;
- `adjoint`;
- `control`;
- `desired_state`;
- `terminal_target`;
- `lower_bound`;
- `upper_bound`.

The file

```text
parabolic_projected_gradient-history.csv
```

contains one line per optimization iteration:

- total cost;
- tracking cost;
- terminal tracking cost;
- control cost;
- gradient norm;
- projected-gradient norm;
- control norm;
- step length.

This CSV file is the simplest way to see whether the method is actually
decreasing the objective.

---

## Code Map

The main methods to inspect are:

- `solve_state()`: forward implicit Euler loop;
- `solve_adjoint()`: backward discrete adjoint loop;
- `compute_gradient()`: mass-scaled control gradient;
- `projected_gradient_step()`: box projection at each time level;
- `compute_projected_gradient_norm()`: constrained residual;
- `choose_step_length()`: currently returns the fixed step length.

The last method is deliberately small.  It is the place where the Armijo
exercise belongs.

---

## Exercise 1: Baseline Run

Run the default parameter file.

Tasks:

1. Open `parabolic_projected_gradient-opt0000.pvd`.
2. Open the last generated `.pvd` file.
3. Compare `state`, `control`, and `desired_state`.
4. Read the CSV file and check whether the cost decreases.

Questions:

- At which times does the control become most active?
- Does the state continue moving after the desired step disappears?
- How does the terminal target influence the final part of the trajectory?

---

## Exercise 2: Regularization

Repeat the run with

$$
\alpha=10^{-1},\quad 10^{-3},\quad 10^{-5}.
$$

Tasks:

1. Compare the control norm in the CSV files.
2. Compare tracking quality in ParaView.
3. Identify when the bounds become active.

Expected conclusion:

- large regularization gives weaker controls and poorer tracking;
- small regularization gives more aggressive controls and more saturation;
- the time-dependent problem makes this trade-off visible not only in space,
  but also in time.

---

## Exercise 3: Time Resolution

Change only the number of time steps:

```text
Number of time steps = 10
Number of time steps = 20
Number of time steps = 40
```

Tasks:

1. Compare the state sequence.
2. Compare the adjoint sequence.
3. Compare the optimization history.

Question:

What changes when the desired state switches on and off between two coarse
time levels?

---

## Exercise 4: Moving Targets

Replace the default desired state by a moving interval or a moving bump.

For example, in 1D:

```text
x > 0.2 + 0.4*t ? (x < 0.4 + 0.4*t ? 1 : 0) : 0
```

or a smooth moving target:

```text
exp(-100*(x-(0.25+0.5*t))^2)
```

Tasks:

1. Visualize the target and state side by side.
2. Identify whether the control anticipates the target motion.
3. Compare discontinuous and smooth targets.

Question:

Why is the adjoint useful for understanding where the control should act
before the mismatch is visible in the state?

---

## Exercise 5: Box Constraints

Change the bounds:

```text
[-2, 2]
[-0.5, 0.5]
[-0.1, 0.1]
```

Tasks:

1. Look at `control`, `lower_bound`, and `upper_bound`.
2. Compare the projected-gradient norm in the CSV file.
3. Decide whether the fixed step length is still reasonable.

Question:

Why is the projected-gradient residual more meaningful than the unconstrained
gradient norm once the control saturates?

---

## Central Exercise: Projected Armijo

The current implementation uses a fixed step length:

```cpp
return fixed_step_length;
```

inside `choose_step_length()`.

Replace this by projected Armijo backtracking.

For a trial step $s$, define

$$
U_s=P_{[U_a,U_b]}(U-sG),
\qquad
D_s=U_s-U.
$$

The Armijo condition should be written using the reduced cost:

$$
j(U_s)
\le
j(U)+c\, (G,D_s)_{L^2(0,T;L^2(\Omega))}.
$$

Notice that $D_s$ already contains the step length.  Do not multiply the
right-hand side by $s$ again.

A practical implementation plan is:

1. start with `armijo_parameters.alpha0`;
2. build `trial_control = control`;
3. call `projected_gradient_step(trial_control, gradient, step)`;
4. solve the state equation for `trial_control`;
5. compute the trial diagnostics;
6. build the direction `trial_control - control`;
7. compute `(gradient, direction)` with `control_inner_product()`;
8. accept the step if the Armijo condition holds;
9. otherwise multiply the step by `armijo_parameters.beta`;
10. stop if the step falls below `armijo_parameters.minimum_alpha`.

The relevant parameters are already present in the file:

```text
subsection Armijo exercise parameters
  set Initial step length    = 1.0
  set Armijo coefficient     = 1e-4
  set Backtracking reduction = 0.5
  set Max backtracks         = 20
  set Minimum step length    = 1e-12
end
```

After implementing Armijo, repeat Exercises 2 and 5.

Questions:

- Does Armijo choose smaller steps when the bounds are tight?
- Does the cost decrease monotonically?
- Does the number of optimization iterations improve?

---

## Final Discussion

This laboratory is the reduced counterpart of the all-at-once structure from
the previous lecture.

The all-at-once matrix is useful for understanding the global saddle-point
problem.  The reduced loop is useful for implementing scalable algorithms:

- one forward parabolic solve;
- one backward parabolic solve;
- one gradient update in the control space.

The important implementation lesson is:

> once the time-discrete adjoint is correct, a full reduced gradient evaluation
> costs one forward solve and one backward solve, independently of the number
> of control degrees of freedom.

That is the computational reason for introducing adjoints in the first place.

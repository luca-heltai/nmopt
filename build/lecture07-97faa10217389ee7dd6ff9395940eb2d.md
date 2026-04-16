# KKT System for PDE-Constrained Optimization

## Recap

Consider the linear-quadratic optimal control problem

$$
\min_{y,u} \frac12 \|y-y_d\|_{L^2(\Omega)}^2 + \frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2
$$

subject to the state equation

$$
-\nabla\cdot(\mu \nabla y) + \beta\cdot\nabla y + \sigma y = u + f
\qquad\text{in }\Omega,
$$

with suitable boundary conditions.

The three unknowns are:

- $y$: the state;
- $p$: the adjoint;
- $u$: the control.

The code assembles and solves the all-at-once KKT system for $(y,p,u)$.
In block form, the discrete system is

$$
\begin{pmatrix}
M & A^T & 0 \\
A & 0 & -B \\
0 & -B^T & \alpha M_u
\end{pmatrix}
\begin{pmatrix}
y\\ p\\ u
\end{pmatrix}
=
\begin{pmatrix}
M y_d\\ F\\ 0
\end{pmatrix}.
$$

Here:

- $A$ is the discretization of diffusion-reaction-transport;
- $M$ is the state mass matrix;
- $M_u$ is the control mass matrix;
- $\alpha$ is the regularization parameter.

The goal of the laboratory is to understand, by computation, how the solution changes when we modify:

- the regularization parameter $\alpha$;
- the PDE coefficients $\mu,\sigma,\beta$;
- the regularity of the desired state $y_d$.

---

## Compile and Run

Assume `deal.II` is already installed and visible to CMake.

From the directory

```bash
cd codes/dealii
```

configure and compile with

```bash
cmake -S . -B build
cmake --build build --target kkt_2d.g
```

This creates the executable

```bash
./build/kkt_2d.g
```

You can also build the non-debug target with

```bash
cmake --build build --target kkt_2d
```

### Parameter File

The executable reads a parameter file of the form

```bash
./build/kkt_2d.g kkt.prm
```

If `kkt.prm` does not exist, the executable stops and writes a default one for you.
So the standard workflow is:

```bash
./build/kkt_2d.g
./build/kkt_2d.g kkt.prm
```

Then modify `kkt.prm` and run again.

The most important entries are:

- `Regularization`
- `State degree`
- `Control degree`
- `Continuous control`
- `Grid generator arguments`
- `Diffusion coefficient / Function expression`
- `Reaction coefficient / Function expression`
- `Transport field / Function expression`
- `Desired state / Function expression`

After each run, the code writes a `.vtu` file that can be visualized with ParaView.

---

## What To Look At

For every experiment, inspect:

- `state`
- `adjoint`
- `control`
- `desired_state`

and answer the same three questions:

1. How well does the state track the desired state?
2. How large or oscillatory is the control?
3. Which feature of the PDE explains what you see?

---

## Exercise 1: Baseline Run

Run the default case.

Tasks:

1. Plot `state`, `control`, and `desired_state`.
2. Describe where the state fails to match the desired state.
3. Explain why the control is not identically zero.

Question:

Why does a discontinuous target already suggest that exact tracking is difficult?

---

## Exercise 2: Role of Regularization

Keep all PDE coefficients fixed and vary only the regularization parameter:

$$
\alpha = 10^{-1},\;10^{-3},\;10^{-6}.
$$

Tasks:

1. Compare the size of the control for the three runs.
2. Compare the quality of the tracking.
3. Decide which solution is the smoothest and which is the most aggressive.

Question:

What is the trade-off introduced by the term

$$
\frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2?
$$

Expected conclusion:

- large $\alpha$ penalizes the control more;
- small $\alpha$ improves tracking but produces a larger and less regular control.

---

## Exercise 3: Diffusion

Keep $\sigma=0$ and $\beta=0$, and vary only the diffusion coefficient:

$$
\mu = 10^{-2},\;10^{-1},\;1,\;10.
$$

Tasks:

1. Compare the sharpness of the state profile.
2. Compare how difficult it is to reproduce the jump in `desired_state`.
3. Compare the amplitude of the control.

Question:

Why does increasing diffusion make the state smoother?

Expected conclusion:

Stronger diffusion damps sharp transitions, so the state cannot easily reproduce a target with jumps.

---

## Exercise 4: Reaction

Fix diffusion and transport, and vary the reaction coefficient:

$$
\sigma = 0,\;1,\;10.
$$

Tasks:

1. Compare the magnitude of the state.
2. Compare how much control is needed.
3. Compare the shape of the adjoint.

Question:

How does the reaction term change the local balance in the PDE?

---

## Exercise 5: Transport

Fix diffusion and reaction, and vary the transport field.

Suggested choices:

- $\beta=(0,0)$
- $\beta=(1,0)$
- $\beta=(5,0)$
- $\beta=(y,-x)$

Tasks:

1. Identify the preferred direction induced by the transport field.
2. Compare the displacement of the state relative to the target.
3. Describe how the control adapts to the flow.

Question:

What is the qualitative difference between diffusion-dominated behavior and transport-dominated behavior?

---

## Exercise 6: Continuous vs Discontinuous Control

Repeat one of the previous experiments with:

- `Continuous control = false`, `Control degree = 0`
- `Continuous control = true`, `Control degree = 1`

Tasks:

1. Compare the visual regularity of the control.
2. Compare the tracking quality.
3. Compare the number of control degrees of freedom.

Question:

What do you gain and what do you lose when you force the control to be continuous?

---

## Exercise 7: Smooth Target vs Discontinuous Target

Compare two desired states:

1. a smooth target, for example
   $$
   y_d(x,y)=e^{-10(x^2+y^2)};
   $$
2. the discontinuous default target
   $$
   y_d(x,y)=
   \begin{cases}
   1 & \text{if } x^2+y^2<0.5^2,\\
   0 & \text{otherwise.}
   \end{cases}
   $$

Tasks:

1. Compare the quality of the tracking.
2. Compare the regularity of the control.
3. Compare where the largest mismatch occurs.

Question:

Which target is more compatible with the regularity naturally imposed by the PDE?

---

## Exercise 8: Desired State Not in $H^1$

Focus on the discontinuous target.

Tasks:

1. Explain why the target has a jump.
2. Explain why such a function is not in $H^1(\Omega)$.
3. Compare this fact with the regularity expected for the state variable.
4. Identify where the state smooths the jump.

Question:

Why should we expect the optimizer to approximate the target rather than reproduce it exactly?

Expected conclusion:

The state is constrained by an elliptic or transport-diffusion PDE and belongs to a smoother space than the discontinuous target. The optimal solution is therefore a compromise between:

- fitting the target;
- satisfying the PDE;
- controlling the cost of the control.

---

## What To Submit

For each exercise, include:

1. the parameter values used;
2. one plot of `state`, `control`, and `desired_state`;
3. a short comment answering the question of the exercise.

In the final summary, explain in a few lines:

- the role of regularization;
- the effect of diffusion, reaction, and transport;
- what changes when the desired state is not in $H^1(\Omega)$.

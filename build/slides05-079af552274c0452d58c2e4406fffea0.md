# Optimality Conditions with Constraints

Numerical Methods for Optimal Control (NMOPT)

Luca Heltai (<luca.heltai@unipi.it>)

----

## Context

The previous lecture already gave:

- the reduced problem
- the operatorial all-at-once formulation
- the reduced gradient
- the unconstrained KKT system

Today:

- add control constraints
- write constrained KKT conditions
- derive projection and active-set structure
- connect to algorithms

----

## Operatorial Setting

Use

$$
V:=H_0^1(\Omega),
\qquad
Q:=L^2(\Omega).
$$

Operators:

$$
A:V\to V',
\qquad
B:Q\to V',
\qquad
M:V\to V',
\qquad
R_Q:Q\to Q'.
$$

State equation:

$$
Ay-Bu=F
\qquad \text{in }V'.
$$

----

## Cost and Unconstrained KKT

Cost functional:

$$
J(y,u)
=
\frac12\langle M(y-y_d),y-y_d\rangle
+\frac{\alpha}{2}\langle R_Q u,u\rangle.
$$

From the previous lecture:

$$
\begin{aligned}
A y - B u &= F,\\
M y - A^* p &= M y_d,\\
\alpha R_Q u + B^* p &= 0.
\end{aligned}
$$

Reduced gradient:

$$
\nabla f(u)=p+\alpha u.
$$

----

## Admissible Controls

Box constraints:

$$
U_{\mathrm{ad}}
=
\{u\in Q: u_{\min}\le u\le u_{\max}\ \text{a.e.}\}.
$$

The reduced problem becomes

$$
\min_{u\in U_{\mathrm{ad}}} f(u).
$$

Only the control equation changes.

----

## Variational Inequality

Optimality is now

$$
f'(\bar u)(u-\bar u)\ge 0
\qquad \forall u\in U_{\mathrm{ad}}.
$$

Using the reduced gradient:

$$
(\bar p+\alpha \bar u,u-\bar u)_Q\ge 0
\qquad \forall u\in U_{\mathrm{ad}}.
$$

This replaces $\nabla f(\bar u)=0$.

----

## Operatorial KKT

Using the normal cone $N_{U_{\mathrm{ad}}}(\bar u)$:

$$
0\in \alpha R_Q \bar u + B^*\bar p + N_{U_{\mathrm{ad}}}(\bar u)
\qquad \text{in }Q'.
$$

Hence

$$
\begin{aligned}
A \bar y - B \bar u &= F &&\text{in }V',\\
M \bar y - A^* \bar p &= M y_d &&\text{in }V',\\
0 &\in \alpha R_Q \bar u + B^* \bar p + N_{U_{\mathrm{ad}}}(\bar u)
&&\text{in }Q'.
\end{aligned}
$$

----

## Saddle-Point Structure

KKT system in dual spaces:

$$
\begin{pmatrix}
M & 0 & -A^*\\
0 & \alpha R_Q & B^*\\
A & -B & 0
\end{pmatrix}
\begin{pmatrix}
\bar y\\
\bar u\\
\bar p
\end{pmatrix}
+
\begin{pmatrix}
0\\
\eta\\
0
\end{pmatrix}
=
\begin{pmatrix}
M y_d\\
0\\
F
\end{pmatrix},
\qquad
\eta\in N_{U_{\mathrm{ad}}}(\bar u).
$$

Linear saddle point + nonlinear control block.

----

## Projection Formula

For box constraints:

$$
\bar u
=
\Pi_{[u_{\min},u_{\max}]}
\left(-\frac1\alpha \bar p\right).
$$

Pointwise:

$$
\bar u(x)=
\begin{cases}
u_{\min}(x), & -\dfrac1\alpha \bar p(x)\le u_{\min}(x),\\[4pt]
u_{\max}(x), & -\dfrac1\alpha \bar p(x)\ge u_{\max}(x),\\[4pt]
-\dfrac1\alpha \bar p(x), & \text{otherwise}.
\end{cases}
$$

----

## Strong-Form Optimality System

For the Poisson model:

$$
\begin{cases}
-\Delta \bar y=\bar u+f,\\[4pt]
-\Delta \bar p=\bar y-y_d,\\[4pt]
\bar u=\Pi_{U_{\mathrm{ad}}}\!\left(-\dfrac1\alpha \bar p\right),\\[4pt]
\bar y=\bar p=0 \text{ on }\partial\Omega.
\end{cases}
$$

----

## Projected Gradient

At iteration $k$:

1. solve state

   $$
   A y^k-Bu^k=F
   $$

2. solve adjoint

   $$
   M y^k-A^*p^k=M y_d
   $$

3. project

   $$
   u^{k+1}
   =
   \Pi_{U_{\mathrm{ad}}}\!\left(u^k-\tau_k(p^k+\alpha u^k)\right)
   $$

For $\tau_k=1/\alpha$:

$$
u^{k+1}=\Pi_{U_{\mathrm{ad}}}\!\left(-\frac1\alpha p^k\right).
$$

----

## Primal-Dual Active Set

Define active sets from $-\alpha^{-1}p^k$:

$$
\mathcal A_-^k=\left\{-\frac1\alpha p^k\le u_{\min}\right\},
\qquad
\mathcal A_+^k=\left\{-\frac1\alpha p^k\ge u_{\max}\right\}.
$$

Inactive set:

$$
\mathcal I^k=\Omega\setminus(\mathcal A_-^k\cup\mathcal A_+^k).
$$

Then solve with sets frozen:

$$
u^{k+1}=u_{\min}\text{ on }\mathcal A_-^k,\qquad
u^{k+1}=u_{\max}\text{ on }\mathcal A_+^k,
$$

$$
\alpha u^{k+1}+p^{k+1}=0
\qquad \text{on }\mathcal I^k.
$$

----

## PDAS Multipliers

Introduce

$$
\lambda\in Q'
$$

through

$$
\alpha R_Q u + B^* p + \lambda = 0.
$$

On frozen active sets, use a Lagrangian with equality constraints:

$$
\mathcal L_k(y,u,p,\lambda_-^k,\lambda_+^k)
=
J(y,u)-\langle Ay-Bu-F,p\rangle
+(\lambda_-^k,u-u_{\min})_{\mathcal A_-^k}
+(\lambda_+^k,u-u_{\max})_{\mathcal A_+^k}.
$$

PDAS = solve this equality-constrained KKT system, update the sets, repeat.

----

## Beyond Box Constraints

If $U_{\mathrm{ad}}\subset Q$ is closed and convex and the projection exists, one still has

$$
u=\Pi_{U_{\mathrm{ad}}}\!\left(u-\sigma \nabla f(u)\right),
$$

equivalently

$$
0\in \nabla f(u)+N_{U_{\mathrm{ad}}}(u).
$$

Active-set ideas still apply, but the active geometry is less explicit than in the box case.

----

## Semismooth Newton

Generalized derivative = Clarke generalized Jacobian:

$$
\partial G(u)
=
\operatorname{co}
\left\{
\lim_{j\to\infty} G'(u_j):
u_j\to u,\;
G \text{ differentiable at }u_j
\right\}.
$$

Semismooth definition:

$$
G(u+s)-G(u)-V_s s=o(\|s\|)
\qquad
\text{for }V_s\in \partial G(u+s).
$$

This replaces classical differentiability in generalized Newton methods.

Define

$$
G(u):=
u-\Pi_{U_{\mathrm{ad}}}\!\left(u-\sigma \nabla f(u)\right).
$$

Solve

$$
G(u)=0
$$

by generalized Newton.

For box constraints:

- generalized derivative of projection is `0` on active points
- generalized derivative is `1` on inactive points

This is why PDAS and semismooth Newton essentially coincide here.

----

## Comparison

- `Projected gradient`: robust, simple, slower
- `PDAS`: active-set identification + linear KKT solve
- `Semismooth Newton`: Newton-type local behavior for the nonsmooth equation

All three come directly from the same continuous optimality system.

----

## Takeaways

- control constraints change only the control optimality condition
- the constrained problem is still a saddle-point problem in

  $$
  V'\times Q'\times V'
  $$

- box constraints give the projection formula

  $$
  u=\Pi_{U_{\mathrm{ad}}}\!\left(-\frac1\alpha p\right)
  $$

- this structure naturally leads to projected gradient, PDAS, and semismooth Newton

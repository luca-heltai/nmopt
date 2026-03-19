# Optimality Conditions for Elliptic Optimal Control with Constraints

## Overview

The previous lecture established the continuous framework for linear-quadratic elliptic
control:

- the reduced problem

  $$
  \min_{u\in Q} f(u),
  $$

- the operatorial all-at-once formulation on

  $$
  V\times Q\times V,
  \qquad
  V:=H_0^1(\Omega),\quad Q:=L^2(\Omega),
  $$

- the adjoint-based gradient formula

  $$
  \nabla f(u)=p+\alpha u,
  $$

- and the unconstrained KKT system in dual spaces.

This lecture extends the framework to the case of
control constraints.

The goal is to derive:

1. the variational inequality for box-constrained controls;
2. the operatorial KKT system with a nonlinear control condition;
3. the pointwise projection formula;
4. the saddle-point interpretation in infinite dimensions.

---

## Notation and Setting

$$
V:=H_0^1(\Omega), \qquad Q:=L^2(\Omega).
$$

The operators are

$$
A:V\to V',
\qquad
B:Q\to V',
\qquad
M:V\to V',
\qquad
R_Q:Q\to Q',
$$

with

$$
\langle Ay,v\rangle_{V',V}=\int_\Omega \nabla y\cdot \nabla v\,dx,
\qquad
\langle Bu,v\rangle_{V',V}=(u,v)_Q,
$$

$$
\langle My,v\rangle_{V',V}=(y,v)_Q,
\qquad
\langle R_Q u,w\rangle_{Q',Q}=(u,w)_Q.
$$

The state equation is

$$
Ay-Bu=F
\qquad \text{in }V',
$$

and the cost functional is

$$
J(y,u)
=
\frac12\langle M(y-y_d),y-y_d\rangle_{V',V}
+\frac{\alpha}{2}\langle R_Q u,u\rangle_{Q',Q}.
$$

Hence the all-at-once problem is

$$
\min_{(y,u)\in V\times Q} J(y,u)
\qquad\text{subject to}\qquad
Ay-Bu=F \text{ in }V'.
$$

From the previous lecture, in the unconstrained case, the first-order system is

$$
\begin{aligned}
A y - B u &= F &&\text{in }V',\\
M y - A^* p &= M y_d &&\text{in }V',\\
\alpha R_Q u + B^* p &= 0 &&\text{in }Q'.
\end{aligned}
$$

Equivalently, in reduced form,

$$
\nabla f(u)=p+\alpha u.
$$

This is the only fact we need to import from the previous lecture.

---

## Admissible Controls

We now impose box constraints on the control:

$$
U_{\mathrm{ad}}
:=
\left\{
u\in Q:
u_{\min}(x)\le u(x)\le u_{\max}(x)
\text{ for a.e. }x\in\Omega
\right\}.
$$

We assume

$$
u_{\min},u_{\max}\in L^\infty(\Omega),
\qquad
u_{\min}(x)\le u_{\max}(x)
\quad \text{a.e. in }\Omega,
$$

so that $U_{\mathrm{ad}}\subset Q$ is nonempty, closed, and convex.

The constrained reduced problem is therefore

$$
\min_{u\in U_{\mathrm{ad}}} f(u).
$$

The state equation is unchanged.
Only the admissible set for the control changes.

---

## Variational Inequality

For minimization of a differentiable functional over a closed convex set in a Hilbert
space, the first-order condition is a variational inequality rather than the equation
$\nabla f(\bar u)=0$.

Thus $\bar u\in U_{\mathrm{ad}}$ is optimal if and only if

$$
f'(\bar u)(u-\bar u)\ge 0
\qquad \forall u\in U_{\mathrm{ad}}.
$$

Using the reduced gradient from the previous lecture, this becomes

$$
(\bar p+\alpha \bar u,u-\bar u)_Q\ge 0
\qquad \forall u\in U_{\mathrm{ad}}.
$$

This is the constrained analogue of stationarity.
If the constraint is inactive, one recovers $\bar p+\alpha \bar u=0$.
If a bound is active, only admissible perturbations are allowed, so equality is replaced
by an inequality.

---

## Operatorial KKT System

The cleanest way to express the constrained problem in the operatorial setting is through
the normal cone to $U_{\mathrm{ad}}$.
For a closed convex set $K\subset Q$, the normal cone at $u\in K$ is

$$
N_K(u)
:=
\left\{
\eta\in Q':
\langle \eta,v-u\rangle_{Q',Q}\le 0
\quad \forall v\in K
\right\}.
$$

Then the variational inequality is equivalent to the inclusion

$$
0\in \alpha R_Q \bar u+B^*\bar p+N_{U_{\mathrm{ad}}}(\bar u)
\qquad \text{in }Q'.
$$

Therefore the full KKT system becomes

$$
\begin{aligned}
A \bar y - B \bar u &= F &&\text{in }V',\\
M \bar y - A^* \bar p &= M y_d &&\text{in }V',\\
0 &\in \alpha R_Q \bar u + B^* \bar p + N_{U_{\mathrm{ad}}}(\bar u)
&&\text{in }Q'.
\end{aligned}
$$

Compared with the unconstrained system from the previous lecture, only the third equation changes:
the linear equation in $Q'$ is replaced by a nonlinear inclusion in $Q'$.

This can also be written as a block saddle-point relation:

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
\eta\in N_{U_{\mathrm{ad}}}(\bar u),
$$

in the dual product space

$$
V'\times Q'\times V'.
$$

So the constrained problem still has the same KKT block structure as in the previous lecture,
but now one block is nonlinear because of the cone term.

---

## Projection Formula

For box constraints, the normal-cone inclusion has a pointwise characterization.
It is equivalent to the projection formula

$$
\bar u
=
\Pi_{[u_{\min},u_{\max}]}
\left(-\frac1\alpha \bar p\right),
$$

where $\Pi_{[u_{\min},u_{\max}]}$ denotes the pointwise projection onto the interval
$[u_{\min}(x),u_{\max}(x)]$.

Explicitly, for almost every $x\in\Omega$,

$$
\bar u(x)=
\begin{cases}
u_{\min}(x)
& \text{if } -\dfrac1\alpha \bar p(x)\le u_{\min}(x),\\[6pt]
u_{\max}(x)
& \text{if } -\dfrac1\alpha \bar p(x)\ge u_{\max}(x),\\[6pt]
-\dfrac1\alpha \bar p(x)
& \text{if } u_{\min}(x)<-\dfrac1\alpha \bar p(x)<u_{\max}(x).
\end{cases}
$$

Thus the control law remains explicit, but it is no longer linear.
The unconstrained formula $u=-\alpha^{-1}p$ is simply clipped onto the admissible interval.

---

## Full Optimality System

Combining state equation, adjoint equation, and control characterization, we obtain

$$
\begin{cases}
A \bar y - B \bar u = F & \text{in }V',\\[4pt]
M \bar y - A^* \bar p = M y_d & \text{in }V',\\[4pt]
\bar u = \Pi_{U_{\mathrm{ad}}}\!\left(-\dfrac1\alpha \bar p\right) & \text{in }Q.
\end{cases}
$$

In the Poisson model of this course, this reads

$$
\begin{cases}
-\Delta \bar y = \bar u + f & \text{in }\Omega,\\[4pt]
-\Delta \bar p = \bar y-y_d & \text{in }\Omega,\\[4pt]
\bar u = \Pi_{U_{\mathrm{ad}}}\!\left(-\dfrac1\alpha \bar p\right) & \text{in }\Omega,\\[4pt]
\bar y=\bar p=0 & \text{on }\partial\Omega.
\end{cases}
$$

The first two equations are linear and exactly the same kind of objects as in the previous lecture.
The third one carries the entire effect of the control constraint.

---

## Interpretation as Infinite-Dimensional KKT

The constrained system is the direct analogue of finite-dimensional KKT conditions:

- the unknowns are functions in $V\times Q\times V$;
- the equations live in the dual product space $V'\times Q'\times V'$;
- the state and adjoint equations are linear operator equations;
- the control condition is a complementarity relation encoded either as a normal-cone inclusion
  or as a projection formula.

So the conceptual picture is:

- previous lecture: linear saddle-point system;
- current lecture: saddle-point system with a nonlinear control block.

This is exactly the structure later exploited by active-set and semismooth Newton methods.

---

## Projected Gradient Methods

The most direct numerical method works on the reduced functional over the convex set
$U_{\mathrm{ad}}$.
Starting from $u^0\in U_{\mathrm{ad}}$, one computes at each iteration:

1. state solve

   $$
   A y^k - B u^k = F;
   $$

2. adjoint solve

   $$
   M y^k - A^* p^k = M y_d;
   $$

3. reduced gradient

   $$
   g^k=\nabla f(u^k)=p^k+\alpha u^k;
   $$

4. projected update

   $$
   u^{k+1}
   =
   \Pi_{U_{\mathrm{ad}}}\!\left(u^k-\tau_k g^k\right).
   $$

This is the Hilbert-space analogue of projected gradient descent for bound-constrained
finite-dimensional optimization.

Two remarks are important:

- the projection enforces admissibility at every iteration;
- each step still requires one state solve and one adjoint solve, exactly as in the unconstrained reduced method.

For the linear-quadratic case, if one chooses

$$
\tau_k=\frac1\alpha,
$$

then

$$
u^k-\tau_k g^k
=
u^k-\frac1\alpha(\alpha u^k+p^k)
=
-\frac1\alpha p^k,
$$

so the update becomes

$$
u^{k+1}
=
\Pi_{U_{\mathrm{ad}}}\!\left(-\frac1\alpha p^k\right).
$$

This explains why the projection formula is both an optimality condition and a natural
fixed-point iteration.

Projected gradient is robust and easy to implement, but it is usually only linearly
convergent and may become slow when the active set is nearly identified but not yet fixed.

---

## Primal-Dual Active Set Methods

The projection formula suggests splitting the domain into active and inactive regions.
Given an iterate $(u^k,p^k)$, define

$$
\mathcal A_-^k
:=
\left\{
x\in\Omega:
-\frac1\alpha p^k(x)\le u_{\min}(x)
\right\},
$$

$$
\mathcal A_+^k
:=
\left\{
x\in\Omega:
-\frac1\alpha p^k(x)\ge u_{\max}(x)
\right\},
$$

and the inactive set

$$
\mathcal I^k
:=
\Omega\setminus(\mathcal A_-^k\cup \mathcal A_+^k).
$$

To make the primal-dual structure explicit, introduce a multiplier

$$
\lambda\in Q',
$$

for the box constraints, and write the control optimality condition as

$$
\alpha R_Q u + B^* p + \lambda = 0
\qquad \text{in }Q'.
$$

For box constraints, $\lambda$ satisfies the complementarity pattern

$$
\lambda(x)\ge 0 \ \text{on the lower active set},
\qquad
\lambda(x)\le 0 \ \text{on the upper active set},
\qquad
\lambda(x)=0 \ \text{on the inactive set}.
$$

Once the active sets are frozen, the inequalities are replaced by equality constraints.
On a given active-set guess $(\mathcal A_-^k,\mathcal A_+^k,\mathcal I^k)$, one may therefore
consider the Lagrangian

$$
\begin{aligned}
\mathcal L_k(y,u,p,\lambda_-^k,\lambda_+^k)
:={}&
\frac12\|y-y_d\|_{L^2(\Omega)}^2
+\frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2 \\
&-\langle Ay-Bu-F,p\rangle_{V',V}\\
&+(\lambda_-^k,u-u_{\min})_{L^2(\mathcal A_-^k)}
+(\lambda_+^k,u-u_{\max})_{L^2(\mathcal A_+^k)},
\end{aligned}
$$

where the equality constraints are imposed only on the active sets.
The multipliers $\lambda_-^k$ and $\lambda_+^k$ are supported on
$\mathcal A_-^k$ and $\mathcal A_+^k$, respectively.

Then the optimality conditions suggest the update rules

$$
u^{k+1}=u_{\min}\quad \text{on }\mathcal A_-^k,
\qquad
u^{k+1}=u_{\max}\quad \text{on }\mathcal A_+^k,
$$

and

$$
\alpha u^{k+1}+p^{k+1}=0
\qquad \text{on }\mathcal I^k.
$$

Equivalently, the multiplier can be recovered from

$$
\lambda^{k+1}
=
-\alpha u^{k+1}-p^{k+1},
$$

with

$$
\lambda^{k+1}=0
\qquad \text{on }\mathcal I^k,
$$

and sign restrictions on $\mathcal A_-^k$ and $\mathcal A_+^k$.

So, once the active sets are guessed, one solves a linear system for
$(y^{k+1},u^{k+1},p^{k+1},\lambda^{k+1})$ with those sets frozen.
In strong form this means:

$$
\begin{cases}
-\Delta y^{k+1}=u^{k+1}+f,\\[4pt]
-\Delta p^{k+1}=y^{k+1}-y_d,\\[4pt]
u^{k+1}=u_{\min}\ \text{on }\mathcal A_-^k,\\[4pt]
u^{k+1}=u_{\max}\ \text{on }\mathcal A_+^k,\\[4pt]
\alpha u^{k+1}+p^{k+1}=0\ \text{on }\mathcal I^k.
\end{cases}
$$

Algorithmically:

1. guess active and inactive sets from the current iterate;
2. solve the equality-constrained KKT system with these sets fixed;
3. update the sets and repeat until they stop changing.

This is called a primal-dual active set method because:

- the primal variable $u$ determines where bounds are active;
- the dual variables $p$ and $\lambda$ decide whether a point should be active or inactive.

In practice, PDAS is often much faster than projected gradient because once the active set
is identified, the remaining step is essentially the solution of the correct linearized KKT system.

The same idea extends beyond box constraints.
Let $U_{\mathrm{ad}}\subset Q$ be any closed convex set for which the projection

$$
\Pi_{U_{\mathrm{ad}}}:Q\to U_{\mathrm{ad}}
$$

is available in explicit or computational form.
Then one may still write the optimality condition as

$$
u=\Pi_{U_{\mathrm{ad}}}\!\left(u-\sigma \nabla f(u)\right),
$$

or equivalently as a normal-cone inclusion

$$
0\in \nabla f(u)+N_{U_{\mathrm{ad}}}(u),
$$

and build active-set type methods by locally identifying the part of the constraint
that behaves as an equality constraint at the current iterate.
For box constraints this identification is pointwise and especially transparent,
which is why PDAS takes its simplest form there.
For more general admissible sets, the same principle survives, but the geometry of
the active manifold and the structure of the projected Newton step become more problem-dependent.

---

## Semismooth Newton Interpretation

The projection operator is not differentiable in the classical sense, but it is
semismooth.
This allows one to apply generalized Newton methods to the nonsmooth optimality system.

Recall the basic definition.
Let $G:Q\to Q$ be locally Lipschitz and directionally differentiable.
Its generalized derivative at $u$ is the Clarke generalized Jacobian

$$
\partial G(u)
:=
\operatorname{co}
\left\{
\lim_{j\to\infty} G'(u_j):
u_j\to u,\;
G \text{ differentiable at }u_j
\right\},
$$

that is, the convex hull of all limits of classical derivatives computed at nearby
points where $G$ is differentiable.

Then $G$ is called semismooth at $u$ if, for every perturbation $s\to 0$ and every
generalized derivative $V_s\in \partial G(u+s)$ in the sense of Clarke,

$$
G(u+s)-G(u)-V_s s=o(\|s\|).
$$

If the stronger estimate

$$
G(u+s)-G(u)-V_s s=O(\|s\|^2)
$$

holds, one speaks of strongly semismoothness.
Semismoothness is the regularity that replaces classical differentiability in
generalized Newton theory.

A convenient reduced equation is

$$
G(u)
:=
u-\Pi_{U_{\mathrm{ad}}}\!\left(u-\sigma \nabla f(u)\right)=0,
$$

with some parameter $\sigma>0$.
At a solution, this is equivalent to the variational inequality.

Indeed, in a Hilbert space the projection onto a closed convex set is characterized by

$$
z=\Pi_{U_{\mathrm{ad}}}(w)
\qquad\Longleftrightarrow\qquad
(w-z,v-z)_Q\le 0
\quad \forall v\in U_{\mathrm{ad}}.
$$

Apply this with

$$
w=u-\sigma \nabla f(u),
\qquad
z=u.
$$

Then

$$
G(u)=0
\qquad\Longleftrightarrow\qquad
u=\Pi_{U_{\mathrm{ad}}}\!\left(u-\sigma \nabla f(u)\right).
$$

Using the projection characterization, this is equivalent to

$$
(u-\sigma \nabla f(u)-u,v-u)_Q\le 0
\qquad \forall v\in U_{\mathrm{ad}},
$$

that is,

$$
-\sigma(\nabla f(u),v-u)_Q\le 0
\qquad \forall v\in U_{\mathrm{ad}}.
$$

Since $\sigma>0$, this is equivalent to

$$
(\nabla f(u),v-u)_Q\ge 0
\qquad \forall v\in U_{\mathrm{ad}},
$$

which is exactly the variational inequality

$$
f'(u)(v-u)\ge 0
\qquad \forall v\in U_{\mathrm{ad}}.
$$

So the fixed-point equation $G(u)=0$ and the constrained first-order condition are equivalent.

For box constraints, the projection acts pointwise, so its generalized derivative is
easy to characterize:

- it is $0$ on points that are projected onto an active bound;
- it is $1$ on inactive points.

Therefore a generalized Newton step for $G(u)=0$ amounts to freezing the current active set
and solving the corresponding linearized state-adjoint-control system.
This is exactly the mechanism behind PDAS.

So, in the linear-quadratic box-constrained case, one can view:

- PDAS as an active-set method;
- semismooth Newton as a generalized Newton method for the projection equation;
- and, in this specific setting, the two viewpoints are essentially equivalent.

This equivalence is one of the main reasons why active-set methods are so effective for
elliptic control problems with box constraints: they inherit Newton-type local behavior
while preserving the simple geometric interpretation of active and inactive regions.

---

## Algorithmic Comparison

The three methods fit the same state-adjoint-control structure, but they use it differently:

- projected gradient uses the projection only as a feasible descent step;
- PDAS uses the projection to identify active sets and then solves a linear KKT system;
- semismooth Newton interprets the projection equation as a nonsmooth root-finding problem.

Typical qualitative behavior:

- projected gradient: cheapest iteration, strongest robustness, slowest asymptotically;
- PDAS: more expensive iteration, often very fast once the active set is close to correct;
- semismooth Newton: Newton-type local convergence, but requires a good formulation and linear solver strategy.

This is why the continuous optimality system matters algorithmically:
it tells us not only what the solution must satisfy, but also what class of numerical
methods is natural for the problem.

---

## Summary

This lecture builds directly on the previous one:

1. the state equation and adjoint equation are unchanged;
2. the unconstrained stationarity condition

   $$
   \alpha R_Q u+B^*p=0
   $$

   is replaced by the constrained inclusion

   $$
   0\in \alpha R_Q u+B^*p+N_{U_{\mathrm{ad}}}(u);
   $$

3. equivalently, the control satisfies the projection formula

   $$
   u=\Pi_{U_{\mathrm{ad}}}\!\left(-\frac1\alpha p\right);
   $$

4. the full KKT system remains a saddle-point problem in dual spaces

   $$
   V'\times Q'\times V',
   $$

   with a nonlinear control block.
5. this structure naturally leads to three algorithmic families:

   - projected gradient methods;
   - primal-dual active set methods;
   - semismooth Newton methods.

This closes the continuous first-order theory for linear-quadratic elliptic control
with box constraints.

## References

- Fredi Troeltzsch, *Optimal Control of Partial Differential Equations*, Chapter 2.
- Alfio Borzi and Volker Schulz, *Computational Optimization of Systems Governed by Partial Differential Equations*.
- Juan Carlos De Los Reyes, *Numerical PDE-Constrained Optimization*.

<!-- FOOTER START -->
<iframe src="https://luca-heltai.github.io/nmopt/slideshow/slides05.html" width="100%" height="800px" style="border: none;"></iframe>

---

```{admonition} 🎬 View Slides
:class: tip

**[Open slides in full screen](https://luca-heltai.github.io/nmopt/slideshow/slides05.html)** for the best viewing experience.
```
<!-- FOOTER END -->

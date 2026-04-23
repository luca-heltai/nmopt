# From ODE-Constrained Optimization to Bochner Spaces

## Overview

The previous lectures developed the continuous and discrete theory for
elliptic optimal control, including:

- reduced and all-at-once formulations;
- adjoint-based gradient representations;
- variational inequalities and normal-cone KKT systems;
- projected-gradient and primal-dual viewpoints.

At this stage, the next genuine change of perspective is the introduction of
**time-dependent dynamics**.
The simplest rigorous entry point is not yet a parabolic PDE, but a linear
control problem governed by an ordinary differential equation.
This lecture has two goals:

1. derive the optimality system for a linear-quadratic control problem governed by an ODE;
2. introduce the vector-valued function spaces needed to pass from ODEs to parabolic PDEs.

The logical path of the lecture is the following:

1. formulate a linear-quadratic control problem on $(0,T)$;
2. prove well-posedness of the state equation and continuity of the control-to-state map;
3. derive the reduced derivative through a linearized state equation;
4. introduce the adjoint ODE and obtain the gradient formula;
5. state the first-order optimality conditions, both unconstrained and box-constrained;
6. isolate the structural difference with the elliptic case: a forward-backward optimality system;
7. introduce strongly measurable vector-valued functions and Bochner spaces;
8. define weak time derivatives in dual spaces and the space
   $$
   W(0,T):=\{y\in L^2(0,T;V): y_t\in L^2(0,T;V')\};
   $$
9. state the basic embedding and integration-by-parts results that make parabolic optimal control possible.

---

## Model Problem

Let $T>0$ be fixed, and let

$$
A\in \mathbb R^{n\times n},
\qquad
B\in \mathbb R^{n\times m},
\qquad
y_0\in \mathbb R^n.
$$

Let also

$$
f\in L^2(0,T;\mathbb R^n),
\qquad
y_d\in L^2(0,T;\mathbb R^n),
\qquad
y_T\in \mathbb R^n,
$$

and let $\alpha>0$, $\beta\ge 0$.
The control space is

$$
U:=L^2(0,T;\mathbb R^m).
$$

The admissible set is a nonempty, closed, convex subset

$$
U_{\mathrm{ad}}\subset U.
$$

For $u\in U$, the state $y=y(u)$ is defined by the linear ODE

$$
\dot y(t)=Ay(t)+Bu(t)+f(t)
\qquad \text{for a.e. }t\in(0,T),
\qquad
y(0)=y_0.
$$

We consider the optimal control problem

$$
\min_{u\in U_{\mathrm{ad}}} j(u):=J(y(u),u),
$$

where

$$
J(y,u)
:=
\frac12\int_0^T |y(t)-y_d(t)|^2\,dt
+\frac\alpha2\int_0^T |u(t)|^2\,dt
+\frac\beta2 |y(T)-y_T|^2.
$$

The three terms are:

- a distributed tracking term in time;
- the Tikhonov regularization on the control;
- an optional terminal observation.

This model already contains all structural ingredients of parabolic optimal
control:

- an evolution equation;
- a cost integrated over time;
- a terminal contribution;
- an adjoint equation with a final condition.

---

## State Equation and Control-to-State Map

We begin with the well-posedness of the state equation.
Since the right-hand side belongs to $L^2(0,T;\mathbb R^n)$, one expects the
solution to belong to $H^1(0,T;\mathbb R^n)$.

Define

$$
Y:=H^1(0,T;\mathbb R^n).
$$

Recall that in finite dimensions

$$
H^1(0,T;\mathbb R^n)\hookrightarrow C([0,T];\mathbb R^n),
$$

hence the terminal value $y(T)$ is meaningful.

**Proposition.**
For every $u\in U$, the state equation admits a unique solution

$$
y\in Y=H^1(0,T;\mathbb R^n).
$$

Moreover, the solution is given by the variation-of-constants formula

$$
y(t)=e^{tA}y_0+\int_0^t e^{(t-s)A}(Bu(s)+f(s))\,ds,
$$

and there exists a constant $C>0$, depending only on $A$, $B$, and $T$, such that

$$
\|y\|_{H^1(0,T;\mathbb R^n)} + \|y\|_{C([0,T];\mathbb R^n)}
\le C\bigl(|y_0|+\|u\|_{L^2(0,T;\mathbb R^m)}+\|f\|_{L^2(0,T;\mathbb R^n)}\bigr).
$$

**Proof.**
Fix $u\in U$ and define

$$
g(t):=Bu(t)+f(t)\in L^2(0,T;\mathbb R^n).
$$

Since $A$ is a constant matrix, the matrix exponential $e^{tA}$ is well defined and continuous.
The formula

$$
y(t):=e^{tA}y_0+\int_0^t e^{(t-s)A}g(s)\,ds
$$

defines a continuous function on $[0,T]$.
Differentiating under the integral sign yields

$$
\dot y(t)=Ae^{tA}y_0+g(t)+\int_0^t A e^{(t-s)A} g(s)\,ds
=Ay(t)+g(t)
$$

for a.e. $t\in(0,T)$, and clearly $y(0)=y_0$.
Thus $y\in H^1(0,T;\mathbb R^n)$ and solves the ODE.

Uniqueness follows from linearity.
Indeed, if $y_1$ and $y_2$ solve the same problem, then $w:=y_1-y_2$ satisfies

$$
\dot w = A w,
\qquad
w(0)=0.
$$

Hence $w(t)=e^{tA}w(0)=0$ for all $t$.

To estimate $y$, let

$$
M_A:=\max_{0\le t\le T} \|e^{tA}\|.
$$

Then

$$
|y(t)|
\le M_A |y_0| + M_A \int_0^T |g(s)|\,ds
\le M_A |y_0| + M_A T^{1/2}\|g\|_{L^2(0,T)}.
$$

Taking the supremum in $t$ gives

$$
\|y\|_{C([0,T])}
\le C\bigl(|y_0|+\|g\|_{L^2(0,T)}\bigr)
\le C\bigl(|y_0|+\|u\|_{L^2(0,T)}+\|f\|_{L^2(0,T)}\bigr).
$$

Since $\dot y=Ay+g$,

$$
\|\dot y\|_{L^2(0,T)}
\le \|A\|\,\|y\|_{L^2(0,T)} + \|g\|_{L^2(0,T)}
\le C\bigl(|y_0|+\|u\|_{L^2(0,T)}+\|f\|_{L^2(0,T)}\bigr).
$$

Combining the bounds for $y$ and $\dot y$ yields the claim. $\square$

---

The proposition allows us to define the control-to-state map

$$
S:U\to Y,
\qquad
u\mapsto y=S(u).
$$

Because the state equation is linear, $S$ is linear and continuous.
The reduced functional is therefore

$$
j(u)=J(S(u),u).
$$

---

## Existence and Uniqueness of an Optimal Control

The previous lectures already established an abstract existence theorem in
Hilbert spaces.
In the present setting the argument becomes particularly transparent.

**Theorem.**
Assume $\alpha>0$ and let $U_{\mathrm{ad}}\subset U$ be nonempty, closed, and convex.
Then the reduced problem

$$
\min_{u\in U_{\mathrm{ad}}} j(u)
$$

admits a unique solution $\bar u\in U_{\mathrm{ad}}$.

**Proof.**
Since $S:U\to Y$ is continuous and $J$ is continuous on $Y\times U$, the map $j:U\to \mathbb R$
 is continuous.
Moreover,

$$
j(u)
=
\frac12\|S(u)-y_d\|_{L^2(0,T;\mathbb R^n)}^2
+\frac\alpha2\|u\|_{L^2(0,T;\mathbb R^m)}^2
+\frac\beta2 |S(u)(T)-y_T|^2
\ge \frac\alpha2\|u\|_U^2.
$$

Hence $j$ is coercive on $U$.
Let $(u_k)$ be a minimizing sequence in $U_{\mathrm{ad}}$.
Coercivity implies that $(u_k)$ is bounded in $U$.
Since $U$ is a Hilbert space, there exists a subsequence, not relabeled, and an element
$\bar u\in U$ such that

$$
u_k \rightharpoonup \bar u
\qquad \text{weakly in }U.
$$

Because $U_{\mathrm{ad}}$ is closed and convex, it is weakly closed, hence $\bar u\in U_{\mathrm{ad}}$.
By continuity and linearity of $S$,

$$
S(u_k) \rightharpoonup S(\bar u)
\qquad \text{weakly in }Y.
$$

Since the trace map

$$
y\mapsto y(T)
$$

is a continuous linear functional on $H^1(0,T;\mathbb R^n)$, we also have

$$
S(u_k)(T) \rightharpoonup S(\bar u)(T)
\qquad \text{in } \mathbb R^n.
$$

Since the norm is weakly lower semicontinuous in Hilbert spaces, each term in $j$ is weakly lower semicontinuous, hence

$$
j(\bar u)
\le \liminf_{k\to\infty} j(u_k),
$$

so $\bar u$ is a minimizer.

To prove uniqueness, let $u_1\neq u_2$ and $\theta\in(0,1)$.
Since $S$ is linear,

$$
S(\theta u_1+(1-\theta)u_2)=\theta S(u_1)+(1-\theta)S(u_2).
$$

The first and third terms of $j$ are convex, while the control term is strictly convex because $\alpha>0$:

$$
\|\theta u_1+(1-\theta)u_2\|_U^2
< \theta \|u_1\|_U^2 + (1-\theta)\|u_2\|_U^2.
$$

Therefore $j$ is strictly convex and the minimizer is unique. $\square$

---

## Linearized State Equation

To differentiate the reduced cost, we perturb the control by $h\in U$.
Since the state equation is linear, the state increment is described by a linear ODE
independent of the base point.

Let $u\in U$ and $h\in U$.
Define

$$
z_h := S'(u)h.
$$

Then $z_h$ solves

$$
\dot z_h(t)=A z_h(t)+B h(t)
\qquad \text{for a.e. } t\in(0,T),
\qquad
z_h(0)=0.
$$

Because the control-to-state map is linear, in fact

$$
S'(u)h = S(h)-S(0)
$$

for every $u\in U$.
Equivalently, one may derive the linearized equation directly from

$$
\dot y = Ay+Bu+f,
$$

by replacing $u$ with $u+\varepsilon h$, subtracting the equation for $u$, dividing by $\varepsilon$,
and passing to the limit.

The directional derivative of the reduced functional therefore reads

$$
j'(u)h
=
\int_0^T (y(t)-y_d(t))\cdot z_h(t)\,dt
+\beta (y(T)-y_T)\cdot z_h(T)
+\alpha \int_0^T u(t)\cdot h(t)\,dt.
$$

The difficulty is that this formula still contains the linearized state $z_h$.
As in the elliptic case, the adjoint variable removes this dependence.

---

## Adjoint Equation

The adjoint equation is the backward-in-time equation dual to the linearized
state equation.
Its role is exactly the same as in the elliptic case, but the time direction is reversed.

Let $u\in U$ be fixed, and let $y=S(u)$.
We define the adjoint state $p\in H^1(0,T;\mathbb R^n)$ by

$$
-\dot p(t)=A^T p(t)+y(t)-y_d(t)
\qquad \text{for a.e. } t\in(0,T),
\qquad
p(T)=\beta\,(y(T)-y_T).
$$

This is again a linear ODE, now with terminal condition prescribed at $t=T$.
Its unique solution is obtained by solving backward in time, or equivalently by the change of variable
$s=T-t$.

**Proposition.**
Let $u\in U$, $y=S(u)$, and let $p$ solve the adjoint equation above.
Then for every $h\in U$, with $z_h$ solving the linearized state equation,

$$
\int_0^T (y(t)-y_d(t))\cdot z_h(t)\,dt
+\beta (y(T)-y_T)\cdot z_h(T)
=
\int_0^T B^T p(t)\cdot h(t)\,dt.
$$

**Proof.**
From the linearized state equation,

$$
\dot z_h - A z_h - B h = 0.
$$

Multiply by $p$ and integrate over $(0,T)$:

$$
\int_0^T p(t)\cdot \bigl(\dot z_h(t)-A z_h(t)-B h(t)\bigr)\,dt = 0.
$$

Using $(p,Az_h)=(A^Tp,z_h)$ and integrating by parts in time,

$$
\int_0^T p\cdot \dot z_h\,dt
= p(T)\cdot z_h(T)-p(0)\cdot z_h(0)-\int_0^T \dot p\cdot z_h\,dt.
$$

Since $z_h(0)=0$,

$$
p(T)\cdot z_h(T)-\int_0^T \bigl(\dot p + A^T p\bigr)\cdot z_h\,dt
-\int_0^T B^T p\cdot h\,dt = 0.
$$

By the adjoint equation,

$$
-(\dot p + A^T p)=y-y_d,
$$

hence

$$
p(T)\cdot z_h(T)+\int_0^T (y-y_d)\cdot z_h\,dt
-\int_0^T B^T p\cdot h\,dt = 0.
$$

Finally, the terminal condition gives

$$
p(T)\cdot z_h(T)=\beta (y(T)-y_T)\cdot z_h(T),
$$

which yields the claim. $\square$

---

## Reduced Gradient Formula

We can now eliminate the linearized state from the derivative.

**Theorem.**
The reduced functional $j:U\to \mathbb R$ is Fréchet differentiable and, for every $u\in U$,

$$
j'(u)h = \int_0^T \bigl(\alpha u(t)+B^T p(t)\bigr)\cdot h(t)\,dt
\qquad \forall h\in U,
$$

where $p$ is the adjoint associated with $u$.
Hence the gradient of $j$ in the Hilbert space $U=L^2(0,T;\mathbb R^m)$ is

$$
\nabla j(u)=\alpha u + B^T p.
$$

**Proof.**
We already computed

$$
j'(u)h
=
\int_0^T (y-y_d)\cdot z_h\,dt
+\beta (y(T)-y_T)\cdot z_h(T)
+\alpha\int_0^T u\cdot h\,dt.
$$

The adjoint identity from the previous proposition gives

$$
\int_0^T (y-y_d)\cdot z_h\,dt
+\beta (y(T)-y_T)\cdot z_h(T)
=
\int_0^T B^T p\cdot h\,dt.
$$

Therefore

$$
j'(u)h = \int_0^T (\alpha u + B^T p)\cdot h\,dt.
$$

Since this is a bounded linear functional of $h$, the Fréchet derivative is represented in $U$
by the function $\alpha u + B^T p$. $\square$

---

## First-Order Optimality System

We can now write the necessary and sufficient optimality conditions.
Because the problem is convex and the reduced functional is strictly convex,
first-order conditions characterize the unique minimizer.

### Unconstrained case

If

$$
U_{\mathrm{ad}}=U=L^2(0,T;\mathbb R^m),
$$

then the stationarity condition is simply

$$
\nabla j(\bar u)=0,
$$

i.e.

$$
\alpha \bar u + B^T \bar p = 0
\qquad \text{in }L^2(0,T;\mathbb R^m).
$$

Thus the optimal control satisfies the explicit relation

$$
\bar u(t) = -\frac1\alpha B^T \bar p(t)
\qquad \text{for a.e. } t\in(0,T).
$$

The optimality system becomes

$$
\begin{cases}
\dot{\bar y}(t)=A\bar y(t)+B\bar u(t)+f(t),
& \text{a.e. } t\in(0,T),\\
\bar y(0)=y_0,\\
-\dot{\bar p}(t)=A^T\bar p(t)+\bar y(t)-y_d(t),
& \text{a.e. } t\in(0,T),\\
\bar p(T)=\beta(\bar y(T)-y_T),\\
\alpha \bar u(t)+B^T\bar p(t)=0,
& \text{a.e. } t\in(0,T).
\end{cases}
$$

This is a **forward-backward system**:

- the state evolves forward from $t=0$;
- the adjoint evolves backward from $t=T$.

This two-sided time structure is the first major structural difference from elliptic control.

### Constrained case

Assume now that $U_{\mathrm{ad}}\subset U$ is closed and convex.
Then $\bar u\in U_{\mathrm{ad}}$ is optimal if and only if

$$
j'(\bar u)(u-\bar u)\ge 0
\qquad \forall u\in U_{\mathrm{ad}}.
$$

Using the adjoint representation of the derivative, this becomes

$$
\int_0^T \bigl(\alpha \bar u(t)+B^T\bar p(t)\bigr)\cdot \bigl(u(t)-\bar u(t)\bigr)\,dt
\ge 0
\qquad \forall u\in U_{\mathrm{ad}}.
$$

Equivalently, in normal-cone form,

$$
0\in \alpha \bar u + B^T \bar p + N_{U_{\mathrm{ad}}}(\bar u)
\qquad \text{in }U.
$$

Hence the full optimality system is

$$
\begin{cases}
\dot{\bar y}=A\bar y+B\bar u+f,
\qquad \bar y(0)=y_0,\\[0.3em]
-\dot{\bar p}=A^T\bar p+\bar y-y_d,
\qquad \bar p(T)=\beta(\bar y(T)-y_T),\\[0.3em]
0\in \alpha \bar u + B^T \bar p + N_{U_{\mathrm{ad}}}(\bar u).
\end{cases}
$$

---

## Box Constraints and Projection Formula

A particularly important case is the box-constrained set

$$
U_{\mathrm{ad}}
:=
\left\{
 u\in L^2(0,T;\mathbb R^m):
 u_a(t)\le u(t)\le u_b(t)
 \text{ for a.e. } t\in(0,T)
\right\},
$$

where $u_a,u_b\in L^\infty(0,T;\mathbb R^m)$ and the inequalities are understood componentwise.

In this case the optimality condition is equivalent to the pointwise projection formula

$$
\bar u(t)=P_{[u_a(t),u_b(t)]}\!\left(-\frac1\alpha B^T\bar p(t)\right)
\qquad \text{for a.e. } t\in(0,T),
$$

where for $\xi\in \mathbb R^m$

$$
P_{[u_a(t),u_b(t)]}(\xi)
=
\min\bigl(\max(\xi,u_a(t)),u_b(t)\bigr)
$$

componentwise.

The proof is identical in structure to the elliptic case:

- the variational inequality is posed in the Hilbert space $L^2(0,T;\mathbb R^m)$;
- the admissible set is a closed convex box;
- the projection is characterized pointwise almost everywhere.

Thus the only genuinely new analytical ingredient introduced by time dependence is not the
control condition, but the forward-backward evolution structure of state and adjoint.

---

## Forward-Backward Interpretation

It is worth isolating the conceptual meaning of the adjoint in the dynamical case.

For a fixed control $u$:

- the state equation propagates the effect of the control from time $0$ to time $T$;
- the adjoint equation propagates sensitivity information from time $T$ backward to time $0$.

The terminal condition

$$
p(T)=\beta(y(T)-y_T)
$$

encodes the derivative of the terminal observation.
If $\beta=0$, then $p(T)=0$ and only the distributed tracking term drives the adjoint.
If instead the cost is purely terminal,

$$
J(y,u)=\frac\alpha2\int_0^T |u(t)|^2\,dt + \frac\beta2 |y(T)-y_T|^2,
$$

then the adjoint satisfies

$$
-\dot p = A^T p,
\qquad
p(T)=\beta(y(T)-y_T).
$$

This is exactly the same phenomenon that one encounters later for parabolic PDEs:

- state equation forward in time;
- adjoint equation backward in time;
- control condition coupling state and adjoint at the same time level.

---

## From ODEs to Evolution Equations

The ODE model can be rewritten abstractly as

$$
y_t + \mathcal A y = \mathcal B u + f,
$$

where now the state $y(t)$ is an element of the finite-dimensional Hilbert space
$H=\mathbb R^n$ and

$$
\mathcal A := -A,
\qquad
\mathcal B := B.
$$

For parabolic PDEs, the same formula remains formally correct, but the state at each time is no longer a vector in $\mathbb R^n$.
Instead:

- $y(t)$ is a function of the space variable, e.g. $y(t)\in H_0^1(\Omega)$ or $L^2(\Omega)$;
- the derivative $y_t$ generally belongs to a dual space rather than to the same space as $y$;
- time integration must be carried out for vector-valued functions.

This is the reason why the usual scalar-valued Lebesgue and Sobolev spaces are not sufficient.
We need function spaces of the form

$$
L^p(0,T;X),
$$

where $X$ is itself a Banach or Hilbert space.
These are the **Bochner spaces**.

---

## Strongly Measurable Vector-Valued Functions

Let $X$ be a Banach space.
A function

$$
y:(0,T)\to X
$$

is called **simple** if it has the form

$$
y(t)=\sum_{k=1}^N x_k\,\chi_{E_k}(t),
$$

where $x_k\in X$ and $E_k\subset(0,T)$ are measurable sets.

A function $y:(0,T)\to X$ is called **strongly measurable** if there exists a sequence of
simple functions $(y_n)$ such that

$$
y_n(t)\to y(t)
\qquad \text{for a.e. } t\in(0,T).
$$

This is the natural notion of measurability for vector-valued functions.
In the Hilbert spaces used in parabolic PDEs, separability holds, so this definition behaves well.

If $y$ is strongly measurable and

$$
\int_0^T \|y(t)\|_X\,dt < \infty,
$$

then $y$ is **Bochner integrable** and one may define

$$
\int_0^T y(t)\,dt \in X
$$

as the limit of the integrals of simple approximations.
This is the vector-valued analogue of the usual Lebesgue integral.

---

## Bochner Spaces $L^p(0,T;X)$

Let $1\le p<\infty$.
We define

$$
L^p(0,T;X)
:=
\left\{
 y:(0,T)\to X:
 y \text{ strongly measurable and }
 \int_0^T \|y(t)\|_X^p\,dt < \infty
\right\}.
$$

The norm is

$$
\|y\|_{L^p(0,T;X)}
:=
\left(\int_0^T \|y(t)\|_X^p\,dt\right)^{1/p}.
$$

Similarly,

$$
L^\infty(0,T;X)
:=
\left\{
 y:(0,T)\to X:
 y \text{ strongly measurable and }
 \operatorname*{ess\,sup}_{t\in(0,T)} \|y(t)\|_X < \infty
\right\}.
$$

Standard facts:

- if $X$ is Banach, then $L^p(0,T;X)$ is Banach;
- if $X$ is Hilbert and $p=2$, then $L^2(0,T;X)$ is Hilbert with inner product

  $$
  (y,z)_{L^2(0,T;X)} := \int_0^T (y(t),z(t))_X\,dt;
  $$

- if $X\hookrightarrow Z$ continuously, then

  $$
  L^p(0,T;X)\hookrightarrow L^p(0,T;Z)
  $$

  continuously.

### Fundamental examples

Let $\Omega\subset \mathbb R^d$ be a bounded Lipschitz domain and set

$$
Q_T:=\Omega\times(0,T).
$$

Then:

1. if $X=L^2(\Omega)$,

   $$
   L^2(0,T;L^2(\Omega)) \cong L^2(Q_T);
   $$

2. if $X=H_0^1(\Omega)$,

   $$
   L^2(0,T;H_0^1(\Omega))
   $$

   consists of functions square integrable in time with values in $H_0^1(\Omega)$;

3. if $X=H^{-1}(\Omega)$,

   $$
   L^2(0,T;H^{-1}(\Omega))
   $$

   is the natural space for weak time derivatives of parabolic states.

Thus, in parabolic theory, the same function $y(x,t)$ is viewed as a map

$$
t\mapsto y(t):=y(\cdot,t)
$$

with values in a spatial function space.

---

## Weak Time Derivatives

For parabolic equations, the time derivative does not usually belong to the same space as the state.
This forces a dual-space formulation.

Let $X$ be a Banach space and let

$$
y\in L^1(0,T;X).
$$

A function

$$
z\in L^1(0,T;X)
$$

is called the **weak time derivative** of $y$ if for every scalar test function
$\varphi\in C_c^\infty(0,T)$ and every $\ell\in X'$,

$$
\int_0^T \langle \ell,y(t)\rangle_{X',X}\, \varphi'(t)\,dt
=
-\int_0^T \langle \ell,z(t)\rangle_{X',X}\, \varphi(t)\,dt.
$$

In this case we write

$$
z = y_t.
$$

If $X$ is Hilbert, one can identify $X\cong X'$ by the Riesz map and recover the familiar scalar definition.

The Sobolev space of $X$-valued functions is then

$$
H^1(0,T;X)
:=
\{y\in L^2(0,T;X): y_t\in L^2(0,T;X)\}.
$$

For ODEs, where $X=\mathbb R^n$, this is the space used for the state and adjoint variables.

For parabolic PDEs, however, the correct setting is generally not $H^1(0,T;X)$ with a single space $X$, but a mixed space involving a Hilbert triple.

---

## Gelfand Triples and the Space $W(0,T)$

Let $V$ and $H$ be Hilbert spaces such that

$$
V \hookrightarrow H
$$

continuously and densely.
By identifying $H$ with its dual $H'$ through the Riesz isomorphism, one obtains the **Gelfand triple**

$$
V \hookrightarrow H \cong H' \hookrightarrow V'.
$$

The last embedding is defined by

$$
\langle h,v\rangle_{V',V} := (h,v)_H
\qquad \forall h\in H,\ \forall v\in V.
$$

The canonical parabolic example is

$$
V=H_0^1(\Omega),
\qquad
H=L^2(\Omega),
\qquad
V'=H^{-1}(\Omega).
$$

The natural energy space for parabolic problems is

$$
W(0,T)
:=
\{y\in L^2(0,T;V): y_t\in L^2(0,T;V')\}.
$$

It is a Hilbert space with norm

$$
\|y\|_{W(0,T)}^2
:=
\|y\|_{L^2(0,T;V)}^2 + \|y_t\|_{L^2(0,T;V')}^2.
$$

This is the parabolic analogue of $H^1(0,T;\mathbb R^n)$ for ODEs.

The key point is the asymmetry:

- the state lives in $V$ with respect to the space variable;
- the time derivative lives in the dual space $V'$.

This is forced by the weak formulation of the PDE.
For the heat equation,

$$
y_t - \Delta y = F,
$$

one expects

$$
y(t)\in H_0^1(\Omega),
\qquad
y_t(t)\in H^{-1}(\Omega).
$$

---

## Fundamental Theorem for $W(0,T)$

The space $W(0,T)$ has a decisive property: its elements possess a continuous representative with values in $H$.
This is what makes initial and terminal conditions meaningful.

**Theorem (Lions-Magenes).**
Let

$$
V \hookrightarrow H \hookrightarrow V'
$$

be a Gelfand triple.
Then:

1. every $y\in W(0,T)$ admits a representative, still denoted by $y$, such that

   $$
   y\in C([0,T];H);
   $$

2. the embedding

   $$
   W(0,T) \hookrightarrow C([0,T];H)
   $$

   is continuous;

3. if $y,v\in W(0,T)$, then the scalar map

   $$
   t\mapsto (y(t),v(t))_H
   $$

   is absolutely continuous and satisfies

   $$
   \frac{d}{dt}(y(t),v(t))_H
   =
   \langle y_t(t),v(t)\rangle_{V',V}
   +
   \langle v_t(t),y(t)\rangle_{V',V}
   $$

   for a.e. $t\in(0,T)$.

In particular, taking $v=y$ gives

$$
\frac12\frac{d}{dt}\|y(t)\|_H^2
=
\langle y_t(t),y(t)\rangle_{V',V}
\qquad \text{for a.e. } t\in(0,T).
$$

Integrating between $s$ and $t$ yields the energy identity

$$
\frac12\|y(t)\|_H^2 - \frac12\|y(s)\|_H^2
=
\int_s^t \langle y_t(\tau),y(\tau)\rangle_{V',V}\,d\tau.
$$

More generally, for $y,v\in W(0,T)$,

$$
(y(t),v(t))_H - (y(s),v(s))_H
=
\int_s^t \langle y_t(\tau),v(\tau)\rangle_{V',V}\,d\tau
+
\int_s^t \langle v_t(\tau),y(\tau)\rangle_{V',V}\,d\tau.
$$

This is the time-integration-by-parts formula needed for parabolic adjoints.

---

## Abstract Parabolic Problem

With the previous tools in place, one can formulate the prototype parabolic state equation.

Let $V\hookrightarrow H\hookrightarrow V'$ be a Gelfand triple.
Let

$$
a:V\times V\to \mathbb R
$$

be a bilinear form satisfying:

1. **continuity**:

   $$
   |a(w,v)|\le M\|w\|_V\|v\|_V
   \qquad \forall w,v\in V;
   $$

2. **coercivity**:

   $$
   a(v,v)\ge c_a \|v\|_V^2
   \qquad \forall v\in V,
   $$

   for some $c_a>0$.

Define the operator

$$
A:V\to V',
\qquad
\langle Ay,v\rangle_{V',V}=a(y,v).
$$

Given

$$
F\in L^2(0,T;V'),
\qquad
y_0\in H,
$$

the weak parabolic problem is:
find $y\in W(0,T)$ such that

$$
\langle y_t(t),v\rangle_{V',V} + a(y(t),v)
= \langle F(t),v\rangle_{V',V}
\qquad \forall v\in V,
\quad \text{for a.e. } t\in(0,T),
$$

with

$$
y(0)=y_0 \quad \text{in } H.
$$

**Theorem.**
Under the assumptions above, the parabolic problem admits a unique solution

$$
y\in W(0,T).
$$

Moreover,

$$
\|y\|_{L^2(0,T;V)}
+
\|y\|_{L^\infty(0,T;H)}
+
\|y_t\|_{L^2(0,T;V')}
\le
C\bigl(\|F\|_{L^2(0,T;V')}+\|y_0\|_H\bigr),
$$

for a constant $C$ depending only on the continuity and coercivity constants and on $T$.

This theorem is the infinite-dimensional analogue of the well-posedness proposition for the ODE state equation.
The analogies are exact:

- matrix $A$ becomes an operator $A:V\to V'$;
- state space $H^1(0,T;\mathbb R^n)$ becomes $W(0,T)$;
- the terminal value is meaningful because $W(0,T)\hookrightarrow C([0,T];H)$.

---

## Heat Equation as Canonical Example

Take

$$
V=H_0^1(\Omega),
\qquad
H=L^2(\Omega),
\qquad
V'=H^{-1}(\Omega),
$$

and define

$$
a(y,v)=\int_\Omega \nabla y\cdot \nabla v\,dx.
$$

Then

$$
\langle Ay,v\rangle = \int_\Omega \nabla y\cdot \nabla v\,dx
$$

corresponds to the operator $A=-\Delta$ in weak form.
Given a control $u\in L^2(0,T;L^2(\Omega))$, one may write

$$
F(t)=u(t)+f(t) \in H \hookrightarrow V'.
$$

The parabolic state equation becomes

$$
\langle y_t(t),v\rangle_{H^{-1},H_0^1}
+
\int_\Omega \nabla y(t)\cdot \nabla v\,dx
=
\int_\Omega (u(t)+f(t))v\,dx
$$

for all $v\in H_0^1(\Omega)$ and a.e. $t\in(0,T)$, with $y(0)=y_0$.

This is the standard weak formulation of the heat equation

$$
y_t - \Delta y = u+f
\qquad \text{in } \Omega\times(0,T),
$$

with homogeneous Dirichlet boundary condition.

---

## What Changes in the Optimality System for Parabolic PDEs?

At the formal level, almost nothing changes.
The parabolic optimal control problem has the structure

$$
\min_{u\in U_{\mathrm{ad}}}
\frac12\int_0^T \|y(t)-y_d(t)\|_H^2\,dt
+
\frac\alpha2 \|u\|_{L^2(0,T;U)}^2
+
\frac\beta2 \|y(T)-y_T\|_H^2
$$

subject to

$$
y_t + Ay = Bu + f,
\qquad
y(0)=y_0.
$$

The optimality system has the same forward-backward pattern as in the ODE case:

- state equation forward in time;
- adjoint equation backward in time;
- control equation or variational inequality at each time.

The new difficulty is not conceptual but functional-analytic:

- the state is an element of $W(0,T)$, not of $H^1(0,T;\mathbb R^n)$;
- the time derivative lives in $V'$;
- the integration-by-parts identity must be understood in the Gelfand triple.

---

## Summary

This lecture introduced the time-dependent side of optimal control in two layers.

### ODE layer

For the linear-quadratic control problem governed by

$$
\dot y = Ay+Bu+f,
\qquad
y(0)=y_0,
$$

we proved:

- well-posedness of the state equation in $H^1(0,T;\mathbb R^n)$;
- existence and uniqueness of the optimal control;
- the adjoint equation

  $$
  -\dot p = A^T p + y - y_d,
  \qquad
  p(T)=\beta(y(T)-y_T);
  $$

- the gradient formula

  $$
  \nabla j(u)=\alpha u + B^T p;
  $$

- the first-order optimality condition

  $$
  \int_0^T (\alpha \bar u + B^T \bar p)\cdot (u-\bar u)\,dt\ge 0.
  $$

Hence time-dependent optimality already appears as a forward-backward system.

### Functional-analytic layer

To pass from ODEs to parabolic PDEs we introduced:

- strongly measurable vector-valued functions;
- Bochner spaces $L^p(0,T;X)$;
- weak time derivatives in dual spaces;
- Gelfand triples $V\hookrightarrow H\hookrightarrow V'$;
- the energy space

  $$
  W(0,T)=\{y\in L^2(0,T;V): y_t\in L^2(0,T;V')\};
  $$

- the embedding

  $$
  W(0,T)\hookrightarrow C([0,T];H),
  $$

  and the corresponding integration-by-parts identity.

These are exactly the tools needed for the next lecture, where the same adjoint-based
optimality machinery will be applied to parabolic PDEs.

---

## References

- F. Tröltzsch, *Optimal Control of Partial Differential Equations*, AMS, 2010.
- J. C. De los Reyes, *Numerical PDE-Constrained Optimization*, Springer, 2015.
- A. Manzoni, A. Quarteroni, S. Salsa, *Optimal Control of Partial Differential Equations*, Springer, 2021.
- J.-L. Lions, *Optimal Control of Systems Governed by Partial Differential Equations*, Springer, 1971.
- J.-L. Lions, E. Magenes, *Non-Homogeneous Boundary Value Problems and Applications*, Springer, 1972.

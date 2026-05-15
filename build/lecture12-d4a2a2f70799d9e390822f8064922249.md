# Time Discretization for Parabolic Optimal Control Problems

## Overview

The previous lectures introduced the continuous optimality system for
time-dependent control problems and the functional-analytic setting needed for
parabolic PDEs:

- Gelfand triples $V\hookrightarrow H\hookrightarrow V'$;
- the energy space
  $$
  W(0,T):=\{y\in L^2(0,T;V): y_t\in L^2(0,T;V')\};
  $$
- forward state equations and backward adjoint equations;
- variational inequalities and projection formulas for box constraints.

This lecture addresses the next numerical question:

> what is the correct discrete optimality system after time discretization?

We focus on one time-stepping method, implicit Euler, and derive the
corresponding discrete optimality system carefully.  The goal is not to cover
all possible time discretizations, but to understand the structure that any
reliable implementation must preserve.

The logical path is:

1. recall the continuous linear-quadratic parabolic control problem;
2. discretize the state equation in time with implicit Euler;
3. define the discrete reduced optimization problem;
4. derive the discrete adjoint equation from a discrete Lagrangian;
5. identify the discrete gradient and the first-order optimality system;
6. compare discretize-then-optimize with optimize-then-discretize;
7. write the fully discrete finite element matrix form;
8. add box constraints through a time-discrete variational inequality.

The central message is:

> the discrete adjoint is not obtained by informally reversing time; it is the
> transpose, with respect to the chosen discrete inner products, of the
> linearized discrete state equation.

This is the principle behind adjoint-consistent time discretizations.

---

## Continuous Problem

Let

$$
V\hookrightarrow H\hookrightarrow V'
$$

be a Gelfand triple, where $V$ and $H$ are Hilbert spaces, the embedding
$V\hookrightarrow H$ is dense and continuous, and $H$ is identified with its
dual.

For the heat equation with homogeneous Dirichlet boundary conditions one may
keep in mind

$$
V=H_0^1(\Omega),
\qquad
H=L^2(\Omega),
\qquad
V'=H^{-1}(\Omega).
$$

Let $U$ be a Hilbert space for the control values, and assume that

$$
B\in \mathcal L(U,V')
$$

is a bounded control operator.  Let

$$
f\in L^2(0,T;V'),
\qquad
y_0\in H,
\qquad
y_d\in L^2(0,T;H),
$$

and let $\alpha>0$.

We consider the state equation

$$
y_t(t)+Ay(t)=Bu(t)+f(t)
\qquad\text{in }V'\text{ for a.e. }t\in(0,T),
$$

with initial condition

$$
y(0)=y_0.
$$

Here $A:V\to V'$ is induced by a bilinear form

$$
a(\cdot,\cdot):V\times V\to \mathbb R,
$$

that is continuous and coercive:

$$
|a(v,w)|\le M\|v\|_V\|w\|_V,
\qquad
a(v,v)\ge \gamma \|v\|_V^2
\qquad \forall v,w\in V,
$$

for constants $M,\gamma>0$.

The weak form of the state equation is:

find $y\in W(0,T)$ such that $y(0)=y_0$ and

$$
\langle y_t(t),v\rangle_{V',V}
+a(y(t),v)
=
\langle Bu(t)+f(t),v\rangle_{V',V}
\qquad \forall v\in V
$$

for almost every $t\in(0,T)$.

The cost functional is

$$
J(y,u)
=
\frac12\int_0^T \|y(t)-y_d(t)\|_H^2\,dt
+\frac\alpha2\int_0^T \|u(t)\|_U^2\,dt + \frac\beta2 \| y(T) - y_T\|_H^2.
$$

The unconstrained optimal control problem is

$$
\min_{u\in L^2(0,T;U)} j(u):=J(y(u),u),
$$

where $y(u)$ denotes the state associated with $u$.

For constrained problems we replace the full control space by a nonempty
closed convex set

$$
U_{\mathrm{ad}}\subset L^2(0,T;U).
$$

---

## Continuous Optimality System

The continuous first-order system has already been derived in the previous
lecture.  We recall it only to fix the notation.

For a given control $u$, the state $y=y(u)$ solves

$$
\begin{cases}
y_t+Ay=Bu+f,\\
y(0)=y_0.
\end{cases}
$$

The corresponding adjoint $p$ solves the backward parabolic problem

$$
\begin{cases}
-p_t+A^*p=M(y-y_d),\\
p(T)=\beta(y(T)-y_T).
\end{cases}
$$

In weak form this means

$$
-\langle p_t(t),v\rangle_{V',V}
+a(v,p(t))
=
(y(t)-y_d(t),v)_H
\qquad \forall v\in V.
$$

Here $A^*$ is the adjoint operator.  If $a$ is symmetric, then $A^*=A$.

The reduced gradient is

$$
\nabla j(u)=\alpha u+B^*p
\qquad\text{in }L^2(0,T;U),
$$

with the usual interpretation of $B^*$ through duality:

$$
(B^*p,w)_U = \langle Bw,p\rangle_{V',V}.
$$

Thus, in the unconstrained case, the optimality system is

$$
\begin{cases}
\bar y_t+A\bar y=B\bar u+f,\\
\bar y(0)=y_0,\\
-\bar p_t+A^*\bar p=M(\bar y-y_d),\\
\bar p(T)=\beta(y(T)-y_T),\\
\alpha \bar u+B^*\bar p=0.
\end{cases}
$$

For a closed convex admissible set $U_{\mathrm{ad}}$, the last equation is
replaced by the variational inequality

$$
\int_0^T
(\alpha \bar u(t)+B^*\bar p(t), u(t)-\bar u(t))_U\,dt
\ge 0
\qquad \forall u\in U_{\mathrm{ad}}.
$$

The purpose of this lecture is to build the discrete analogue of this system
from the discrete optimization problem itself.

---

## Implicit Euler Discretization of the State

Let

$$
0=t_0<t_1<\cdots<t_N=T,
\qquad
\tau:=t_k-t_{k-1}=\frac{T}{N}.
$$

We use one control value $u^k\in U$ on each interval

$$
I_k:=(t_{k-1},t_k],
\qquad k=1,\ldots,N.
$$

The data are approximated by values

$$
f^k\in V',
\qquad
y_d^k\in H.
$$

For instance, one may take time averages

$$
f^k:=\frac1\tau\int_{t_{k-1}}^{t_k} f(t)\,dt,
\qquad
y_d^k:=\frac1\tau\int_{t_{k-1}}^{t_k} y_d(t)\,dt.
$$

The implicit Euler approximation of the state equation is:

find $y^1,\ldots,y^N\in V$, with $y^0=y_0$, such that

$$
\left(\frac{y^k-y^{k-1}}{\tau},v\right)_H
+a(y^k,v)
=
\langle Bu^k+f^k,v\rangle_{V',V}
\qquad \forall v\in V,
$$

for $k=1,\ldots,N$.

Equivalently,

$$
(y^k,v)_H+\tau a(y^k,v)
=
(y^{k-1},v)_H
+\tau\langle Bu^k+f^k,v\rangle_{V',V}.
$$

This is a forward time-stepping scheme.  Once $y^{k-1}$ and $u^k$ are known,
the equation for $y^k$ is an elliptic problem with bilinear form

$$
(y,v)_H+\tau a(y,v).
$$

Coercivity follows immediately:

$$
(v,v)_H+\tau a(v,v)
\ge
\|v\|_H^2+\tau\gamma\|v\|_V^2.
$$

Therefore each time step is well posed by the Lax-Milgram theorem.

**Proposition.**
For every sequence of controls

$$
u_\tau=(u^1,\ldots,u^N)\in U^N
$$

there exists a unique discrete state sequence

$$
y_\tau=(y^1,\ldots,y^N)\in V^N.
$$

Moreover, the map

$$
S_\tau:U^N\to V^N,
\qquad
u_\tau\mapsto y_\tau,
$$

is affine and continuous.

**Proof.**
For fixed $k$, assume $y^{k-1}$ is already known.  The left-hand side is the
coercive bilinear form

$$
m_\tau(y,v):=(y,v)_H+\tau a(y,v),
$$

and the right-hand side

$$
\ell_k(v):=(y^{k-1},v)_H
+\tau\langle Bu^k+f^k,v\rangle_{V',V}
$$

is a bounded linear functional on $V$.  Hence $y^k$ exists and is unique.
Induction over $k=1,\ldots,N$ gives the full discrete state.

Linearity in the variables $u^1,\ldots,u^N$ and affine dependence on
$(y_0,f^1,\ldots,f^N)$ follow directly from the recursion.  Continuity follows
from the stability estimate obtained by applying Lax-Milgram at each step and
iterating over finitely many time levels. $\square$

---

## The Discrete Optimization Problem

The implicit Euler state equation defines a finite-dimensional-in-time
optimization problem.

The discrete cost functional is

$$
J_\tau(y_\tau,u_\tau)
:=
\frac{\tau}{2}\sum_{k=1}^N \|y^k-y_d^k\|_H^2
+\frac{\alpha\tau}{2}\sum_{k=1}^N \|u^k\|_U^2.
$$

The reduced discrete functional is

$$
j_\tau(u_\tau):=J_\tau(S_\tau u_\tau,u_\tau).
$$

The unconstrained discrete problem is

$$
\min_{u_\tau\in U^N} j_\tau(u_\tau).
$$

For constrained controls, a natural time-discrete admissible set is

$$
U_{\mathrm{ad},\tau}
\subset U^N.
$$

For box constraints this will be specified later level by level.

Since the state equation is linear and the cost functional is quadratic, the
reduced functional $j_\tau$ is quadratic.  Because $\alpha>0$, it is strictly
convex in the control.  Therefore the unconstrained problem has a unique
minimizer, and the constrained problem has a unique minimizer whenever
$U_{\mathrm{ad},\tau}$ is nonempty, closed, and convex.

The remaining task is to compute the gradient of $j_\tau$ without explicitly
forming the derivative of $S_\tau$.

---

## Discrete Lagrangian

We now derive the discrete adjoint equation directly from the time-discrete
optimization problem.  This is the discretize-then-optimize approach.

For each time level define the state residual

$$
R^k(y_\tau,u_\tau)\in V'
$$

by

$$
\langle R^k(y_\tau,u_\tau),v\rangle_{V',V}
:=
(y^k-y^{k-1},v)_H
+\tau a(y^k,v)
-\tau\langle Bu^k+f^k,v\rangle_{V',V}.
$$

The discrete state equation is

$$
R^k(y_\tau,u_\tau)=0,
\qquad k=1,\ldots,N.
$$

Introduce Lagrange multipliers

$$
p^1,\ldots,p^N\in V.
$$

We choose the discrete Lagrangian

$$
\mathscr L_\tau(y_\tau,u_\tau,p_\tau)
:=
J_\tau(y_\tau,u_\tau)
-\sum_{k=1}^N
\langle R^k(y_\tau,u_\tau),p^k\rangle_{V',V}.
$$

That is,

$$
\begin{aligned}
\mathscr L_\tau
=&
\frac{\tau}{2}\sum_{k=1}^N \|y^k-y_d^k\|_H^2
+\frac{\alpha\tau}{2}\sum_{k=1}^N \|u^k\|_U^2
\\
&-
\sum_{k=1}^N
\Big[
(y^k-y^{k-1},p^k)_H
+\tau a(y^k,p^k)
-\tau\langle Bu^k+f^k,p^k\rangle_{V',V}
\Big].
\end{aligned}
$$

This sign convention is chosen so that the discrete adjoint equation has the
same sign as the continuous adjoint equation and the reduced gradient is
$\alpha u^k+B^*p^k$.

---

## Variation with Respect to the State

Let $z^1,\ldots,z^N\in V$ be arbitrary variations of the state.  The initial
value $y^0=y_0$ is fixed, so $z^0=0$.

The derivative of the tracking term is

$$
\tau\sum_{k=1}^N (y^k-y_d^k,z^k)_H.
$$

The derivative of the residual part in the Lagrangian is

$$
-
\sum_{k=1}^N
\Big[
(z^k-z^{k-1},p^k)_H
+\tau a(z^k,p^k)
\Big].
$$

The term containing $z^k$ appears in two neighboring residuals, with the
overall minus sign from the Lagrangian:

- from $R^k$, through $(z^k,p^k)_H+\tau a(z^k,p^k)$;
- from $R^{k+1}$, through $-(z^k,p^{k+1})_H$.

For $k=1,\ldots,N-1$, the coefficient of $z^k$ is therefore

$$
\tau(y^k-y_d^k,z^k)_H
-(z^k,p^k)_H
+(z^k,p^{k+1})_H
-\tau a(z^k,p^k).
$$

For $k=N$, there is no residual $R^{N+1}$.  It is convenient to introduce the
terminal value

$$
p^{N+1}:=\beta(y^N-y_T).
$$

Then all time levels can be written uniformly as

$$
\tau(y^k-y_d^k,z^k)_H
-(z^k,p^k-p^{k+1})_H
-\tau a(z^k,p^k).
$$

The stationarity condition with respect to $y^k$ is therefore:

find $p^1,\ldots,p^N\in V$, with $p^{N+1}=\beta(y^N-y_T)$, such that

$$
(z,p^k-p^{k+1})_H
+\tau a(z,p^k)
=
\tau (y^k-y_d^k,z)_H
\qquad \forall z\in V,
$$

for $k=N,N-1,\ldots,1$.

Dividing by $\tau$ gives

$$
\left(\frac{p^k-p^{k+1}}{\tau},z\right)_H
+a(z,p^k)
=
(y^k-y_d^k,z)_H.
$$

This is an implicit Euler scheme run backward in time for the adjoint
equation

$$
-p_t+A^*p=y-y_d,
\qquad
p(T)=\beta(y(T)-y_T).
$$

In boxed form:

$$
\boxed{
\left(\frac{p^k-p^{k+1}}{\tau},z\right)_H
+a(z,p^k)
=
(y^k-y_d^k,z)_H
\quad \forall z\in V,
\qquad
p^{N+1}=\beta(y^N-y_T).
}
$$

This equation is solved backward:

$$
k=N,N-1,\ldots,1.
$$

---

## Variation with Respect to the Control

We now compute the derivative with respect to $u^k$.

Using the same Lagrangian convention as above, the control derivative at time
level $k$ is

$$
\tau\alpha(u^k,w)_U
+\tau\langle Bw,p^k\rangle_{V',V}
$$

for every control variation $w\in U$.

By definition of $B^*$,

$$
\langle Bw,p^k\rangle_{V',V}
=
(B^*p^k,w)_U.
$$

Hence

$$
\frac{\partial j_\tau}{\partial u^k}(u_\tau)
=
\tau(\alpha u^k+B^*p^k)
$$

as an element of $U$ when the product space $U^N$ is equipped with the
unweighted product inner product.

It is often more natural to equip the discrete control space with the
time-weighted inner product

$$
(u_\tau,w_\tau)_{U_\tau}
:=
\tau\sum_{k=1}^N (u^k,w^k)_U.
$$

With respect to this inner product, the gradient is

$$
\nabla j_\tau(u_\tau)^k
=
\alpha u^k+B^*p^k,
\qquad k=1,\ldots,N.
$$

This distinction is important in implementations:

- the algebraic derivative contains the quadrature weight $\tau$;
- the $L^2(0,T;U)$-consistent gradient does not.

---

## Discrete Optimality System

The unconstrained discrete optimality system is:

find

$$
(y^1,\ldots,y^N),
\qquad
(p^1,\ldots,p^N),
\qquad
(u^1,\ldots,u^N),
$$

with $y^0=y_0$ and $p^{N+1}=\beta(y^N-y_T)$, such that for all $v,z\in V$ and
$w\in U$:

$$
\boxed{
\left(\frac{y^k-y^{k-1}}{\tau},v\right)_H
+a(y^k,v)
=
\langle Bu^k+f^k,v\rangle_{V',V}
}
$$

for $k=1,\ldots,N$,

$$
\boxed{
\left(\frac{p^k-p^{k+1}}{\tau},z\right)_H
+a(z,p^k)
=
(y^k-y_d^k,z)_H
}
$$

for $k=N,\ldots,1$, and

$$
\boxed{
\alpha u^k+B^*p^k=0
}
$$

for $k=1,\ldots,N$.

The structure is forward-backward:

- the state equation propagates $y^0\mapsto y^1\mapsto\cdots\mapsto y^N$;
- the adjoint equation propagates $p^{N+1}\mapsto p^N\mapsto\cdots\mapsto p^1$;
- the control equation couples state and adjoint at the same time level.

This is the time-discrete analogue of the continuous parabolic optimality
system.

---

## Discretize-Then-Optimize and Optimize-Then-Discretize

There are two common routes to a numerical optimality system.

**Discretize then optimize (DTO).**
First discretize the state equation and the objective functional.  This gives
a finite-dimensional or finite-time-dimensional optimization problem.  Then
derive its first-order optimality conditions.

This is what we did above.

**Optimize then discretize (OTD).**
First derive the continuous optimality system:

$$
\begin{cases}
y_t+Ay=Bu+f,\\
-p_t+A^*p=y-y_d,\\
\alpha u+B^*p=0.
\end{cases}
$$

Then discretize this coupled forward-backward system.

For implicit Euler, the DTO adjoint equation is

$$
\left(\frac{p^k-p^{k+1}}{\tau},z\right)_H
+a(z,p^k)
=
(y^k-y_d^k,z)_H.
$$

This is a backward implicit Euler discretization of the adjoint equation.
Thus, for this simple linear-quadratic problem, DTO and OTD can be made to
match if the adjoint equation is discretized with the time-reversed scheme
that is algebraically adjoint to the state scheme.

The important warning is that this agreement is not automatic.
It depends on:

- the quadrature rule used in the objective;
- the time locations of $u^k$, $y^k$, and $p^k$;
- the discrete inner products used to identify gradients;
- the treatment of terminal terms;
- the choice of time-stepping method.

DTO guarantees that the computed adjoint is the exact adjoint of the discrete
state equation.  OTD gives the same result only if the discretization is chosen
in an adjoint-consistent way.

For gradient-based optimization of a discretized problem, DTO is usually the
safer point of view: the gradient being used is genuinely the gradient of the
discrete objective.

---

## Fully Discrete Finite Element Form

We now add a standard finite element discretization in space.

Let

$$
V_h\subset V
$$

be a finite element space with basis $\{\varphi_i\}_{i=1}^{n_y}$.
Let the control be represented in a finite-dimensional space $U_h$ with basis
$\{\psi_j\}_{j=1}^{n_u}$.

Write the coefficient vectors as

$$
Y^k\in\mathbb R^{n_y},
\qquad
P^k\in\mathbb R^{n_y},
\qquad
U^k\in\mathbb R^{n_u}.
$$

Define the mass and stiffness matrices

$$
M_{ij}:=(\varphi_j,\varphi_i)_H,
\qquad
K_{ij}:=a(\varphi_j,\varphi_i).
$$

Define also the discrete control operator $B_h$ by

$$
(B_h)_{ij}:=\langle B\psi_j,\varphi_i\rangle_{V',V}.
$$

The state equation becomes

$$
\boxed{
(M+\tau K)Y^k
=
MY^{k-1}
+\tau B_h U^k
+\tau F^k.
}
$$

The adjoint equation becomes

$$
\boxed{
(M+\tau K)^T P^k
=
M^T P^{k+1}
+\tau M(Y^k-Y_d^k).
}
$$

If $M$ and $K$ are symmetric, this simplifies to

$$
(M+\tau K)P^k
=
MP^{k+1}
+\tau M(Y^k-Y_d^k).
$$

Let $M_u$ be the control mass matrix.  The discrete control equation is

$$
\boxed{
\alpha M_u U^k + B_h^T P^k = 0.
}
$$

Equivalently, after applying $M_u^{-1}$,

$$
\alpha U^k + M_u^{-1}B_h^T P^k=0.
$$

This last form is the coefficient representation of the $L^2$-gradient in the
control space.  As in the elliptic case, one must distinguish:

- the algebraic residual $\alpha M_u U^k+B_h^T P^k$;
- the Riesz-represented gradient $\alpha U^k+M_u^{-1}B_h^T P^k$.

This distinction becomes important for projected-gradient and active-set
methods.

---

## All-at-Once Structure

The time-stepping form suggests solving the state forward and the adjoint
backward.  However, the whole system can also be assembled as one global
space-time KKT system.

Define

$$
A_\tau:=M+\tau K.
$$

The state equations are

$$
A_\tau Y^1-\tau B_hU^1 = MY^0+\tau F^1,
$$

and, for $k=2,\ldots,N$,

$$
-MY^{k-1}+A_\tau Y^k-\tau B_hU^k=\tau F^k.
$$

Thus the global state operator is block lower bidiagonal:

$$
E
=
\begin{pmatrix}
A_\tau \\
-M & A_\tau \\
& -M & A_\tau \\
& & \ddots & \ddots \\
& & & -M & A_\tau
\end{pmatrix}.
$$

The adjoint operator is the transpose block structure:

$$
E^T
=
\begin{pmatrix}
A_\tau^T & -M^T \\
& A_\tau^T & -M^T \\
& & \ddots & \ddots \\
& & & A_\tau^T & -M^T \\
& & & & A_\tau^T
\end{pmatrix}.
$$

This is the algebraic reason why the adjoint runs backward in time.

The full unconstrained KKT system has the schematic form

$$
\begin{pmatrix}
\tau \mathcal M_y & E^T & 0\\
E & 0 & -\tau\mathcal B\\
0 & -\tau\mathcal B^T & \alpha\tau\mathcal M_u
\end{pmatrix}
\begin{pmatrix}
Y\\
P\\
U
\end{pmatrix}
=
\begin{pmatrix}
\tau\mathcal M_yY_d\\
F_{\tau,y_0}\\
0
\end{pmatrix}.
$$

Here $\mathcal M_y$, $\mathcal M_u$, and $\mathcal B$ are block-diagonal
space-time matrices containing $M$, $M_u$, and $B_h$ at each time level.

This all-at-once formulation is useful because it exposes the saddle-point
structure of the full problem.  It is also the natural form for space-time
preconditioners.  In this lecture, however, the main point is more basic:

> the backward adjoint equation is the transpose of the forward time-stepping
> operator.

---

## Box Constraints

Let the continuous admissible set be a pointwise box in time:

$$
U_{\mathrm{ad}}
=
\{u\in L^2(0,T;U): u_a(t)\le u(t)\le u_b(t)
\text{ for a.e. }t\in(0,T)\}.
$$

After time discretization we obtain bounds

$$
u_a^k\le u^k\le u_b^k,
\qquad k=1,\ldots,N,
$$

and the discrete admissible set

$$
U_{\mathrm{ad},\tau}
=
\{(u^1,\ldots,u^N)\in U^N:
u_a^k\le u^k\le u_b^k,\ k=1,\ldots,N\}.
$$

The discrete first-order condition is the variational inequality

$$
\tau\sum_{k=1}^N
(\alpha \bar u^k+B^*\bar p^k,
u^k-\bar u^k)_U
\ge 0
\qquad
\forall u_\tau\in U_{\mathrm{ad},\tau}.
$$

Since $\tau>0$, this is equivalent to

$$
\sum_{k=1}^N
(\alpha \bar u^k+B^*\bar p^k,
u^k-\bar u^k)_U
\ge 0
\qquad
\forall u_\tau\in U_{\mathrm{ad},\tau}.
$$

For a box constraint, the projection formula is level-wise:

$$
\boxed{
\bar u^k
=
P_{[u_a^k,u_b^k]}
\left(
-\frac1\alpha B^*\bar p^k
\right),
\qquad
k=1,\ldots,N.
}
$$

In coefficient form, if $U^k$ denotes the vector of control coefficients, the
Riesz-represented gradient is

$$
G_u^k
=
\alpha U^k+M_u^{-1}B_h^T P^k.
$$

This is the quantity that should be compared with primal bound gaps in
projected-gradient or active-set methods.  The raw algebraic residual

$$
\alpha M_uU^k+B_h^TP^k
$$

lives in the dual coefficient space and has different scaling.

Thus the active-set logic from the elliptic case applies independently at
each time level, but the state and adjoint equations still couple all time
levels through forward and backward propagation.

---

## Reduced Numerical Strategy

The reduced viewpoint eliminates the state and adjoint from the optimization
variables.

Given a control sequence $u_\tau$:

1. solve the state equation forward for $y^1,\ldots,y^N$;
2. solve the adjoint equation backward for $p^N,\ldots,p^1$;
3. assemble the gradient
   $$
   g^k=\alpha u^k+B^*p^k;
   $$
4. update the control with a gradient, conjugate-gradient, quasi-Newton, or
   projected method.

The cost of one gradient evaluation is essentially:

- one forward parabolic solve;
- one backward parabolic solve;
- one application of $B^*$ at each time level.

This is why the adjoint method is indispensable: the cost of the gradient does
not scale with the number of control degrees of freedom.

For box-constrained problems the same loop is used, but the update step is
projected or active-set based.

---

## What Should Be Remembered

Time discretization is not just a technical detail.  It determines the exact
algebraic optimization problem being solved.

For implicit Euler:

- the state equation is a forward recurrence;
- the discrete adjoint is a backward recurrence;
- the adjoint recurrence is the transpose of the state recurrence;
- the discrete gradient at time level $k$ is
  $$
  \alpha u^k+B^*p^k;
  $$
- in finite elements the control gradient requires the control mass matrix;
- box constraints give a projection formula at each time level.

The key DTO lesson is:

> derive the adjoint from the discrete equations if you want the exact
> gradient of the discrete objective.

The OTD viewpoint is still valuable, because it explains what the discrete
system approximates.  But for implementation and optimization, DTO is the
most reliable way to avoid sign errors, index shifts, and inconsistent
gradients.

---

## References

- F. Troeltzsch, *Optimal Control of Partial Differential Equations*, AMS, 2010.
- J. C. De los Reyes, *Numerical PDE-Constrained Optimization*, Springer, 2015.
- M. Hinze, R. Pinnau, M. Ulbrich, and S. Ulbrich,
  *Optimization with PDE Constraints*, Springer, 2009.
- A. Manzoni, A. Quarteroni, S. Salsa,
  *Optimal Control of Partial Differential Equations*, Springer, 2021.

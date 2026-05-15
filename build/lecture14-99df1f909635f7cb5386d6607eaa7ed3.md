# Boundary Control and Parameter Estimation

## Overview

The previous lectures focused mostly on distributed controls, where the
control acts as a volume source term.  This lecture moves to two situations
where the same optimal-control structure remains visible, but the functional
analytic details become more delicate:

- controls acting on the boundary;
- unknown coefficients or material parameters.

The main new point is that the control-to-state map no longer has the same
regularity as in the distributed case.  Boundary controls interact with trace
operators, while parameter controls enter the differential operator itself.
Both effects change the correct space in which gradients and optimality
conditions should be interpreted.

The logical path is:

- formulate Neumann boundary control in weak form;
- prove well-posedness by the trace theorem and Lax-Milgram;
- derive the adjoint equation and the reduced gradient using a Lagrangian;
- discuss Dirichlet boundary control through liftings;
- formulate coefficient identification as a PDE-constrained inverse problem;
- prove the sensitivity equation for the coefficient-to-state map;
- derive the adjoint gradient for parameter estimation from a Lagrangian;
- interpret box constraints as variational inequalities and projection
  formulas when enough regularity is available.

Throughout the lecture, let $\Omega\subset\mathbb R^d$ be a bounded Lipschitz
domain, and let $\Gamma:=\partial\Omega$.

---

## Boundary Controls and Trace Spaces

For a distributed control, the state equation typically has the form

$$
\mathcal A y = Bu+f
\qquad\text{in }\Omega,
$$

where $B$ maps the control space into a volume dual space such as
$H^{-1}(\Omega)$.

For a boundary control, the action of $u$ is instead mediated by the trace
operator.  The fundamental mapping is

$$
\gamma:H^1(\Omega)\to H^{1/2}(\Gamma),
$$

which is continuous.  Since $H^{1/2}(\Gamma)$ embeds continuously into
$L^2(\Gamma)$ on bounded Lipschitz boundaries, there is a constant
$C_{\mathrm{tr}}>0$ such that

$$
\|\gamma v\|_{L^2(\Gamma)}
\le C_{\mathrm{tr}}\|v\|_{H^1(\Omega)}
\qquad \forall v\in H^1(\Omega).
$$

This estimate is the key reason why Neumann data in $L^2(\Gamma)$ define a
bounded linear functional on $H^1(\Omega)$:

$$
\left|\int_\Gamma u\,\gamma v\,ds\right|
\le
\|u\|_{L^2(\Gamma)}\|\gamma v\|_{L^2(\Gamma)}
\le
C_{\mathrm{tr}}\|u\|_{L^2(\Gamma)}\|v\|_{H^1(\Omega)}.
$$

Thus Neumann boundary control is naturally compatible with weak formulations.
Dirichlet boundary control is subtler, because prescribing $y=u$ on $\Gamma$
requires $u$ to be a trace of an $H^1(\Omega)$ function, hence

$$
u\in H^{1/2}(\Gamma)
$$

is the natural energy space.

---

## A Coercive Neumann Boundary Control Model

To avoid the compatibility condition of the pure Neumann Laplacian, we first
consider the coercive elliptic equation

$$
\begin{cases}
-\Delta y+\sigma y=f & \text{in }\Omega,\\
\partial_n y=u & \text{on }\Gamma,
\end{cases}
$$

where $\sigma>0$.  The pure Neumann case $\sigma=0$ is discussed below.

Let

$$
V:=H^1(\Omega),
\qquad
U:=L^2(\Gamma),
$$

and define

$$
a(y,v):=
\int_\Omega \nabla y\cdot\nabla v\,dx
+
\sigma\int_\Omega yv\,dx,
$$

$$
F(v):=\int_\Omega f v\,dx,
\qquad
B(u,v):=\int_\Gamma u\,\gamma v\,ds.
$$

The weak state equation is:

find $y\in V$ such that

$$
a(y,v)=F(v)+B(u,v)
\qquad \forall v\in V.
$$

For a desired state $y_d\in L^2(\Omega)$ and a regularization parameter
$\alpha>0$, the boundary control problem is

$$
\min_{u\in U_{\mathrm{ad}}}
J(y,u)
:=
\frac12\|y-y_d\|_{L^2(\Omega)}^2
+
\frac\alpha2\|u\|_{L^2(\Gamma)}^2
$$

subject to the weak state equation.

The admissible set is assumed to be nonempty, closed, and convex in
$L^2(\Gamma)$.  A standard example is the box-constrained set

$$
U_{\mathrm{ad}}
:=
\{u\in L^2(\Gamma): u_a\le u\le u_b \text{ a.e. on }\Gamma\},
$$

where $u_a,u_b\in L^\infty(\Gamma)$ and $u_a\le u_b$ a.e.

---

## Well-Posedness of the Neumann State Equation

The bilinear form $a$ is continuous on $V\times V$ because

$$
|a(y,v)|
\le
\|\nabla y\|_{L^2(\Omega)}\|\nabla v\|_{L^2(\Omega)}
+
\sigma\|y\|_{L^2(\Omega)}\|v\|_{L^2(\Omega)}
\le
C\|y\|_{H^1(\Omega)}\|v\|_{H^1(\Omega)}.
$$

It is also coercive on $H^1(\Omega)$:

$$
a(v,v)
=
\|\nabla v\|_{L^2(\Omega)}^2
+
\sigma\|v\|_{L^2(\Omega)}^2
\ge
\min\{1,\sigma\}\|v\|_{H^1(\Omega)}^2.
$$

The right-hand side is bounded on $V$.  Indeed,

$$
|F(v)|
\le
\|f\|_{L^2(\Omega)}\|v\|_{L^2(\Omega)}
\le
\|f\|_{L^2(\Omega)}\|v\|_{H^1(\Omega)},
$$

and the trace estimate gives

$$
|B(u,v)|
\le
C_{\mathrm{tr}}\|u\|_{L^2(\Gamma)}\|v\|_{H^1(\Omega)}.
$$

By Lax-Milgram, for every $u\in L^2(\Gamma)$ there exists a unique state
$y(u)\in H^1(\Omega)$ and

$$
\|y(u)\|_{H^1(\Omega)}
\le
C\left(
\|f\|_{L^2(\Omega)}
+
\|u\|_{L^2(\Gamma)}
\right).
$$

Moreover, the control-to-state map

$$
S:L^2(\Gamma)\to H^1(\Omega),
\qquad
S(u)=y(u),
$$

is affine and continuous.  If $f$ is fixed, its derivative is the linear map
$S'(u)h=z$, where $z\in H^1(\Omega)$ solves

$$
a(z,v)=B(h,v)
\qquad \forall v\in H^1(\Omega).
$$

For the pure Neumann equation $\sigma=0$, constants belong to the kernel of
the operator.  One either works in the quotient space $H^1(\Omega)/\mathbb R$,
or fixes the representative by imposing zero mean:

$$
V_0:=
\left\{
v\in H^1(\Omega):
\int_\Omega v\,dx=0
\right\}.
$$

On $V_0$, Poincare's inequality gives coercivity of
$a(v,v)=\|\nabla v\|_{L^2(\Omega)}^2$.  If the equation is required to hold
for all test functions in $H^1(\Omega)$, the data must also satisfy the
compatibility condition

$$
\int_\Omega f\,dx+\int_\Gamma u\,ds=0.
$$

This is why adding a positive reaction term is often the cleaner model for
unconstrained Neumann boundary control.

---

## Reduced Functional and Existence of an Optimal Control

Since the state is uniquely determined by $u$, define the reduced functional

$$
j(u):=J(S(u),u).
$$

The estimate above implies that $j$ is continuous on $L^2(\Gamma)$.  The term
$\frac\alpha2\|u\|_{L^2(\Gamma)}^2$ makes $j$ coercive:

$$
j(u)\to+\infty
\qquad\text{as}\qquad
\|u\|_{L^2(\Gamma)}\to+\infty.
$$

Because $S$ is affine and continuous, the tracking term
$\frac12\|S(u)-y_d\|_{L^2(\Omega)}^2$ is weakly lower semicontinuous.  The
regularization term is weakly lower semicontinuous as well.  Therefore, by
the direct method of the calculus of variations, $j$ attains a minimizer on
every nonempty closed convex set $U_{\mathrm{ad}}\subset L^2(\Gamma)$.

Since $\alpha>0$, the reduced functional is strictly convex.  Hence the
optimal control is unique.

---

## Adjoint Equation for Neumann Boundary Control

We derive the adjoint through the Lagrangian formalism.  The weak state
equation is the constraint

$$
a(y,v)-F(v)-B(u,v)=0
\qquad \forall v\in V.
$$

Introduce the Lagrangian

$$
\mathcal L(y,u,p)
=
\frac12\|y-y_d\|_{L^2(\Omega)}^2
+
\frac\alpha2\|u\|_{L^2(\Gamma)}^2
-
\bigl[
a(y,p)-F(p)-B(u,p)
\bigr],
$$

where $p\in V$ is the Lagrange multiplier.  The sign in front of the
constraint is a convention; this choice gives the gradient
$\alpha u+\gamma p$.

Stationarity with respect to $p$ gives the state equation:

$$
D_p\mathcal L(y,u,p)\,\psi
=
-\bigl[
a(y,\psi)-F(\psi)-B(u,\psi)
\bigr]
=0
\qquad \forall \psi\in V.
$$

Stationarity with respect to $y$ gives the adjoint equation.  For every
state variation $\eta\in V$,

$$
D_y\mathcal L(y,u,p)\,\eta
=
(y-y_d,\eta)_{L^2(\Omega)}
-
a(\eta,p).
$$

Thus

$$
a(\eta,p)=(y-y_d,\eta)_{L^2(\Omega)}
\qquad \forall \eta\in V.
$$

The adjoint problem is well posed by the same Lax-Milgram argument.  If the
solution is smooth, this weak equation corresponds to

$$
\begin{cases}
-\Delta p+\sigma p=y-y_d & \text{in }\Omega,\\
\partial_n p=0 & \text{on }\Gamma.
\end{cases}
$$

Finally, stationarity with respect to the control gives the derivative in
the control direction.  For $h\in L^2(\Gamma)$,

$$
D_u\mathcal L(y,u,p)\,h
=
\alpha(u,h)_{L^2(\Gamma)}
+
B(h,p)
=
\int_\Gamma
\left(
\alpha u+\gamma p
\right)h\,ds.
$$

Because the state and adjoint stationarity equations cancel all implicit
state variations, this is the derivative of the reduced functional:

$$
j'(u)h
=
\int_\Gamma
\left(
\alpha u+\gamma p
\right)h\,ds.
$$

Thus the $L^2(\Gamma)$-gradient of the reduced functional is

$$
\nabla j(u)=\alpha u+\gamma p.
$$

If one puts the opposite sign in front of the constraint in the Lagrangian,
the adjoint variable changes sign and the same optimality system can be
written with $\alpha u-\gamma p$ instead.

---

## Variational Inequality and Projection Formula

Let $\bar u\in U_{\mathrm{ad}}$ be the optimal control, with state $\bar y$
and adjoint $\bar p$.  Since $U_{\mathrm{ad}}$ is closed and convex, the
first-order necessary and sufficient condition is

$$
j'(\bar u)(v-\bar u)\ge0
\qquad \forall v\in U_{\mathrm{ad}}.
$$

Using the gradient formula, this becomes

$$
(\alpha\bar u+\gamma\bar p,v-\bar u)_{L^2(\Gamma)}
\ge0
\qquad \forall v\in U_{\mathrm{ad}}.
$$

For box constraints, the variational inequality is equivalent to the
pointwise projection formula

$$
\bar u
=
P_{[u_a,u_b]}
\left(
-\frac1\alpha\gamma\bar p
\right)
\qquad\text{a.e. on }\Gamma.
$$

Indeed, for almost every boundary point $x\in\Gamma$, the scalar condition is

$$
(\alpha\bar u(x)+\gamma\bar p(x))(v-\bar u(x))\ge0
\qquad
\forall v\in [u_a(x),u_b(x)].
$$

This says exactly that

$$
-\bigl(\alpha\bar u(x)+\gamma\bar p(x)\bigr)
\in
N_{[u_a(x),u_b(x)]}(\bar u(x)),
$$

where $N_C$ denotes the normal cone to the convex set $C$.  The normal-cone
identity for projections gives the formula above.

Equivalently:

- if $u_a<\bar u<u_b$, then $\alpha\bar u+\gamma\bar p=0$;
- if $\bar u=u_a$, then $\alpha\bar u+\gamma\bar p\ge0$;
- if $\bar u=u_b$, then $\alpha\bar u+\gamma\bar p\le0$.

---

## Dirichlet Boundary Control

For Dirichlet boundary control, the state equation is

$$
\begin{cases}
-\Delta y+\sigma y=f & \text{in }\Omega,\\
y=u & \text{on }\Gamma.
\end{cases}
$$

Now the control cannot be an arbitrary element of $L^2(\Gamma)$ if the state
is required to lie in $H^1(\Omega)$.  The natural space is

$$
U=H^{1/2}(\Gamma).
$$

Let $V_0:=H_0^1(\Omega)$.  By the trace theorem, there exists a continuous
lifting operator

$$
R:H^{1/2}(\Gamma)\to H^1(\Omega)
$$

such that

$$
\gamma(Ru)=u.
$$

Write

$$
y=w+Ru,
\qquad
w\in H_0^1(\Omega).
$$

The weak equation for $w$ is

$$
a(w,v)=F(v)-a(Ru,v)
\qquad \forall v\in H_0^1(\Omega).
$$

Since $a$ is coercive on $H_0^1(\Omega)$, Lax-Milgram gives a unique
$w\in H_0^1(\Omega)$, hence a unique state $y=w+Ru$.

The estimate is

$$
\|y\|_{H^1(\Omega)}
\le
C\left(
\|f\|_{H^{-1}(\Omega)}
+
\|u\|_{H^{1/2}(\Gamma)}
\right).
$$

Thus the Dirichlet control-to-state map is continuous from
$H^{1/2}(\Gamma)$ into $H^1(\Omega)$.

---

## Dirichlet Gradient and Regularity Warning

For Dirichlet control, it is useful to use the lifted variable
$y=w+Ru$, with $w\in H_0^1(\Omega)$.  The weak constraint is

$$
a(w,v)-F(v)+a(Ru,v)=0
\qquad \forall v\in H_0^1(\Omega).
$$

The Lagrangian is

$$
\mathcal L(w,u,p)
=
\frac12\|w+Ru-y_d\|_{L^2(\Omega)}^2
+
\frac\alpha2\|u\|_U^2
-
\bigl[
a(w,p)-F(p)+a(Ru,p)
\bigr],
$$

where $p\in H_0^1(\Omega)$ is the Lagrange multiplier.  Stationarity with
respect to $p$ gives the lifted state equation.  Stationarity with respect to
$w$ gives, for every $\eta\in H_0^1(\Omega)$,

$$
D_w\mathcal L(w,u,p)\,\eta
=
(y-y_d,\eta)_{L^2(\Omega)}
-
a(\eta,p)
=0.
$$

Therefore the adjoint equation is

$$
a(\eta,p)=(y-y_d,\eta)_{L^2(\Omega)}
\qquad \forall \eta\in H_0^1(\Omega).
$$

The control derivative is obtained by differentiating the same Lagrangian
with respect to $u$.  For $h\in U=H^{1/2}(\Gamma)$,

$$
D_u\mathcal L(w,u,p)\,h
=
\alpha\langle R_U u,h\rangle_{U',U}
+
(y-y_d,Rh)_{L^2(\Omega)}
-
a(Rh,p).
$$

The subtle point is that $Rh$ is generally not an element of
$H_0^1(\Omega)$, because its trace is $h$.  Therefore one cannot simply insert
$Rh$ as a test function in the adjoint equation, which is valid only for
zero-trace variations.

Assume for a moment that $p$ is smooth enough for Green's identity to be used
classically.  Since the adjoint satisfies

$$
-\Delta p+\sigma p=y-y_d
\qquad\text{in }\Omega,
\qquad
p=0
\qquad\text{on }\Gamma,
$$

we compute

$$
a(Rh,p)
=
\int_\Omega \nabla(Rh)\cdot\nabla p\,dx
+
\sigma\int_\Omega Rh\,p\,dx.
$$

Integrating the first term by parts gives

$$
a(Rh,p)
=
\int_\Omega Rh(-\Delta p+\sigma p)\,dx
+
\int_\Gamma \gamma(Rh)\,\partial_n p\,ds.
$$

Using $\gamma(Rh)=h$ and the strong adjoint equation, this becomes

$$
a(Rh,p)
=
(y-y_d,Rh)_{L^2(\Omega)}
+
\int_\Gamma h\,\partial_n p\,ds.
$$

Hence the two non-regularization terms in the Lagrangian control derivative
reduce to the boundary term

$$
(y-y_d,Rh)_{L^2(\Omega)}-a(Rh,p)
=
-\langle \partial_n p,h\rangle_{H^{-1/2},H^{1/2}}.
$$

Hence the control stationarity condition is naturally written as

$$
\alpha R_U \bar u-\partial_n \bar p=0
\qquad\text{in }H^{-1/2}(\Gamma),
$$

where $R_U:U\to U'$ is the Riesz map of $H^{1/2}(\Gamma)$.

This is the key distinction from Neumann control:

- Neumann control with $u\in L^2(\Gamma)$ gives an $L^2(\Gamma)$ gradient
  involving the trace $\gamma p$;
- Dirichlet control gives a boundary gradient involving the normal derivative
  $\partial_n p$, which is generally only an element of $H^{-1/2}(\Gamma)$.

An $L^2(\Gamma)$ projection formula for Dirichlet control is therefore a
formal or additional-regularity statement.  It is valid, for instance, if
$\partial_n p\in L^2(\Gamma)$ and the control is regularized in
$L^2(\Gamma)$:

$$
\bar u
=
P_{[u_a,u_b]}
\left(
\frac1\alpha\partial_n\bar p
\right).
$$

Without such regularity, the correct condition is the variational inequality
in the duality pairing between $H^{-1/2}(\Gamma)$ and $H^{1/2}(\Gamma)$.

---

## Parameter Estimation Problem

We now consider an inverse problem in which the unknown is a coefficient in
the PDE.  Let

$$
0<q_{\min}\le q(x)\le q_{\max}<\infty
\qquad\text{a.e. in }\Omega.
$$

For each admissible coefficient $q$, the state $y=y(q)\in H_0^1(\Omega)$ is
defined by

$$
\int_\Omega q\nabla y\cdot\nabla v\,dx
=
\langle f,v\rangle_{H^{-1},H_0^1}
\qquad \forall v\in H_0^1(\Omega).
$$

In strong form this corresponds to

$$
\begin{cases}
-\nabla\cdot(q\nabla y)=f & \text{in }\Omega,\\
y=0 & \text{on }\Gamma.
\end{cases}
$$

The coefficient identification problem is to choose $q$ so that $y(q)$ fits
observed data $y_d$.  A standard Tikhonov functional is

$$
j(q)
=
\frac12\|y(q)-y_d\|_{L^2(\Omega)}^2
+
\frac\alpha2\|q-q_{\mathrm{ref}}\|_{L^2(\Omega)}^2.
$$

The admissible set is often taken as

$$
Q_{\mathrm{ad}}
=
\{q\in L^\infty(\Omega):
q_a\le q\le q_b \text{ a.e. in }\Omega\},
$$

with $0<q_a\le q_b$ a.e.

For every $q\in Q_{\mathrm{ad}}$, the bilinear form

$$
a_q(y,v):=\int_\Omega q\nabla y\cdot\nabla v\,dx
$$

is continuous and coercive on $H_0^1(\Omega)$:

$$
|a_q(y,v)|
\le
\|q\|_{L^\infty(\Omega)}
\|\nabla y\|_{L^2(\Omega)}
\|\nabla v\|_{L^2(\Omega)},
$$

$$
a_q(v,v)
\ge
q_{\min}\|\nabla v\|_{L^2(\Omega)}^2.
$$

Lax-Milgram therefore gives a unique state $y(q)$ and the stability estimate

$$
\|y(q)\|_{H_0^1(\Omega)}
\le
\frac{C_P}{q_{\min}}\|f\|_{H^{-1}(\Omega)}.
$$

---

## Sensitivity Equation for the Coefficient

The coefficient-to-state map is nonlinear.  We derive its derivative.

Let $\delta q\in L^\infty(\Omega)$ be a perturbation such that
$q+t\delta q$ remains uniformly positive for small $t$.  Let

$$
y_t:=y(q+t\delta q),
\qquad
y:=y(q).
$$

Subtracting the two weak state equations gives

$$
\int_\Omega q\nabla(y_t-y)\cdot\nabla v\,dx
=
-t\int_\Omega \delta q\nabla y_t\cdot\nabla v\,dx
\qquad \forall v\in H_0^1(\Omega).
$$

Divide by $t$ and pass formally to the limit $t\to0$.  The sensitivity

$$
z:=y'(q)\delta q
$$

solves

$$
\int_\Omega q\nabla z\cdot\nabla v\,dx
=
-\int_\Omega \delta q\nabla y\cdot\nabla v\,dx
\qquad \forall v\in H_0^1(\Omega).
$$

In strong form,

$$
\begin{cases}
-\nabla\cdot(q\nabla z)=\nabla\cdot(\delta q\nabla y) & \text{in }\Omega,\\
z=0 & \text{on }\Gamma.
\end{cases}
$$

The equation is well posed because the right-hand side is bounded on
$H_0^1(\Omega)$:

$$
\left|
\int_\Omega \delta q\nabla y\cdot\nabla v\,dx
\right|
\le
\|\delta q\|_{L^\infty(\Omega)}
\|\nabla y\|_{L^2(\Omega)}
\|\nabla v\|_{L^2(\Omega)}.
$$

Thus

$$
\|z\|_{H_0^1(\Omega)}
\le
C
\|\delta q\|_{L^\infty(\Omega)}
\|y\|_{H_0^1(\Omega)}.
$$

This proves that the derivative maps bounded coefficient perturbations to
state perturbations continuously.

---

## Adjoint Gradient for Parameter Estimation

We again derive the adjoint and the gradient through the Lagrangian.  The weak
state constraint is

$$
\int_\Omega q\nabla y\cdot\nabla v\,dx
-
\langle f,v\rangle_{H^{-1},H_0^1}
=0
\qquad \forall v\in H_0^1(\Omega).
$$

Introduce

$$
\mathcal L(y,q,p)
=
\frac12\|y-y_d\|_{L^2(\Omega)}^2
+
\frac\alpha2\|q-q_{\mathrm{ref}}\|_{L^2(\Omega)}^2
-
\left[
\int_\Omega q\nabla y\cdot\nabla p\,dx
-
\langle f,p\rangle_{H^{-1},H_0^1}
\right].
$$

Stationarity with respect to $p$ gives the state equation.  Stationarity with
respect to $y$ gives, for every $\eta\in H_0^1(\Omega)$,

$$
D_y\mathcal L(y,q,p)\,\eta
=
(y-y_d,\eta)_{L^2(\Omega)}
-
\int_\Omega q\nabla \eta\cdot\nabla p\,dx
=0.
$$

Therefore

$$
\int_\Omega q\nabla \eta\cdot\nabla p\,dx
=
(y-y_d,\eta)_{L^2(\Omega)}
\qquad \forall \eta\in H_0^1(\Omega).
$$

Equivalently, in strong form,

$$
\begin{cases}
-\nabla\cdot(q\nabla p)=y-y_d & \text{in }\Omega,\\
p=0 & \text{on }\Gamma.
\end{cases}
$$

The coefficient derivative is

$$
D_q\mathcal L(y,q,p)\,\delta q
=
\int_\Omega
\left[
\alpha(q-q_{\mathrm{ref}})
-
\nabla y\cdot\nabla p
\right]
\delta q\,dx.
$$

Since the state and adjoint equations are precisely the stationarity
conditions in $p$ and $y$, this is the derivative of the reduced functional:

$$
j'(q)\delta q
=
D_q\mathcal L(y,q,p)\,\delta q.
$$

If the product $\nabla y\cdot\nabla p$ belongs to $L^2(\Omega)$, then the
$L^2(\Omega)$-gradient is

$$
\nabla j(q)
=
\alpha(q-q_{\mathrm{ref}})
-
\nabla y\cdot\nabla p.
$$

In the minimal energy setting, $\nabla y,\nabla p\in L^2(\Omega)^d$, so the
product is a priori only in $L^1(\Omega)$.  Therefore, the formula is always
valid as a duality expression, while interpreting it as an $L^2$ gradient
requires additional regularity or a discrete finite element setting.

---

## Box Constraints for Parameter Estimation

Let $\bar q\in Q_{\mathrm{ad}}$ be a local minimizer, with state $\bar y$ and
adjoint $\bar p$.  The first-order condition is the variational inequality

$$
j'(\bar q)(r-\bar q)\ge0
\qquad \forall r\in Q_{\mathrm{ad}}.
$$

Using the adjoint formula, this reads

$$
\int_\Omega
\left[
\alpha(\bar q-q_{\mathrm{ref}})
-
\nabla\bar y\cdot\nabla\bar p
\right]
(r-\bar q)\,dx
\ge0
\qquad \forall r\in Q_{\mathrm{ad}}.
$$

When the gradient is an $L^2(\Omega)$ function, this is equivalent to the
pointwise projection formula

$$
\bar q
=
P_{[q_a,q_b]}
\left(
q_{\mathrm{ref}}
+
\frac1\alpha\nabla\bar y\cdot\nabla\bar p
\right)
\qquad\text{a.e. in }\Omega.
$$

The corresponding active-set interpretation is:

- on the inactive set $q_a<\bar q<q_b$,

  $$
  \alpha(\bar q-q_{\mathrm{ref}})
  -
  \nabla\bar y\cdot\nabla\bar p=0;
  $$

- on the lower active set $\bar q=q_a$,

  $$
  \alpha(\bar q-q_{\mathrm{ref}})
  -
  \nabla\bar y\cdot\nabla\bar p\ge0;
  $$

- on the upper active set $\bar q=q_b$,

  $$
  \alpha(\bar q-q_{\mathrm{ref}})
  -
  \nabla\bar y\cdot\nabla\bar p\le0.
  $$

---

## Existence and Regularization in Inverse Problems

For distributed and Neumann controls, the regularization term
$\frac\alpha2\|u\|_{L^2}^2$ gives coercivity in the same Hilbert space in
which the control is optimized.

For coefficient estimation, the situation is more subtle:

- the PDE requires $q\in L^\infty(\Omega)$ with a positive lower bound;
- the Tikhonov term above controls only the $L^2(\Omega)$ distance from
  $q_{\mathrm{ref}}$;
- weak convergence of coefficients does not automatically imply strong
  convergence of the corresponding fluxes $q\nabla y$.

In a finite element discretization, existence is automatic because the
admissible set is finite-dimensional and closed.  In an infinite-dimensional
analysis, one usually strengthens the compactness by choosing, for example,

$$
\frac\alpha2\|q-q_{\mathrm{ref}}\|_{H^1(\Omega)}^2
$$

or a bounded-variation regularization.  These choices improve compactness and
make the inverse problem stable with respect to noisy data.

This is the essential role of regularization: it is not merely a numerical
penalty, but part of the mathematical formulation that turns an ill-posed
inverse problem into a well-posed optimization problem.

---

## Numerical Structure

Both classes of problems lead to the same computational loop:

- solve the state equation;
- solve the adjoint equation;
- compute the reduced gradient;
- apply a descent, Newton, or active-set update;
- enforce the admissible constraints.

For Neumann boundary control, the discrete gradient contains a boundary mass
matrix because the control inner product is on $\Gamma$:

$$
M_\Gamma G
=
\alpha M_\Gamma U
+
M_\Gamma P_\Gamma.
$$

For Dirichlet boundary control, the discrete gradient involves the discrete
normal derivative of the adjoint.  This is more sensitive to mesh quality and
boundary approximation.

For parameter estimation, the gradient at a quadrature point has the form

$$
g_h
=
\alpha(q_h-q_{\mathrm{ref},h})
-
\nabla y_h\cdot\nabla p_h.
$$

Thus the quality of the gradient depends on the accuracy of both state and
adjoint gradients.  This is why coefficient identification problems often
require stabilization, mesh refinement, and careful treatment of noisy data.

---

## Summary

Boundary control and parameter estimation fit into the same reduced
optimization framework as distributed control, but the spaces change:

- Neumann control uses trace terms and gives gradients involving $\gamma p$;
- Dirichlet control uses liftings and gives dual boundary gradients involving
  $\partial_n p$;
- parameter estimation differentiates the PDE operator itself and gives
  gradients involving $\nabla y\cdot\nabla p$.

The main mathematical lessons are:

- the weak formulation determines the correct control space;
- Lax-Milgram gives the state and adjoint equations once coercivity and
  boundedness are checked;
- the reduced gradient is obtained by differentiating the Lagrangian with
  respect to the control or parameter after imposing state and adjoint
  stationarity;
- projection formulas are pointwise versions of variational inequalities;
- inverse problems require regularization not only for algorithms, but also
  for stability and well-posedness.

---

## Suggested Reading

### Tröltzsch

- Boundary control for elliptic PDEs
- Dirichlet and Neumann control
- Variational inequalities

### Manzoni-Quarteroni-Salsa

- Boundary observation and control
- Parameter identification
- Inverse problems

### De Los Reyes

- Reduced formulation
- Adjoint equations
- Regularization

# Reduced Elliptic OCPs

Numerical Methods for Optimal Control (NMOPT)

Luca Heltai (<luca.heltai@unipi.it>)

----

## Context

Lecture 3 gave an optimization toolbox for reduced problems

$$\min_u f(u).$$

Now the variable $u$ is a function, and each evaluation of $f$ or $\nabla f$
passes through a PDE solve.

----

## Model Problem

Distributed control of Poisson:

$$\min_{(y,u)} \frac12\|y-y_d\|_{L^2(\Omega)}^2 + \frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2$$

subject to

$$\begin{cases} -\Delta y=u+f & \text{in }\Omega,\\ y=0 & \text{on }\partial\Omega. \end{cases}$$

Data: $f,y_d\in L^2(\Omega)$, parameter $\alpha>0$.

----

## Strong to Weak

Start from

$$\int_\Omega (-\Delta y)v\,dx = \int_\Omega (u+f)v\,dx.$$

Integrate by parts:

$$\int_\Omega \nabla y\cdot \nabla v\,dx - \int_{\partial\Omega}\frac{\partial y}{\partial n}v\,ds = \int_\Omega (u+f)v\,dx.$$

If $v=0$ on $\partial\Omega$, the boundary term vanishes.

----

## Weak Formulation

Find $y\in H_0^1(\Omega)$ such that

$$\int_\Omega \nabla y\cdot \nabla v\,dx = \int_\Omega (u+f)v\,dx \qquad \forall v\in H_0^1(\Omega).$$

Why this space?

- first weak derivatives in $L^2$
- homogeneous Dirichlet condition built into the trace
- right-hand side meaningful for $u,f\in L^2(\Omega)$

----

## State Equation

Bilinear form and right-hand side:

$$a(y,v):=\int_\Omega \nabla y\cdot \nabla v\,dx, \qquad \ell_u(v):=\int_\Omega (u+f)v\,dx.$$

Lax-Milgram gives a unique $y\in H_0^1(\Omega)$ for each $u\in L^2(\Omega)$.

This defines the control-to-state map

$$S:L^2(\Omega)\to H_0^1(\Omega),\qquad u\mapsto y.$$

----

## Reduced Functional

Eliminate the state:

$$y=S(u).$$

Define

$$f(u)=J(S(u),u)=\frac12\|S(u)-y_d\|_{L^2(\Omega)}^2 + \frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2.$$

Reduced problem:

$$\min_{u\in L^2(\Omega)} f(u).$$

----

## Existence in Infinite Dimension

Finite-dimensional compactness is gone:

- closed and bounded is not strongly compact in $L^2$

What replaces it?

- coercivity
- weak compactness in reflexive spaces
- weak lower semicontinuity

----

## Coercivity

For this problem,

$$f(u)\ge \frac{\alpha}{2}\|u\|_{L^2(\Omega)}^2.$$

So minimizing sequences are bounded:

- if $f(u_n)\le C$ for large $n$, then

$$\frac{\alpha}{2}\|u_n\|_{L^2(\Omega)}^2\le C \quad\Longrightarrow\quad \|u_n\|_{L^2(\Omega)}\le \sqrt{\frac{2C}{\alpha}}.$$

----

## Weak Compactness

Since $L^2(\Omega)$ is reflexive, bounded sequences admit weakly convergent subsequences:

$$u_n\rightharpoonup \bar u \qquad \text{in }L^2(\Omega).$$

Weak convergence is enough if $f$ is weakly lower semicontinuous:

$$f(\bar u)\le \liminf_{n\to\infty} f(u_n).$$

For a minimizing sequence, this gives $f(\bar u)=\inf f$.

----

## Uniqueness

Uniqueness comes from strict convexity:

- $S$ is linear
- the tracking term is convex
- $\frac{\alpha}{2}\|u\|^2$ is strictly convex

Hence $f$ has a unique minimizer.

----

## Frechet Derivative

Compute from

$$f'(u)h=\lim_{t\to 0}\frac{f(u+th)-f(u)}{t}.$$

Since $S$ is linear,

$$S(u+th)=S(u)+tS(h).$$

After expanding the squares:

$$f'(u)h=(S(u)-y_d,S(h))_{L^2(\Omega)}+\alpha(u,h)_{L^2(\Omega)}.$$

----

## Why the Adjoint?

The problem is the term

$$(S(u)-y_d,S(h))_{L^2(\Omega)}.$$

The direction $h$ appears only through the state variation $S(h)$.

We want to rewrite this as an explicit inner product with $h$.

----

## Adjoint Equation

Find $p\in H_0^1(\Omega)$ such that

$$\int_\Omega \nabla p\cdot \nabla v\,dx = \int_\Omega (y-y_d)v\,dx \qquad \forall v\in H_0^1(\Omega).$$

Important point:

- in general this is the equation of the adjoint operator
- here the Poisson operator is self-adjoint, so the adjoint coincides with $-\Delta$

----

## State Variation

Let $z=S(h)$ be the solution of

$$\int_\Omega \nabla z\cdot \nabla v\,dx = \int_\Omega hv\,dx \qquad \forall v\in H_0^1(\Omega).$$

This is the first-order variation of the state in direction $h$.

----

## Reduced Gradient

Using the adjoint equation with test function $z$,

$$(y-y_d,z)=\int_\Omega \nabla p\cdot \nabla z\,dx.$$

Using the variation equation with test function $p$,

$$\int_\Omega \nabla z\cdot \nabla p\,dx=\int_\Omega hp\,dx.$$

Hence

$$f'(u)h=(p+\alpha u,h)_{L^2(\Omega)}.$$

----

## Gradient Formula

The reduced gradient is

$$\nabla f(u)=p+\alpha u.$$

This is the key computational identity:

- 1 state solve
- 1 adjoint solve
- 1 gradient evaluation

----

## Optimality System

For the unconstrained problem:

$$f'(\bar u)h=0 \qquad \forall h\in L^2(\Omega).$$

Equivalently,

$$\nabla f(\bar u)=0 \qquad\Longleftrightarrow\qquad \alpha \bar u+\bar p=0.$$

Together with state and adjoint:

$$\begin{cases} -\Delta \bar y=\bar u+f,\\ -\Delta \bar p=\bar y-y_d,\\ \alpha \bar u+\bar p=0. \end{cases}$$

----

## Algorithmic View

At iteration $k$:

1. solve state for $y_k$
2. solve adjoint for $p_k$
3. form $g_k=\alpha u_k+p_k$
4. update

$$u_{k+1}=u_k-\tau_k g_k$$

The methods of Lecture 3 now operate in function space.

----

## Transition

Next step:

1. add box constraints on the control
2. derive variational inequalities
3. obtain projection formulas and active/inactive sets

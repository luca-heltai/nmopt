# KKT Conditions in Finite Dimensions

Numerical Methods for Optimal Control (NMOPT)

Luca Heltai (<luca.heltai@unipi.it>)

----

## Problem Setting

We consider

$$\min_{u\in\mathbb{R}^n} f(u)\quad\text{s.t.}\quad \varphi_j(u)\ge 0,\; j=1,\dots,m.$$ 

Feasible set:

$$\mathcal U_{\mathrm{ad}}:=\{u\in\mathbb R^n:\varphi_j(u)\ge 0\ \forall j\}.$$ 

Active set at $\bar u$:

$$\mathcal A(\bar u):=\{j:\varphi_j(\bar u)=0\}.$$ 

----

## First-Order Admissible Directions

Linearized admissibility at $\bar u$:

$$\nabla\varphi_j(\bar u)\cdot d\ge 0\quad \forall j\in\mathcal A(\bar u).$$

Interpretation:

- negative product violates constraint $j$ at first order
- nonnegative product keeps constraint $j$ admissible at first order

----

## First-Order Necessary Condition

At a local constrained minimizer:

$$\nabla f(\bar u)\cdot d\ge 0\quad\text{for all admissible }d.$$ 

Why:

- if an admissible $d$ had $\nabla f(\bar u)\cdot d<0$
- then $\bar u+t d$ would decrease $f$ for small $t>0$
- contradiction with local optimality

----

## Geometric Example A

Single active inequality: halfspace tangent set and polar ray.

<img class="fit-figure" src="assets/02_polar_single_active_constraint.png" alt="Tangent halfspace and polar ray">

----

## Geometric Example B

First-quadrant cone and opposite-quadrant polar.

<img class="fit-figure" src="assets/02_polar_first_quadrant.png" alt="First quadrant cone and polar">

----

## Geometric Example C

Ray and its polar halfspace.

<img class="fit-figure" src="assets/02_polar_ray_halfspace.png" alt="Ray and polar halfspace">

----

## KKT System

For regular active constraints, there exist multipliers $\lambda_j\ge 0$ such that

$$\nabla f(\bar u)-\sum_{j=1}^m\lambda_j\nabla\varphi_j(\bar u)=0.$$ 

Together with feasibility/complementarity:

$$\varphi_j(\bar u)\ge 0,\quad \lambda_j\ge 0,\quad \lambda_j\varphi_j(\bar u)=0.$$ 

----

## Box Constraints

$$a_i\le u_i\le b_i\quad\Longleftrightarrow\quad u_i-a_i\ge 0,\; b_i-u_i\ge 0.$$ 

Componentwise first-order condition:

$$\partial_i f(\bar u)\begin{cases}\ge 0 & \bar u_i=a_i\\ \le 0 & \bar u_i=b_i\\ =0 & a_i<\bar u_i<b_i\end{cases}$$

This is the prototype for bound-constrained controls in PDE settings.

----

## Second Order

Lagrangian:

$$\mathcal L(u,\lambda):=f(u)-\sum_{j=1}^m\lambda_j\varphi_j(u).$$

Critical directions (minimal form):

$$\mathcal C(\bar u):=\{d:\nabla\varphi_j(\bar u)\cdot d\ge 0\ \forall j\in\mathcal A(\bar u),\ \nabla f(\bar u)\cdot d=0\}.$$ 

Necessary:

$$d^T\nabla^2_{uu}\mathcal L(\bar u,\lambda)d\ge 0\quad\forall d\in\mathcal C(\bar u).$$

----

## Summary

1. No admissible descent directions at a constrained minimum.
2. This gives gradient balance with active constraint normals.
3. That balance is KKT with nonnegative multipliers.
4. Second order checks Lagrangian curvature on critical directions.

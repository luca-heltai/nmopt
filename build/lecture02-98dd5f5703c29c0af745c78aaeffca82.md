# KKT Conditions in Finite Dimensions

## Overview

In this lecture we show the core ideas needed to understand inequality constrained first-order and second-order optimality in finite dimensions.

We focus on the problem
$$
\min_{u\in\mathbb{R}^n} f(u)
\quad\text{s.t.}\quad
\varphi_j(u)\ge 0,\; j=1,\dots,m,
$$
with $f,\varphi_j\in C^1(\mathbb{R}^n)$.

The take-home message is simple:

1. at a constrained minimum there is no admissible descent direction;
2. this forces the objective gradient to balance active constraint normals;
3. this balance is exactly the KKT system;
4. second order is curvature of the Lagrangian along directions that remain tangent to active constraints.
5. with the convention $\mathcal L=f-\sum_j\lambda_j\varphi_j$ and $\varphi_j\ge 0$, active multipliers are nonnegative, and they vanish for inactive constraints (complementarity).

---

## 1. Problem Setting and Active Constraints

Define the feasible set
$$
\mathcal U_{\mathrm{ad}}:=\{u\in\mathbb R^n:\varphi_j(u)\ge 0\ \forall j\}.
$$

Let $\bar u\in\mathcal U_{\mathrm{ad}}$ be a local minimizer.

Not all constraints matter equally at $\bar u$.
Only the constraints that are exactly at the boundary can influence first-order optimality.
So we define the active set
$$
\mathcal A(\bar u):=\{j\in\{1,\dots,m\}:\varphi_j(\bar u)=0\}.
$$

Inactive constraints satisfy $\varphi_j(\bar u)>0$ and do not contribute at first order. For example: if we have a contact constraint between two bodies (that cannot penetrate), then only the constraints corresponding to points in contact are active, while the others are inactive.

---

## 2. Admissible First-Order Directions

Near $\bar u$, a direction $d\in\mathbb R^n$ is first-order admissible if it does not decrease active constraints.
Linearizing active constraints gives
$$
\nabla\varphi_j(\bar u)\cdot d\ge 0
\quad\forall j\in\mathcal A(\bar u).
$$

This condition is the local geometric ingredient we need.

Interpretation:

- if $\nabla\varphi_j(\bar u)\cdot d<0$, then moving along $d$ pushes outside feasibility for constraint $j$;
- if $\nabla\varphi_j(\bar u)\cdot d\ge 0$, then constraint $j$ is not violated at first order.

---

## 3. First-Order Necessary Condition

At a local minimizer, every admissible first-order direction must be non-descent:
$$
\nabla f(\bar u)\cdot d\ge 0
\quad\text{for every } d \text{ such that }
\nabla\varphi_j(\bar u)\cdot d\ge 0\ \forall j\in\mathcal A(\bar u).
$$

Why this is unavoidable:

- if there existed one admissible $d$ with $\nabla f(\bar u)\cdot d<0$,
- then a small step $\bar u+t d$ (with $t>0$ small) would decrease the objective while remaining feasible at first order,
- contradicting local minimality.

This is the true origin of KKT.

---

## 4. From First-Order Geometry to KKT

For regular active constraints (the standard finite-dimensional nondegenerate case),
there exist multipliers $\lambda_j\ge 0$ such that
$$
\nabla f(\bar u)-\sum_{j\in\mathcal A(\bar u)}\lambda_j\nabla\varphi_j(\bar u)=0.
$$

Extending by $\lambda_j=0$ on inactive constraints yields the standard stationarity form
$$
\nabla f(\bar u)-\sum_{j=1}^m\lambda_j\nabla\varphi_j(\bar u)=0.
$$

Together with primal feasibility and complementarity, we obtain
$$
\begin{aligned}
\nabla f(\bar u)-\sum_{j=1}^m\lambda_j\nabla\varphi_j(\bar u) &= 0,\\
\varphi_j(\bar u) &\ge 0,\quad j=1,\dots,m,\\
\lambda_j &\ge 0,\quad j=1,\dots,m,\\
\lambda_j\,\varphi_j(\bar u) &= 0,\quad j=1,\dots,m.
\end{aligned}
$$

Complementarity meaning:

- if $\varphi_j(\bar u)>0$ (inactive), then $\lambda_j=0$;
- if $\lambda_j>0$, then necessarily $\varphi_j(\bar u)=0$ (active).

```{admonition} Regularity Note
:class: note

In finite dimensions, this KKT form is valid under standard nondegeneracy assumptions on active constraints.
For box constraints, this regularity is automatic in the usual formulation.
```

---

## 5. Box Constraints: The Prototype Case

Consider
$$
\min_{u\in\mathbb R^n} f(u)
\quad\text{s.t.}\quad
a_i\le u_i\le b_i,\quad i=1,\dots,n.
$$

Equivalent inequality form:
$$
u_i-a_i\ge 0,
\qquad
b_i-u_i\ge 0.
$$

The componentwise first-order condition becomes
$$
\begin{cases}
\partial_i f(\bar u)\ge 0 & \text{if } \bar u_i=a_i,\\
\partial_i f(\bar u)\le 0 & \text{if } \bar u_i=b_i,\\
\partial_i f(\bar u)=0 & \text{if } a_i<\bar u_i<b_i.
\end{cases}
$$

This is the exact finite-dimensional model behind bound-constrained control variables in PDE-constrained optimization.

---

## 6. Second-Order Conditions

Define the Lagrangian
$$
\mathcal L(u,\lambda):=f(u)-\sum_{j=1}^m\lambda_j\varphi_j(u).
$$

After first-order conditions are satisfied, the next question is curvature along first-order feasible/tangent directions.
A practical critical set is
$$
\mathcal C(\bar u):=\{d:\nabla\varphi_j(\bar u)\cdot d\ge 0\ \forall j\in\mathcal A(\bar u),\ \nabla f(\bar u)\cdot d=0\}.
$$

Second-order necessary condition:
$$
d^T\nabla^2_{uu}\mathcal L(\bar u,\lambda)d\ge 0
\quad\forall d\in\mathcal C(\bar u).
$$

Second-order sufficient condition (strict local minimality):
$$
d^T\nabla^2_{uu}\mathcal L(\bar u,\lambda)d>0
\quad\forall d\in\mathcal C(\bar u)\setminus\{0\}.
$$

Interpretation:

- first order kills linear descent along admissible directions;
- second order checks that curvature is positive on directions where first order is neutral.

---

## 7. Summary

At a constrained local minimizer:

1. no admissible first-order direction can decrease $f$;
2. this forces a balance between objective gradient and active constraint normals;
3. that balance is KKT with nonnegative multipliers and complementarity;
4. second order is positivity of Lagrangian curvature on critical directions.

This is the essential finite-dimensional picture we need before moving to PDE-constrained settings.

<!-- FOOTER START -->
<iframe src="../slideshow/slides02.html" width="100%" height="800px" style="border: none;"></iframe>

---

```{admonition} 🎬 View Slides
:class: tip

**[Open slides in full screen](../slideshow/slides02.html)** for the best viewing experience.
```
<!-- FOOTER END -->

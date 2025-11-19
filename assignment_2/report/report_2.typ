#import "@preview/physica:0.9.6": dv
#import "@preview/physica:0.9.4": curl, grad, tensor, pdv, dv, eval, TT
#import "preamble.typ": mref, mathformatter


//! Styling
#set page(
  numbering: "1/1",
  margin: (
    top: 1.5cm,
    bottom: 1.5cm,
    left: 2.0cm,
    right: 3.0cm,
  )
)

#set heading(numbering: "1.a.1)")
#set math.equation(numbering: "(1)")


//! Declarations
#let wider = h(3em)
#let vv = mathformatter(underbars: 0, bold: true, upright: true)
#let mm = mathformatter(underbars: 0, bold: true, upright: true)
#set math.mat(delim:"[")
#set math.vec(delim:"[")

#let Span = math.op("Span")

//! Cover Page

#align(center, text(17pt)[
  *01617 Introduction to Dynamical Systems (E25) - Assignment 2*
])

#strong[Authors]:
- Jeppe Klitgaard `<s250250@dtu.dk>`
- Yunis Wirkus `<s250700@dtu.dk>`

Date: 2025-11-28

#pagebreak()


//! Content
= Problem 1 — The Duffing Equation <sec:1>

We are given another flavour of the Duffing equation:
$
  dot.double(x)_1 &= f(x_1, x_2) = mu x_1 - x_1^3 \
$

or equivalently as a system of first order ODEs:

$
  dot(x)_1 &= f_1(x_1, x_2) = x_2 \
  dot(x)_2 &= f_2(x_1, x_2) = mu x_1 - x_1^3
$ <eq:1_duffing_original>

with $mu$ and $(x_1, x_2) in RR^2$. We use the suppressing notation for $x_(1,2)(t): RR^+ arrow RR$, without indicating the dependence on time $t$ explicitly. Further, we use the common dot notation $dot(x)_(1,2) = frac(d x_(1,2), d t)$, analogously for higher order derivatives.

This is a common form of the Duffing equation, describing a dampened and driven oscillator with a non-linear restoring force @ REF. 

== Finding critical points <sec:1a>

We wish to find values of $x_1, x_2$ such that $f_1(x_1, x_2) = dot(x)_1 = 0$ and $f_2(x_1, x_2) = dot(x)_2 = 0$. Thus, we find that: 

$
  dot(x)_1 &= 0 quad &&arrow.double quad x_2 = 0 \
  dot(x)_2 &= x_1(mu - x_1^2) = 0 quad &&arrow.double quad x_1 in {0, +sqrt(mu), -sqrt(mu)} \
$

From this, we can conclude that the critical points of the system are given by:

$
  x^* in cases(
    {(0,0)} quad &forall mu in RR,
    {(-sqrt(mu), 0), (0,0), (sqrt(mu), 0)} quad &forall mu > 0,
    {(0,0), (-i sqrt(|mu|), 0), (i sqrt(|mu|), 0)} quad &forall mu < 0,
  )
$<eq:1a_all_critical_points>

We will only consider real critical points in the following analysis. Hence we will investigate the system behaviour around the critical points 

$
  x^* in cases(
    {(-sqrt(mu), 0), (0,0), (sqrt(mu), 0)} quad &mu > 0,
    {(0,0)} quad &mu ≤ 0
  )
$<eq:1a_critical_points>

for easier reference later, we enumerate the critical points as follows: 
$
  x^*_0 = (0,0) quad x^*_1 = (sqrt(mu), 0) quad x^*_2 = (-sqrt(mu), 0)
$

\
== Stability of critical points<sec:1b>

We will analyse the stability of the critical points found in @eq:1a_critical_points by computing the Jacobian matrix of the system at each critical point and evaluating its eigenvalues. Here we will exclude the case for $mu=0$.

The Jacobian matrix $J = D_f$ of the system/vectorfield $bold(f) = (f_1, f_2)^T$ from @eq:1_duffing_original given by:

$
  D_f(x_1, x_2) = 
  mat(
    partial(f_1, x_1), partial(f_1, x_2); 
    partial(f_2, x_1), partial(f_2, x_2)
  ) 
  = mat(
    0, 1; 
    mu - 3 x_1^2, 0
  )
$<eq:1b_jacobian>

At the *critical point $x^*_0 = (0,0)$*, we find the Jacobian:
$
  D_f(x^*_0) = 
  mat(
    0, 1; 
    mu, 0
  ) 
  quad arrow.double quad 
  p_(x^*_0)(λ) = 
  det(
    mat(
      λ, -1
      ; -mu, λ
    )
  ) = λ^2 - mu = 0
  quad arrow.double quad
  λ_(1,2) = ± sqrt(mu)
$

As the eigenvalues are depended on the parameter $mu$, we consider two cases. As afore mentioned, excluding the degenerate case $mu=0$, where $λ_(1,2) = 0$ is a algebraic multiplicity of 2, thus a critical case with non-hyperbolic behaviour @ REF, hence initially inconclusive.

1. If $mu>0$: eigenvalues are real, opposite signs $lambda_1=+ sqrt(mu), lambda_2=-sqrt(mu)$. That is a *saddle point* with one stable (negative eigenvalue $lambda_2$) and one unstable direction. Saddles are hyperbolic but unstable as trajectories are repelled in one direction and attracted in another.

2. If $mu < 0$: eigenvalues are pure imaginary $lambda = plus.minus i sqrt(|mu|)$. This indicates that the linearization via the Jacobian only gives a _linear center_ (center for the linearization). For conclusions relating back to the non-linear system we need more information.

Purely imaginary eigenvalues indicate _neutral stability_ in the linearized system, meaning that trajectories neither converge nor diverge from the critical point, but rather orbit around it. However, this does not guarantee stability in the non-linear system, as non-linear effects can lead to different behaviours.

At the *critical points $x^*_(1,2) = (plus.minus sqrt(mu), 0)$*, which only exists for $mu > 0$, we find the Jacobian:

$
  D_f(x^*_(1,2)) &= 
  mat(
    0, 1; 
    mu - 3(mu), 0
  ) =
  mat(
    0, 1; 
    -2 mu, 0
  ) \
  &arrow.double quad 
  p_(x^*_(1,2))(λ) = 
  det(
    mat(
      λ, -1
      ; 2 mu, λ
    )
  ) = λ^2 + 2 mu = 0 \
  &arrow.double quad
  λ_(1,2) = ± i sqrt(2 mu)
$

Again, these eigenvalues are purely imaginary, indicating _neutral stability_ in the linearized system. As before, this does not guarantee stability in the non-linear system, and further analysis would be required to determine the actual stability of these critical points.

Summing up, we can only definitely classify one of the critical points $x^*_0$ in @eq:1a_critical_points as a saddle point (unstable) for $mu > 0$. The other critical points exhibit purely imaginary eigenvalues, indicating _neutral stability_ in the linearized system, but their actual stability in the non-linear system remains inconclusive without further analysis.
For $mu < 0$, only $x^*_0$ exists and exhibits purely imaginary eigenvalues, again indicating _neutral stability_ in the linearized system, but inconclusive for the non-linear system.

\
== Sketching the Phase Portrait & Flow <sec:1c>

\
*Clarify: Sketch flow 'in' the phase portrait? Meaning plotting trajectories of the solution in the phase portrait?*

#figure(
  grid(
    columns: 2,
    gutter: 1em,
    grid.cell(
      image("../img/1_c_mu_gt_1.png",      
      height: 7cm)
    ),
    grid.cell(
      image("../img/1_c_mu_lt_1.png",      
      height: 7cm)
    )
  ),
  caption: [Sketch of the phase plane of vectorfield $bold(f)$ from @eq:1_duffing_original for $mu > 0$ (left) and $mu < 0$ (right). Critical points are marked in red (saddles) and blue (centers). Arrows indicate the flow of the system.],
)<fig:1c>

- insert hand-drawn sketches here
- comment on the plot?!

\
== Analytical Homoclinic Orbit (1st Quadrant) <sec:1d>

Before we begin discussions, it helps to formulate a more intuitive definition of _homoclinic orbits_ (@ REF). A homoclinic orbit is a trajectory in the phase space of a dynamical system that starts and ends at the same saddle point. In other words, it is a path that leaves a saddle point along its unstable manifold and then returns to the same saddle point along its stable manifold.
In this case, it does a tour around other critical points in between.

In the phase portrait for $mu > 0$ (left figure in @fig:1c), we can observe homoclinic orbits associated with the saddle point $x^*_0$ at the origin $(0,0)$. These orbits leave the saddle point along its unstable manifold and return to it along its stable manifold, forming closed loops around the other critical points $x^*_(1,2)$ at $(plus.minus sqrt(mu), 0)$.

Now we can determine a homeomorphic function $h in C^1$, that can identify / trace these homoclinic orbits in the phase space. It is important to look at the function by each quadrant of the phase plane to prevserve injectivity. This functiion $h: RR arrow RR$ will describe a graph in the phase plane, i.e. $x_2 = h(x_1)$.

To derive a first order differential equation for the solution of the graph $x_2 = h(x_1)$, along that graph curve, it must satisfy:

$
  frac(d x_2, d x_1) = frac(dot(x)_2, dot(x)_1) = frac(f_2(x_1, x_2), f_1(x_1, x_2)) = frac(mu x_1 - x_1^3, x_2) = frac(mu x_1 - x_1^3, h(x_1)) = frac(d h, d x_1)
$<eq:1d_homoclinic_graph_condition>

Hence, we can derive the first order ODE for $h(x_1)$:

$
  frac(d h, d x_1) = dot(h)(x_1) = (mu x_1 - x_1^3) / h(x_1)
$<eq:1d_homoclinic_ode>

This is homoclinic at the origin, as $h(0) = 0$ and the first quadrant satisfies $h(x_1) > 0$ for $x_1 > 0$.

@eq:1d_homoclinic_ode is a separable first order ODE for $h(x_1)$, which we can solve by separating variables (multiplying LHS by $h(x_1)$ and RHS by $d x_1$ respectively) and integrating both sides:

$
  frac(d h, d x_1) &= (mu x_1 - x_1^3) / h(x_1) \
  arrow.l.r.double quad quad quad quad 
  h(x_1) d h &= (mu x_1 - x_1^3) d x_1 \
  arrow.l.r.double quad quad quad 
  ∫ h(x_1) d h &= ∫ (mu x_1 - x_1^3) d x_1 \
  arrow.l.r.double quad quad quad quad 
  1/2 h(x_1)^2 &= 1/2 mu x_1^2 - 1/4 x_1^4 + C \
  arrow.l.r.double quad quad quad quad 
  h(x_1)^2 &= mu x_1^2 - 1/2 x_1^4 + 2C
$<eq:1d_homoclinic_solution_integral>

Where $C in RR$ is the integration constant.

Since we already established that the homoclinic orbit passes through the origin, we know that $h(0)=0$. And as $x_1$ for all $x^*_(0,1,2)$, we plug in $h(0)^2 = 2C = 0$ from @eq:1d_homoclinic_solution_integral to solve for $C=0$. Therefore:

$
  h(x_1)^2 &= mu x_1^2 - 1/2 x_1^4 \
  h(x_1) &= plus.minus x_1 sqrt(mu - frac(x_1^2, 2)) \
$<eq:1d_homoclinic_solution>

For the first quadrant, we take the positive root, yielding the homoclinic orbit function:

$
  h(x_1)^+ = x_1 sqrt(mu - frac(x_1^2, 2)) quad quad forall x_1 in [0, sqrt(2 mu)]
$

Note that we can interpret the maximum possible value for $x_1 = sqrt(2 mu)$ (for the root to ne non-negative) as the turning point of the homoclinic orbit in the first quadrant, after which it returns to the origin along the stable manifold along the orbit on the fourth quadrant ($x_1 >0, x_2 < 0$). Thus it starts at the origin $h(x_1=0)=0$ and transitions into the fourth quadrant at $h(x_1 = sqrt(2 mu))=0$, see @fig:1c.

\
== Plotting the Homoclinic Orbit (1st quadrant) <sec:1e>

#figure(
  image("../img/1_e_homoclinic_orbit.png",      
  height: 8cm),
  caption: [Plot of the homoclinic orbit $x_2 = h(x_1)$ from @eq:1d_homoclinic_solution for $mu>0$. The orbit starts and ends at the saddle point at the origin $(0,0)$, forming a closed loop around the other critical points.],
)<fig:1e>

- comment on the plot?!

\
== Analytical Homoclinic Orbit (4th Quadrant) <sec:1f>

Given the nature of @eq:1d_homoclinic_solution, for the fourth quadrant of the phase plot ($x_1 > 0, x_2 < 0$), it is merely taking the negative root of the homoclinic orbit function, i.e.: 

$
  h(x_1)^- = - x_1 sqrt(mu - frac(x_1^2, 2)) quad quad forall x_1 in [0, sqrt(2 mu)]
$

The resulting graph of $h$ can be seen in @fig:1e in light-green, which is a mirror of the first quadrant homoclinic orbit across the $x_1$-axis. All statements concerning turning points and behaviour apply analogously.


#pagebreak()
= Problem 2 — Bifurcation <sec:2>
\
== Bifurcation parameter & diagram <sec:2a>

A bifurcation parameter is a parameter in a dynamical system that, when varied, causes a change in the system's behavior or structure. Bifurcations often lead to the emergence or disappearance of critical points, periodic orbits, or changes in stability at those points.

As our analysis in @sec:1 showed, the parameter $mu$ in the Duffing equation from @eq:1_duffing_original acts as a bifurcation parameter. Specifically, we observed an emergence of two additional critical points $x^*_(1,2)$ at $(plus.minus sqrt(mu), 0)$ when $mu$ transitions from negative to positive values. This indicates a bifurcation occurring at the critical value $mu=0$.

#figure(
  image("../img/2_a_bifurcation_param_diagram.png",      
  height: 8cm),
  caption: [Bifurcation diagram illustrating the emergence of critical points in the Duffing equation from @eq:1_duffing_original as the bifurcation parameter $mu$ varies.],
)<fig:2a>

- comment on the plot?!

\
== Type of Bifurcation <sec:2b>
The bifurcation observed in @fig:2a from Duffing equation from @eq:1_duffing_original as the parameter $mu$ crosses zero is the prototype of a *supercritical pitchfork bifurcation* (@ REF). 
In this supercritical pitchfork bifurcation, as the bifurcation parameter $mu$ passes through a transition value (here, $mu=0$), a single center stable critical point at the origin $(0,0)$ for $mu<0$ becomes an unstable saddle point for $mu>0$, and two new stable critical points (at $(plus.minus sqrt(mu), 0)$) emerge symmetrically on either side of the original point.

\
== Hyperbolicity for bifurcation transition $mu=0$ <sec:2c>
At $mu=0$ all eigenvalues of the Jacobian at the critical point 

$
  x^*_0 arrow D_f(x^*_0) = D_f(0,0) = mat(0, 1; 0, 0)
$

are zero, i.e. $lambda_(1,2) = 0$, it has algebraic multiplicity of two at $0$. This indicates that the critical point is non-hyperbolic at the bifurcation point $mu=0$. Hence, the linearization theorem that requires hyperbolicity does not apply here, and we cannot draw conclusions about the stability of the critical point based solely on the linearized system.

\
== Lyapunov Stability for non-hyperbolic case <sec:2d>
As we determined in @sec:2c, the critical point $x^*_0$ at the bifurcation point $mu=0$ is non-hyperbolic and we cannot conclude stability from the linearized system alone. Hence, we can use a Lyapunov function to analyse the stability of this critical point. A very common first guess is the simple quadratic form $V(x_1, x_2) = x_1^2 + x_2^2$. Another common choice is the energy function; in this case a natural candidate, since the system can be interpreted as a mechanical oscillator.

$
  V(x_1, x_2) = frac(1,2) x_2^2 + frac(1,4) x_1^4
$<eq:2d_lyapunov_function>

The Lyapunov function needs to be positive definite, i.e. $V(x_1, x_2) > 0$ for all $(x_1, x_2) ≠ (0,0)$ and $V(0,0) = 0$. This is satisfied for @eq:2d_lyapunov_function, as both terms are squares and thus non-negative, and the only point where $V$ is zero is at the origin.

Next, we compute the time derivative of $V$ along the trajectories of the system $bold(f)$ for $mu=0$:

$
  dot(V)(x_1, x_2) &= partial(V, x_1) dot(x)_1 + partial(V, x_2) dot(x)_2 \
  &= x_1^3 dot(x)_1 + x_2 dot(x)_2 \
  &= x_1^3 x_2 + x_2 (mu - x_1^3) \
  &= x_1^3 x_2 + x_2 (- x_1^3) \
  &equiv 0
$

Since $V$ is positive definite and $dot(V)≡0$, we can interpret that the energy in the system stays constant along each trajectory. That means, trajectories are contained within their (energy-)level-set of the Lyapunov function. 
Thus, the origin $(x_1, x_2) = (0,0)$ at $mu=0$ is Lyapunov (stable): small initial conditions remain small for all time. 
However, because $dot(V)$ is identically zero (not negative definite), one cannot conclude asymptotic stability. In fact, the origin is not asymptotically stable (solutions near origin are periodic orbits—center behavior, similar to the linear centers $x^*_(1,2)$ with purely imaginary Eigenvalues in the linearized system).

So, the origin is a stable center (Lyapunov stable) but not asymptotically stable.


\
== Contour Plot of Lyapunov Function <sec:2e>

#figure(
  image("../img/2_e_contourV_phase.png",      
  height: 8cm),
  caption: [Contour plot of the Lyapunov function $V(x_1, x_2)$ from @eq:2d_lyapunov_function for $mu=0$. The contours represent level sets of constant energy in the system. Overlay over the phase plot with vector field arrows.],
)<fig:2e>

As clearly visible in @fig:2e, the contours of the Lyapunov function $V(x_1, x_2)$ from @eq:2d_lyapunov_function represent not only level sets of constant energy but also trajectories as solutions of the system @eq:1_duffing_original for $mu=0$. This is consistent with our earlier analysis in @sec:2d, where we found that the time derivative of the Lyapunov function $dot(V)$ is identically zero along trajectories, indicating that energy remains constant.

*REVIEW*
Further we observe, that the level sets along with trajectories of the vector field form elliptical closed orbits around the origin. The elliptical shape could be attributed to the stretch of the system along the $x_1$-axis, due to the quartic term in $V$ (@eq:2d_lyapunov_function), compared to only the quadratic term for $x_2$.


\
== Stability for $(0,0)$ at $mu=0$<sec:2f>
The required Theorem from @perko[Chapther 2.11, p.151, Theorem 2] states:

For when matrix $A$ of a linearized system has two zero eigenvalues, i.e. $det(A) = tr(A) = 0$, but $A≠0$, the system can be put into a 'normal' form:

$
  dot(x) &= y \
  dot(y) &= a_k x^k [1+h(x)] + b_n x^n y [1+g(x)] + y^2 R(x,y) \
$<eq:2f_perko_normal_form>

where $h(x), g(x), R(x,y)$ are analytic in the neighborhood of the origin with $h(0) = g(0) = 0$ and $k ≥ 2, n ≥ 1$.

#quote(attribution: [Perko])[
  *Theorem 2* \
  Let $k= 2m+1$ with $m ≥1$ in @eq:2f_perko_normal_form and let $lambda = b_n^2 + 4(m+1) a_k$. Then if $a_k > 0$, the origin is a topological saddle. If $a_k < 0$, the origin is \
  (1) a focus or a center if $b_n=0$ and also if $b_n≠0$ and $n>m$ or if $n=m$ and $lambda < 0$, \
  (2) a node if $b_n ≠ 0$, $n$ is an even number and $n < m$ and also if $b_n ≠ 0$, $n$ is an even number, $n =m$ and $lambda ≥0$ and \
  (3) a critical point with an elliptic domain if $b_n ≠ 0$, $n$ is an odd number and $n < m$ and also if $b_n ≠ 0$, $n$ is an odd number, $n =m$ and $lambda ≥ 0$.
]

From this theorem, we can conclude that our system @eq:1_duffing_original at the critical point $(0,0)$ for $mu=0$ can be transformed into the normal form @eq:2f_perko_normal_form with: $k=3, a_k=-1, b_n=0, h(x)=0, R(x,y)=0$

Thus, we have a focus or a center according to case (1) in the theorem, which is consistent with our earlier findings in @sec:2d using the Lyapunov function, where we concluded that the origin is a stable center (Lyapunov stable) but not asymptotically stable.

#bibliography("refs.bib")

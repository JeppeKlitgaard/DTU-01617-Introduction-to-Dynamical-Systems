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

#show heading: set block(above: 1.5em)

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
$ <eq:1_duffing>

or equivalently as a system of first order ODEs:

$
  dot(x)_1 &= f_1(x_1, x_2) = x_2 \
  dot(x)_2 &= f_2(x_1, x_2) = mu x_1 - x_1^3
$ <eq:1_duffing_original>

with $mu∈ℝ$ and $(x_1, x_2) in E = RR^2$. We use the suppressed notation for $x_(1,2)(t): RR^+ arrow RR$ such that $dot(x)_(1,2) = dv(x_(1, 2), t)$, with higher order derivatives being indicated by additional dots.

// TODO: REF
This is a common form of the Duffing equation, describing a damped and driven oscillator with a non-linear restoring force @ REF.

== Finding critical points <sec:1a>

We wish to find values of $x_1, x_2$ such that $f_1(x_1, x_2) = dot(x)_1 = 0$ and $f_2(x_1, x_2) = dot(x)_2 = 0$. Thus, we find that:

$
  dot(x)_1 &= 0 quad &&arrow.double quad x_2 = 0 \
  dot(x)_2 &= x_1(mu - x_1^2) = 0 quad &&arrow.double quad x_1 in {0, +sqrt(mu), -sqrt(mu)} \
$

From this, we can conclude that the critical points of the system are given by:

$
  vv(x^*) in cases(
    {(0,0)} quad &forall mu in RR,
    {(-sqrt(mu), 0), (0,0), (sqrt(mu), 0)} quad &forall mu > 0,
    {(0,0), cancel((-i sqrt(|mu|), 0)), cancel((i sqrt(|mu|), 0))} quad &forall mu < 0,
  )
$<eq:1a_all_critical_points>

Where the $μ<0$ case produces two complex critical points which are discarded due to the domain restriction $(x_1, x_2) ∈ ℝ^2$. Hence we will investigate the system behaviour around the critical points

$
  vv(x^*) in cases(
    {(-sqrt(mu), 0), (0,0), (sqrt(mu), 0)} quad &mu > 0,
    {(0,0)} quad &mu ≤ 0
  )
$<eq:1a_critical_points>

for easier reference later, we enumerate the critical points as follows:
$
  vv(x^*_0) = (0,0) quad vv(x^*_1) = (sqrt(mu), 0) quad vv(x^*_2) = (-sqrt(mu), 0)
$ <eq:1a_critical_points_ref>


== Stability of critical points<sec:1b>

We will analyse the stability of the critical points found in @eq:1a_critical_points by computing the Jacobian matrix of the system at each critical point and evaluating its eigenvalues. Here we will exclude the case for $mu=0$.

The Jacobian matrix $A = D vv(f(vv(x)))$ of the system/vector field $vv(f) = (f_1, f_2)^TT ∈ C^1 (E)$ from @eq:1_duffing_original given by:

$
  D vv(f)(x_1, x_2) =
  mat(
    pdv(f_1, x_1), pdv(f_1, x_2);
    pdv(f_2, x_1), pdv(f_2, x_2)
  )
  = mat(
    0, 1;
    mu - 3 x_1^2, 0
  )
$<eq:1b_jacobian>

=== Case: $μ < 0$ <sec:1b_mu_neg>

We consider the case $μ < 0$ first, which only has a single critical point at the origin $vv(x^*_0) = (0,0)$,
at which we linearise the system by finding the Jacobian matrix:
$
  A = D vv(f(vv(x^*_0))) = mat(
    0, 1;
    mu, 0
  )
$

Which we obtain the eigenvalues $λ_i$ of by solving the characteristic polynomial:
$
  p_(x^*_0)(λ) = det(A - λ I) = λ^2 - μ = 0 ⇒ λ_(1,2) = ± sqrt(μ) = ± i sqrt(|μ|)
$

Which form a pure imaginary complex conjugate pair, indicating that the critical point is a _center_ for the linearised system. We seek to employ @perko[Chapter 2.10, Theorems 5 & 6], which requires that we investigate the symmetry of the vector field $vv(f)$. Using the definition of symmetry from @perko[Chapter 2.10, Definition 6], we consider a transformation $(t, x_1) → (-t, -x_1)$ applied to @eq:1_duffing:
$

  dot.double(x)_1 &= μ x_1 - x_1^3 ⇔ - dot.double(x)_1 &= - μ x_1 + x_1^3
$

Which indicates that the vector field is symmetric with respect to the $x_2$ axis. Thus, by @perko[Chapter 2.10, Theorem 6], we conclude that the critical point $vv(x^*_0) = (0,0)$ is also a _center_ of the non-linear system.

In order to gauge the stability of the critical point, we employ the stability theory presented in @perko[Chapter 2.9], where Definition 1 defines stability as:

*Definition 1*: Let $vv(ϕ)_t$ denote denote the flow of the differential equation $dot(vv(x)) = vv(f(vv(x)))$ defined for all $t ∈ ℝ$. An equilibrium point $vv(x)_0$ is _stable_ if for all $ε > 0$ there exists a $δ > 0$ such that for all $vv(x) ∈ N_ε (vv(x)_0)$ we have
$
  vv(ϕ)_t (vv(x)) ∈ N_ε (vv(x)_0)
$

The equilibrium point $vv(x)_0$ is _unstable_ if it is not stable. And $vv(x)_0$ is _asymptotically stable_ if it is stable and if there exists a $δ>0$ such that for all $vv(x) ∈ N_δ (vv(x)_0)$ we have
$
  lim_(t → ∞) vv(ϕ)_t (vv(x)) = vv(x_0)
$

We start by constructing a Liapunov function $V(vv(x))$ for the system, which by @perko[Chapter 2.9, Theorem 3]
must satisfy:
$
  V(vv(x_0)) = 0 wider and wider V(vv(x)) > 0 "if" vv(x) ≠ vv(x)_0
$ <eq:1_liapunov_requirements>

We construct such a function by considering energy conservation of the system and production of the Hamiltonian, $H$, which by definition satisfies:
$
  pdv(H, x_2) = dot(x_1) = x_2 wider pdv(H, x_1) = dot(x_2) = -μ x_1 + x_1^3
$
Such that:
$
  H(x_1, x_2) &= ∫ x_2 dif x_2 = 1/2 x_2^2 + h_1 (x_1)\
  H(x_1, x_2) &= ∫ -μ x_1 + x_1^3 dif x_1 = - 1/2 μ x_1^2 + 1/4 x_1^4 + h_2 (x_2)\
$

Allowing the construction of the Hamiltonian, which we let be our Liapunov function:
$
  V (vv(x)) = H(x_1, x_2) = 1/2 x_2^2 - 1/2 μ x_1^2 + 1/4 x_1^4 + cancel(c)
$ <eq:1_hamiltonian>

Which satisfies the necessary conditions for a Liapunov function given in @eq:1_liapunov_requirements by inspection. We consider the temporal derivative of $V (vv(x))$:
$
  dot(V) (vv(x)) = 0
$

Which by @perko[Chapter 2.9, Theorem 3] implies that the critical point is _stable_, but not _asymptotically stable_. This is sometimes referred to as _Liapunov stablility_ or _neutral stability_.

=== Case: $μ > 0$

We then consider the case $μ > 0$, which has three critical points at $vv(x^*_0) = (0,0)$, $vv(x^*_1) = (sqrt(mu), 0)$, and $vv(x^*_2) = (-sqrt(mu), 0)$.

We first compute the Jacobian at $vv(x^*_0) = (0, 0)$ and the eigenvalues of the linearised system:
$
  A = D vv(f(vv(x^*_0))) =
  mat(
    0, 1;
    mu, 0
  )
$

$
  p_(x^*_0)(λ) = det(A - λ I) = λ^2 - μ = 0 ⇒ λ_(1,2) = ± sqrt(μ)
$

Where the eigenvalues are real and of opposite signs, indicating the critical point is a _saddle point_ for the linearised system. Additionally, we note that saddle points are hyperbolic as their eigenvalues are real and non-zero, and that they are _unstable_ by @perko[Chapter 2.9, Theorem 2].

As such, we are able to leverage @perko[Chapter 2.10, Theorem 3] to conclude that the critical point $vv(x^*_0) = (0,0)$ is also a _saddle point_ of the non-linear system.

We may now consider the other two critical points $vv(x^*_(1,2)) = (±sqrt(mu), 0)$, for which we compute the Jacobian and eigenvalues of the linearised system:
$
  A = D vv(f(vv(x^*_(1,2)))) =
  mat(
    0, 1;
     μ - 3 x_1^2, 0
  ) =
  mat(
    0, 1;
    μ - 3μ, 0
  ) =
  mat(
    0, 1;
    -2 μ, 0
  )
$
$
  p_(x^*_(1,2))(λ) = det(A - λ I) = λ^2 + 2 μ = 0 ⇒ λ_(1,2) = ± i sqrt(2 μ)
$

By the same logic as in @sec:1b_mu_neg, we find that these are _centers_ of the linearised system.

In order to determine their stability, seek to prove that they are local minimisers by computing the Hessian:
$
  D^2 V (vv(x)) = mat(
    pdv(V, x_1, 2), pdv(V, x_1, x_2);
    pdv(V, x_2, x_1), pdv(V, x_2, 2)
  ) =
  mat(
    -1μ + 3μ, 0;
    0, 1
  ) =
  mat(
    2μ, 0;
    0, 1
  )
$

Which is positive definite for $μ > 0$, implying that these critical points are local minimisers of the Liapunov function from @eq:1_hamiltonian. Thus, by @perko[Chapter 2.9, Theorem 3], we conclude that these points are also _Liapunov stable_, but not _asymptotically stable_.

#figure(
  table(
    columns: 4,
    table.header([*$μ$*], [*$vv(x)^*$*], [*Class*], [*Stability*]),
    $μ < 0$, $(0,0)$, [Center], [Liapunov Stable],
    $μ > 0$, $(0,0)$, [Saddle], [Unstable],
    $μ > 0$, $(sqrt(μ),0)$, [Center], [Liapunov Stable],
    $μ > 0$, $(-sqrt(μ),0)$, [Center], [Liapunov Stable],
  ),
  caption: [
    Classification and stability summary  of @sec:1b
  ]
)

#line(length: 100%)
--- Jeppe ↑ -- Yunis ↓

At the critical point $vv(x^*_0) = (0,0)$, we linearise the system by finding the Jacobian:
$
  D vv(f(vv(x^*_0))) =
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

== Sketching the Phase Portrait & Flow <sec:1c>

Using the results from @sec:1b, we can sketch the phase portraits for both cases $μ > 0$ and $μ < 0$ as seen in @fig:1c_sketch.

The homoclinic orbit shown in purple can be observed encircling the stable centres at $(±sqrt(μ), 0)$ and crossing at the saddle point in $(0, 0)$.

#figure(
  image("../img/1_c_sketch.jpg", height: 6cm),
  caption: [Sketch of the phase plane of vector field $vv(f)$ from @eq:1_duffing_original for $μ > 0$ (left) and $μ < 0$ (right).],
)<fig:1c_sketch>

We additionally plot the numerically computed phase portraits rendered using Python's `matplotlib` library as shown in @fig:1c_numerical, which confirms our qualitative sketches from @fig:1c_sketch.

#figure(
  image("../img/1_c_numerical.png"),
  caption: [Sketch of the phase portrait of vector field $vv(f)$ from @eq:1_duffing_original for $mu > 0$ (left) and $mu < 0$ (right). Critical points are marked in red (saddles) and blue (centers). Arrows indicate the flow of the system.],
)<fig:1c_numerical>

== Analytical Homoclinic Orbit (1st Quadrant) <sec:1d>

Before we begin discussions, it helps to formulate a more intuitive definition of _homoclinic orbits_~@perko[Chapter 3.3]. A homoclinic orbit is a trajectory in the phase space of a dynamical system that starts and ends at the same saddle point. In other words, it is a path that leaves a saddle point along its unstable manifold and then returns to the same saddle point along its stable manifold.
// JK: Not in this case, always, I think. In this case, it does a tour around other critical points in between.

In the phase portrait for $mu > 0$ shown on the right panel of @fig:1c_numerical, we can observe a homoclinic orbit associated with the saddle point $x^*_0$ at the origin $(0,0)$. These orbits leave the saddle point along its unstable manifold and return to it along its stable manifold, forming closed loops around the other critical points $x^*_(1,2)$ at $(plus.minus sqrt(mu), 0)$.
We note that the purple cycle shown in the right panel of @fig:1c_sketch is determined by the equilibrium point at the origin, which in combination with the homoclinic orbit forms a _separatrix cycle_~@perko[Chapter 3.3]. This cycle separates different types of motion in the phase space, distinguishing between trajectories that remain close to the center points and those that diverge away.

Now we can determine a function $h in C^1$, that prescribes these homoclinic orbits in the phase space, where notably the treatment of $h$ may require piecewise definitions to avoid construction of a many-valued function. This function $h: RR arrow RR$ will describe a graph in the phase plane according to
$
x_2 = h(x_1)
$

We require that a solution trajectory lying on the curve formed by $h(x_1)$ must remain on that curve for all time $t$, which requires that the derivatives of $h(x_1)$ and $x_2$ are equal along that curve, yielding the following ordinary differential equation involving $h$:
$
  dot(x)_2 = dv(h(x_1), t) = dv(h, x_1) dot(x)_1
  ⇔
  dv(h, x_1) = dot(x)_2 / dot(x)_1 = (μ x_1 - x_1^3) / x_2
$ <eq:1d_invariance_condition>

Which we may solve by separation of variables:
$
  dv(h, x_1) &= (μ x_1 - x_1^3) / h(x_1)\
  &⇕\
  h(x_1) dif h &= (μ x_1 - x_1^3) dif x_1\
  &⇕\
  h(x_1)^2 &= μ x_1^2 - 1/2 x_1^4 + C^2\
  &⇕\
  h(x_1) &= cancel(±)sqrt(μ x_1^2 - 1/2 x_1^4 + C^2)\
$

Where $C in RR$ is the integration constant and the negative root is discarded on the grounds that $x_2 > 0$ if we restrict ourselves to the first quadrant of the phase plane.

Since we require that the homoclinic orbit passes through the origin, we impose $h(0) = 0 = sqrt(0+0+C^2) = 0$, thus yielding the particular solution:
$
h(x_1) = x_1 sqrt(μ - 1/2 x_1^2) wider x_1 in [0, sqrt(2 μ)]
$ <eq:1d_homoclinic_solution>

Where we have found the domain of $h$ by requiring that the radicand is non-negative to obtain the interval:
$
μ - 1/2 x_1^2 ≥ 0 ⇔ x_1^2 ≤ 2 μ ⇔ x_1 in [0, sqrt(2 μ)]
$

Regarding the map $h(x_1)$ when considered only in the first quadrant of the phase plane, ($x_1 > 0, x_2 > 0$), we note that $h$ is _not_ a one-to-one mapping between $x_1$ and $x_2$, though it is non-negative in this quadrant and as such the we find $0 < x_2 = dot(x)_1$, which implies that the position along $x_1$ is monotonically increasing with time $t$, which could be what the question hints at.
As such, we conclude that both $x_1$ and $(x_1, x_2)$ are one-to-one mappings with respect to time $t$
in this quadrant.

To prove that $h(x_1)$ cannot be one-to-one, we may simply refer to the graph in @fig:1c_numerical and pick a counter example of two points in the first quadrant for which the second coordinate $x_2$ is the same. Alternatively, the derivative of $h(x_1)$ may be investigated to find that both negative and positive slopes exist in the first quadrant, indicating that $h$ is not monotonic and thus not one-to-one.

Note that we can interpret the maximum possible value for $x_1 = sqrt(2 mu)$ as the turning point of the homoclinic orbit in the first quadrant, after which it returns to the origin along the stable manifold along the orbit on the fourth quadrant ($x_1 > 0, x_2 < 0$). Thus it starts at the origin $h(x_1=0)=0$ and transitions into the fourth quadrant at $h(x_1 = sqrt(2 mu))=0$, see @fig:1c_numerical.

// ---


// And as $x_1$ for all $x^*_(0,1,2)$, we plug in $h(0)^2 = 2C = 0$ from @eq:1d_homoclinic_solution_integral to solve for $C=0$. Therefore:

// To derive a first order differential equation for the solution of the graph $x_2 = h(x_1)$, along that graph curve, it must satisfy:

// $
//   frac(d x_2, d x_1) = frac(dot(x)_2, dot(x)_1) = frac(f_2(x_1, x_2), f_1(x_1, x_2)) = frac(mu x_1 - x_1^3, x_2) = frac(mu x_1 - x_1^3, h(x_1)) = frac(d h, d x_1)
// $<eq:1d_homoclinic_graph_condition>

// Hence, we can derive the first order ODE for $h(x_1)$:

// $
//   frac(d h, d x_1) = dot(h)(x_1) = (mu x_1 - x_1^3) / h(x_1)
// $<eq:1d_homoclinic_ode>

// This is homoclinic at the origin, as $h(0) = 0$ and the first quadrant satisfies $h(x_1) > 0$ for $x_1 > 0$.

// @eq:1d_homoclinic_ode is a separable first order ODE for $h(x_1)$, which we can solve by separating variables (multiplying LHS by $h(x_1)$ and RHS by $d x_1$ respectively) and integrating both sides:

// $
//   frac(d h, d x_1) &= (mu x_1 - x_1^3) / h(x_1) \
//   arrow.l.r.double quad quad quad quad
//   h(x_1) d h &= (mu x_1 - x_1^3) d x_1 \
//   arrow.l.r.double quad quad quad
//   ∫ h(x_1) d h &= ∫ (mu x_1 - x_1^3) d x_1 \
//   arrow.l.r.double quad quad quad quad
//   1/2 h(x_1)^2 &= 1/2 mu x_1^2 - 1/4 x_1^4 + C \
//   arrow.l.r.double quad quad quad quad
//   h(x_1)^2 &= mu x_1^2 - 1/2 x_1^4 + 2C
// $<eq:1d_homoclinic_solution_integral>

// Where $C in RR$ is the integration constant.

// Since we require that the homoclinic orbit passes through the origin, we know that $h(0)=0$. And as $x_1$ for all $x^*_(0,1,2)$, we plug in $h(0)^2 = 2C = 0$ from @eq:1d_homoclinic_solution_integral to solve for $C=0$. Therefore:

// $
//   h(x_1)^2 &= mu x_1^2 - 1/2 x_1^4 \
//   h(x_1) &= plus.minus x_1 sqrt(mu - frac(x_1^2, 2)) \
// $<eq:1d_homoclinic_solution>

// For the first quadrant, we take the positive root, yielding the homoclinic orbit function:

// $
//   h(x_1)^+ = x_1 sqrt(mu - frac(x_1^2, 2)) quad quad forall x_1 in [0, sqrt(2 mu)]
// $

// Note that we can interpret the maximum possible value for $x_1 = sqrt(2 mu)$ (for the root to ne non-negative) as the turning point of the homoclinic orbit in the first quadrant, after which it returns to the origin along the stable manifold along the orbit on the fourth quadrant ($x_1 >0, x_2 < 0$). Thus it starts at the origin $h(x_1=0)=0$ and transitions into the fourth quadrant at $h(x_1 = sqrt(2 mu))=0$, see @fig:1c_numerical.

== Plotting the Homoclinic Orbit (1st quadrant) <sec:1e>

#figure(
  image("../img/1_e.png",
  height: 8cm),
  caption: [Plot of the homoclinic orbit $x_2 = h(x_1)$ from @eq:1d_invariance_condition for $mu>0$. The orbit starts and ends at the saddle point at the origin $(0,0)$, forming a closed loop around the other critical points.],
)<fig:1e>

In @fig:1e we observe that alongside the homoclinic orbital path in the fourth quadrant, the the orbit path in the first quadrant encircles a stable center at $(sqrt(mu), 0)$ with each trajectory within the encircled area corresponding to a different energy level of the system, or equivalently a unique bounded orbit around the center. Each of these trajectories are invariant sets of the system. This aligns with our intuition of a conserved energy system and is consistent with @fig:1c_sketch and @fig:1c_numerical.

We note that all trajectories on the orbit are tangent to to it as expected and required by the invariance condition in @eq:1d_invariance_condition.

Lastly, we again note that $h(x_1)$ forms a boundary between different types of motion in the phase space, where trajectories below the epigraph correspond to bounded orbits around the center at $(sqrt(mu), 0)$, while those above it form cycles around both centers.

== Analytical Homoclinic Orbit (4th Quadrant) <sec:1f>

Given the nature of @eq:1d_homoclinic_solution, for the fourth quadrant of the phase plot ($x_1 > 0, x_2 < 0$), it is merely taking the negative root of the homoclinic orbit function, i.e.:

$
  h(x_1)^- = - x_1 sqrt(mu - frac(x_1^2, 2)) quad quad forall x_1 in [0, sqrt(2 mu)]
$

The resulting graph of $h(x_1)$ can be seen in @fig:1e in orange, where it mirrors the first quadrant homoclinic orbit across the $x_1$-axis. All statements concerning turning points and behaviour apply analogously.

#pagebreak()
= Problem 2 — Bifurcation <sec:2>

== Bifurcation parameter & diagram <sec:2a>

A bifurcation parameter is a parameter in a dynamical system that, when varied, causes a change in the system's behavior or structure. Bifurcations often lead to the emergence or disappearance of critical points, periodic orbits, or changes in stability at those points.

As our analysis in @sec:1 showed, the parameter $mu$ in the Duffing equation from @eq:1_duffing_original acts as a bifurcation parameter. Specifically, we observed an emergence of two additional critical points $x^*_(1,2)$ at $(plus.minus sqrt(mu), 0)$ when $mu$ transitions from negative to positive values. This indicates a bifurcation occurring at the critical value $mu=0$.

#figure(
  image("../img/2_a_bifurcation_param_diagram.png",
  height: 8cm),
  caption: [Bifurcation diagram illustrating the emergence of critical points in the Duffing equation from @eq:1_duffing_original as the bifurcation parameter $mu$ varies.],
)<fig:2a>



- comment on the plot?!

== Type of Bifurcation <sec:2b>
The bifurcation observed in @fig:2a from Duffing equation from @eq:1_duffing_original as the parameter $mu$ crosses zero is the prototype of a *supercritical pitchfork bifurcation* (@ REF).
In this supercritical pitchfork bifurcation, as the bifurcation parameter $mu$ passes through a transition value (here, $mu=0$), a single center stable critical point at the origin $(0,0)$ for $mu<0$ becomes an unstable saddle point for $mu>0$, and two new stable critical points (at $(plus.minus sqrt(mu), 0)$) emerge symmetrically on either side of the original point.

== Hyperbolicity for bifurcation transition $mu=0$ <sec:2c>
At $mu=0$ all eigenvalues of the Jacobian at the critical point

$
  x^*_0 arrow D_f(x^*_0) = D_f(0,0) = mat(0, 1; 0, 0)
$

are zero, i.e. $lambda_(1,2) = 0$, it has algebraic multiplicity of two at $0$. This indicates that the critical point is non-hyperbolic at the bifurcation point $mu=0$. Hence, the linearization theorem that requires hyperbolicity does not apply here, and we cannot draw conclusions about the stability of the critical point based solely on the linearized system.

== Lyapunov Stability for non-hyperbolic case <sec:2d>
As we determined in @sec:2c, the critical point $x^*_0$ at the bifurcation point $mu=0$ is non-hyperbolic and we cannot conclude stability from the linearized system alone. Hence, we can use a Lyapunov function to analyse the stability of this critical point. A very common first guess is the simple quadratic form $V(x_1, x_2) = x_1^2 + x_2^2$. Another common choice is the energy function; in this case a natural candidate, since the system can be interpreted as a mechanical oscillator.

$
  V(x_1, x_2) = frac(1,2) x_2^2 + frac(1,4) x_1^4
$<eq:2d_lyapunov_function>

The Lyapunov function needs to be positive definite, i.e. $V(x_1, x_2) > 0$ for all $(x_1, x_2) ≠ (0,0)$ and $V(0,0) = 0$. This is satisfied for @eq:2d_lyapunov_function, as both terms are squares and thus non-negative, and the only point where $V$ is zero is at the origin.

Next, we compute the time derivative of $V$ along the trajectories / solutions $x_1(t), x_2(t)$ of the system $bold(f)$ for $mu=0$:

$
  dot(V)(x_1(t), x_2(t)) &= nabla_t V(x_1(t), x_2(t)) (dot(x)_1(t), dot(x)_2(t)) quad quad "via chain rule"\
  &= partial V_(x_1) dot(x)_1 + partial V_(x_2) dot(x)_2
  quad quad quad quad quad quad quad quad
  "inserting into function"\
  &= x_1^3 dot(x)_1 + x_2 dot(x)_2 \
  &= x_1^3 x_2 + x_2 (mu - x_1^3) \
  &= x_1^3 x_2 + x_2 (- x_1^3) \
  &equiv 0
$

Since $V$ is positive definite and $dot(V)≡0$, we can interpret that the energy in the system stays constant along each trajectory. That means, trajectories are contained within their (energy-)level-set of the Lyapunov function.
Thus, the origin $(x_1, x_2) = (0,0)$ at $mu=0$ is Lyapunov (stable): small initial conditions remain small for all time.
However, because $dot(V)$ is identically zero (not negative definite), one cannot conclude asymptotic stability. In fact, the origin is not asymptotically stable (solutions near origin are periodic orbits—center behavior, similar to the linear centers $x^*_(1,2)$ with purely imaginary Eigenvalues in the linearized system).

So, the origin is a stable center (Lyapunov stable) but not asymptotically stable.


== Contour Plot of Lyapunov Function <sec:2e>

#figure(
  image("../img/2_e_contourV_phase.png",
  height: 8cm),
  caption: [Contour plot of the Lyapunov function $V(x_1, x_2)$ from @eq:2d_lyapunov_function for $mu=0$. The contours represent level sets of constant energy in the system. Overlay over the phase plot with vector field arrows.],
)<fig:2e>

As clearly visible in @fig:2e, the contours of the Lyapunov function $V(x_1, x_2)$ from @eq:2d_lyapunov_function represent not only level sets of constant energy but also trajectories as solutions of the system @eq:1_duffing_original for $mu=0$. This is consistent with our earlier analysis in @sec:2d, where we found that the time derivative of the Lyapunov function $dot(V)$ is identically zero along trajectories, indicating that energy remains constant.

*REVIEW*
Further we observe, that the level sets along with trajectories of the vector field form elliptical closed orbits around the origin. The elliptical shape could be attributed to the stretch of the system along the $x_1$-axis, due to the quartic term in $V$ (@eq:2d_lyapunov_function), compared to only the quadratic term for $x_2$.


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

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
  *02417 Times Series Analysis - Assignment 4*
])

#strong[Authors]:
- Jeppe Klitgaard `<s250250@dtu.dk>`
- Yunis Wirkus `<s250700@dtu.dk>`

Date: 2025-10-08

#pagebreak()


//! Content
= Problem 1 — Rescaling the Duffing Equation <sec:1>

We are given the Duffing equation:
$
  M dv(X, T, 2) + α dv(X, T) = k X - κ X^3 + A sin(Ω T)
$ <eq:1_duffing_original>

Into which we introduce the following linear transformations of the coordinates $X$ and $T$:
$
  T = k_t t wider X(T) = k_x x(t)
$ <eq:1_substitutions>

We consider the variable transformation of the first and second derivatives using the chain rule:
$
  dv(X, T) &= dv(X, x) dv(x, t) dv(t, T) = k_x/k_t dv(x, t)\
  dv(X, T, 2) &= dv(, T) [dv(X, T)] = dv(, t) dv(t, T)[k_x/k_t dv(x, t)] = k_x/k_t^2 dv(x, t, 2)
$ <eq:1_derivative_transformations>

Plugging @eq:1_substitutions and @eq:1_derivative_transformations into @eq:1_duffing_original, we obtain:

$
  M k_x/k_t^2 dv(x, t, 2) + α k_x/k_t dv(x, t) = k k_x x - κ k_x^3 x^3 + A sin(Ω k_t t)
$ <eq:1_duffing_substituted_1>

== Finding the Scaling Factors <sec:1a>
We seek to find rescaling factors $k_x, k_t$ such that we obtain the following form, correponding to normalisation of the force term $f(X) = k X - κ X^3$:
$
  dv(x, t, 2) + β dv(x, t) = x - x^3 + a sin(ω t)
$ <eq:1_duffing_target>

We impose the condition that the coefficient on $dv(x, t, 2)$ must be $1$ by dividing @eq:1_duffing_substituted_1 through by $M k_x/k_t^2$:

$
  dv(x, t, 2) + (α k_t)/M dv(x, t) = (k k_t^2)/M x - (κ k_t^2 k_x^2)/M x^3 + (A k_t^2)/(M k_x) sin(Ω k_t t)
$ <eq:1_duffing_substituted_2>

Our remaining normalisation conditions imply:
$
  (k k_t^2)/M &= 1 quad &&<=> quad k_t &&= ±sqrt(M/k)\
  (κ k^2_x k_t^2)/M &= 1 quad &&<=> quad k_x &&= ±sqrt(M/(κ k_t^2))
$ <eq:1_conditions_1>

Where we shall take the positive roots for both $k_t$ and $k_x$ for reasons discussed below in @sec:1b,
such that:
$
  k_t = sqrt(M/k) wider and wider k_x = sqrt(k/κ)
$ <eq:1_scaling_factors>

Where we note that $k_x$ appears to parameterise the relative strength of the linear and non-linear components of the force term $f(X)$, while $k_t$ appears to parameterise the relative strength of the inertial term and force term $f(X)$.

Plugging in our scaling factors from @eq:1_scaling_factors into @eq:1_duffing_substituted_2, we obtain:
$
  dv(x, t, 2) + α/sqrt(k M) dv(x, t) = x - x^3 + A sqrt(κ/k^3) sin(Ω sqrt(M/k) t)
$

Which when compared against @eq:1_duffing_target reveals the rescales parameters as:
$
β &= α/sqrt(k M) wider a = A sqrt(κ/k^3) wider ω = Ω sqrt(M/k)
$

== Choice of Scaling Factors <sec:1b>

There is some ambiguity in the question regarding choice of scaling factors $k_x, k_t$.

Firstly, one could have chosen a different _target equation_ than @eq:1_duffing_target, for examply by choosing to normalise the driving term $A sin(Ω T)$ instead of the restoring force term $k X - κ X^3$.

This choice of scaling may benefit subsequent analysis of the system in some cases.

Alternatively, the question may refer to the ambiguity that arose in the choice of sign for $k_x, k_t$ in @eq:1_conditions_1. We consider that the sign choice does not affect those terms where the scaling factors appear squared, but does affect the linear dampening term $α κ_t / M dv(x, t)$ and the driving term $(A k_t^2)/(M k_x) sin(Ω k_t t)$.

We note that a choice of $k_t < 0$ would lead to a sign change in the dampening term, which would imply
that the system is being driven rather than dampened. This is would be considered unphysical for most systems under consideration, so we may discard this choice.

Similarly, a choice of $k_x < 0$ would lead to a sign change in the driving term, which would imply a phase shift of $π$ in the driving force. While this is not unphysical, it is a trivial phase shift that does not affect the dynamics of the system. For this reason, it is reasonable to discard this choice as well.

= Problem 2 — Dynamical Analysis on Scalar First Order ODE

We are given the following first order ODE:
$
dv(x, t) = f(x) = (x+1)^2(1-x)
$ <eq:2_ode>

Where $x : ℝ → ℝ$.

== Sketch and Stationary Points

By inspection $f(x)$ has a double root at $x_0 = -1$ and a single root at $x_0 = 1$,
which are thus the stationary points of the system where $f(x) = 0$.

Instead of computing the derivative of $f(x)$, we compute the sign of $f(x)$ in the intervals defined by the stationary points:
$
  f(-2) &= +3\
  f(-0) &= +1\
  f(2) &= -9\
$

We sketch the function $f(x)$ and its stationary points, as shown below in @fig:2a

#figure(
  image("diagram_2a.jpg", height: 6cm),
  caption: [Sketch of $f(x)$ with stationary points at $x∈ {-1, 1}$.],
) <fig:2a>

== Solution Flow and Stability

We compute the derivative $f'(x_0)$:
$
  f'(x) = -(x+1)(3x-1)
$

Which we evaluate at the stationary points:
$
f'(-1) = 0 wider and wider f'(1) = -4
$

We consider the single root at $x_0 = 1$, which has a negative derivative corresponding to a stable node.

Upon inspecting the stationary point at $x_0 =-1$, we find a vanishing derivative as expected from the double root, meaning we cannot determine stability from the linear approximation. Instead, we must consider the sign of the second derivative, $f''(x)$, in order to determine stability for this point:

$
  f''(x) &= -6x -2\
  f''(-1) &= +4
$

To aide our understanding, we start at Ansatz $x(t) = x_0 + η(t)$ with temporal derivative $dv(x, t) =f(x)= dv(η, t)$ and perform a Taylor expansion of $f(x)$ around $x_0$, noting that we have already found $f(x_0) = f'(x_0) = 0$:
$
  f(x) &= cancel(f(x_0)) + cancel(f'(x_0) η) + 1/2 f''(x_0) η^2 + O(η^3)
$

Thus we find:
$
  dv(η, t) &≈ 1/2 f''(x_0) η^2\
  eval(dv(η, t))_(x_0=-1) &≈ 2 η^2\
$

We may solve this ODE using separation of variables:
$
  (dif η) /η^2 &= 2 dif t\
  ∫ 1/η^2 dif η &= ∫ 2 dif t\
  1/η &= -2 t + C\
  η(t) &= 1/(-2 t + C)
$

We solve for the integration constant $C$ using initial condition $η(0) = η_0$:
$
  η(0) = η_0 = 1/(-2(0) + C) = 1/C wider => wider C = 1/η_0
$

Thus obtaining:
$
  η(t) = 1/(1/η_0 - 2t) = η_0 / (1 - 2 η_0 t)
$

This allows us to consider the stability of the point $x_0 = -1$ for time $t > 0$ based on the
sign of the initial perturbation $η_0$:
$
  η_0 > 0 quad&⇒quad η(t) → +∞ &&quad "as" t → 1/(2 η_0)_- quad &&⇒ quad "unstable"\
  η_0 < 0 quad&⇒quad η(t) → 0_- &&quad "as" t → ∞ quad&&⇒quad"stable from below"
$

Thus we have _half-stability_ at the point $x_0 = -1$, where the system is stable for perturbations from the left, but unstable for perturbations from the right.

We sketch the solution flow of the system in a phase line diagram using the information we have obtained
as shown in @fig:2b.

#figure(
  image("diagram_2b.jpg", height: 4cm),
  caption: [Solution flow of the system along the $x$-axis, showing stability of the stationary points.],
) <fig:2b>

As noted above, it is not possible to determine the stability of the point $x_0 = -1$ from the linear approximation alone, as $f'(-1) = 0$.

== Sketching Solution Trajectories <sec:2c>

We sketch the solutions $ϕ_t (x_0)$ as a function of time $t > 0$ as shown in @fig:2c.

#figure(
  image("diagram_2c.jpg", height: 9cm),
  caption: [Sketch of solution trajectories $ϕ_t(x_0)$ for various initial conditions $x_0$.],
) <fig:2c>

Where the half-stability of the point $x_0 = -1$ is evident from the trajectories approaching the point from below, but diverging away from it for initial conditions above the point.

== Interval Solution <sec:2d>

We will find a solution to @eq:2_ode using separation of variables and decomposition of rational functions.

First, we separate the variables to obtain the following integrals:
$
  ∫ 1/((x+1)^2(1-x)) dif x = ∫ dif t
$

After which we decompose the spatial integrand, $1/f(x)$, into partial fractions:
$
  1/((x+1)^2(1-x)) = A/(x+1) + B/(x+1)^2 + C/(1-x)
$

We clear the denominators by multiplying through by $(x+1)^2(1-x)$ and collect in powers of $x$:
$
  1
  &= A(x+1)(1-x) + B(1-x) + C(x+1)^2\
  &= (C-A) x^2 + (2C-B) x + (A + B + C)\
$

From which we obtain the following system of equations:
$
  C - A &= 0\
  2C - B &= 0\
  A + B + C &= 1\
$

Which we solve as:
$
  A = 1/4 wider B = 1/2 wider C = 1/4
$

We thus have the following integral to solve:
$
  ∫ 1/(4(x+1)) + 1/(2(x+1)^2) + 1/(4(1-x)) dif x = ∫ dif t
$

We do so by integrating the following terms and summing with the appropriate coefficients:
$
  ∫ 1/(x+1) dif x &= ln|x+1| + C_1\
  ∫ 1/(x+1)^2 dif x &= -1/(x+1) + C_2\
  ∫ 1/(1-x) dif x &= -ln|1-x| + C_3\
  ∫ dif t &= t + C_4\
$

Yielding solution:
$
  1/4 ln((1+x)/(1-x)) - 1/2 1/(1+x) = t - t_0
$ <eq:2_solution_interval>

Where we have combined the integration constants into a single constant $t_0$ on the right hand side
and exploited the fact that the arguments of the logarithms are positive in the interval $x ∈ ]-1; 1[$
to drop the absolute values.

== Limits of Interval Solution <sec:2e>

We consider the solution given in @eq:2_solution_interval in the limits $x → -1_+$ and $x → 1_-$.

We first consider the limit $x → -1_+$, in which case the argument of the logarithm approaches zero from above, $0_+$. We know that $ln(0_+) → -∞$ and conclude that the logarithmic term goes to $-∞$ in this limit.
Similarly, by inspection, $lim_(x → -1_+) (1/(1+x)) = +∞$ such that $-1/2 1/(1+x) → -∞$.
Thus, as both terms approach the same limit, we find that $t→-∞$ as $x → -1_+$.

Secondly, we consider the limit $x → 1_-$, in which case the argument of the logarithm approaches $+∞$ such that $lim_(x → 1_-) ln((1+x)/(1-x)) = +∞$. The rational term evaluates to $-1/2 1/(1+1)=-1/4$ in the limit
and is thus trivially dominated by the logarithmic term. As such, we conclude that $t→ +∞$ as $x → 1_-$.

= Problem 3 — Dynamical Analysis on a Two-Dimensional System

We are given a system of differential equations:
$
  dv(vv(x), t) = mm(B) vv(x) wider "where" wider mm(B) = mat(-1, 2; 0, -1), quad vv(x)(0) = vv(x)_0
$ <eq:3_system>

With $vv(x) : ℝ → ℝ^2$ and $vv(x)_0 = (x_01, x_02)^TT$.

== Solution using Matrix Exponentials <sec:3a>

We invoke the Fundamental Theorem of Linear Systems as given by @perko[Chapter 1.4],
which states that a system of the form @eq:3_system has a unique solution given by:
$
  vv(x)(t) = e^(mm(B) t) vv(x)_0
$

As such, our problem has been reduced to finding the matrix exponential $e^(mm(B) t)$.

=== Solution using Coordinate Transformation <sec:3a_transform>
While this problem is easily solved using @perko[Chapter 1.4, Corollary 4],
we will instead opt for the idiomatic approach using a change of basis to cast the problem into its Jordan normal form. We do this by first determining the eigenvalues of $mm(B)$ to be $λ_1 = λ_2 = -1$ using the known result that the eigenvalues of an upper triangular matrix are given by its diagonal elements.

We note that the eigenvalues are equal and that we thus have an _algebraic multiplicity_ of 2, but
we only have a _geometric multiplicity_ of 1. As such, we need to find a _generalized eigenvector_ in addition to the _regular eigenvector_.

The eigenvectors, $vv(v)$ are found by solving the familiar equation:
$
  (mm(B) - λ mm(I)) vv(v) = 0
$

Which trivially yields the eigenvector:
$
vv(v)_1 = (1, 0)^TT
$

As dictated by @perko[Chapter 1.7], we may find the generalized eigenvector $vv(v)_2$ which satisfies:
$
  (mm(B) - λ mm(I))^k vv(v) = 0
$

Where $k$ is the _degree of nilpotency_ of the matrix $(mm(B) - λ mm(I))$, which we find as:
$
  mat(0, 2; 0, 0)^2 = mat(0, 0; 0, 0) wider ⇒ wider k = 2
$

We obtain $vv(v)_2$ using the Jordan chain conditions, which can be used to recursively find generalized eigenvectors @LectureNotesWeek3:

$
  (mm(B) - λ mm(I)) vv(v)_k = vv(v)_(k-1)
$

Thus, we find $vv(v)_2$ by solving:
$
  (mm(B) - λ mm(I)) vv(v)_2 &= vv(v)_1\
  mat(0, 2; 0, 0) vv(v)_2 &= vec(1, 0)\
$

From which we find a general solution, from which we may pick $c=0$ to obtain a particular solution:
$
  vv(v)_2 = (c, 1\/2)^TT = (0, 1\/2)^TT
$

As explained in @perko[Chapter 1.6], we may now construct an invertible matrix $mm(P)$ using our eigenvectors as columns:
$
  mm(P) = mat(vv(v)_1, vv(v)_2) = mat(1, 0; 0, 1\/2)
$

Which transforms $mm(B)$ into its Jordan normal form $mm(C)$:
$
  mm(C) = mm(P)^(-1) mm(B) mm(P) = mat(1, 0; 0, 2) mat(-1, 2; 0, -1) mat(1, 0; 0, 1\/2) = mat(-1, 1; 0, -1)
$

Which we confirm as a Jordan normal form matching _Case II_ of @perko[Chapter 1.5], allowing us to characterise the system as an _improper stable node_ given $λ_1 = λ_2 < 0$.

We may now leverage the properties described in @perko[Chapter 1.3] to find:
$
  e^(mm(B) t) = mm(P) e^(mm(C) t) mm(P)^(-1)
$

Where in the case of $mm(C) = mat(λ, 1; 0, λ)$ we have the known result that:
$
  e^(mm(C) t) = e^(λ t) mat(1, t; 0, 1)
$

Which yield the solution:
$
  e^(mm(B) t) = mm(P) e^(mm(C) t) mm(P)^(-1) = mat(1, 0; 0, 1\/2) e^(-t) mat(1, t; 0, 1) mat(1, 0; 0, 2) = e^(-t) mat(1, 2t; 0, 1)
$

Such that the solution to @eq:3_system:
$
  vv(x)(t) = e^(-t) mat(1, 2t; 0, 1) vv(x)_0
$

=== Alternative Approach <sec:3a_alt>

The treatment above was rather cumbersome, but was undertaken to practice the method of transforming to Jordan normal form and back.

Given the rather simple form of $mm(B)$, we may have instead opted for a more direct approach to finding the matrix exponential.

As before, we will utilise the Fundamental Theorem of Linear Systems as given by @perko[Chapter 1.4],
such that the problem has been reduced to finding the matrix exponential $e^(mm(B) t)$.

We may decompose $mm(e^(mm(B) t))$ as follows:
$
  mm(e^(mm(B) t)) = exp(mat(-t, 2t; 0, -t))
  &= exp(mat(-t, 0; 0, -t) + mat(0, 2t; 0, 0))\
  &= exp(-t mm(I) + mat(0, 2t; 0, 0))\
$

We seek to prove that $-t mm(I)$ and $mat(0, 2t; 0, 0)$ commute, which is trivial upon noting that $t mm(I)$ is a _scalar matrix_ and thus commutes with all matrices of the same dimension. This enables the decomposition of the exponential of a sum into a product of exponentials as per @perko[Chapter 1.3, Proposition 2]:
$
  mm(e^(mm(B) t))
  &= exp(-t mm(I) + mat(0, 2t; 0, 0))\
  &= e^(-t mm(I)) exp(mat(0, 2t; 0, 0))\
  &= e^(-t) mm(I) exp(mat(0, 2t; 0, 0))\
$

We can now evaluate $exp(mat(0, 2t; 0, 0))$ using the definition of the matrix exponential as a power series:
$
  e^(mm(S)) = lim_(n → ∞) sum_(k=0)^n mm(S)^k/k!
$

Where we note that $mm(S) = mat(0, 2t; 0, 0)$ is _nilpotent_ of degree 2 as $mm(S)^2 = mm(0)$, such that the power series truncates after the first order term:
$
  exp(mat(0, 2t; 0, 0)) = mm(I) + mat(0, 2t; 0, 0) + mm(0) = mat(1, 2t; 0, 1)
$

Combining these results, we find:
$
  mm(e^(mm(B) t)) = e^(-t) mm(I) mat(1, 2t; 0, 1) = e^(-t) mat(1, 2t; 0, 1)
$

Which is consistent with the treatment in @sec:3a_transform, yielding the same solution to @eq:3_system:
$
  vv(x)(t) = e^(-t) mat(1, 2t; 0, 1) vv(x)_0
$

== Phase Portrait and Stability <sec:3b>

Our solutions to @sec:3a can aid us in sketching the phase portrait of the system.

The degenerate eigenvalue $λ_1 = λ_2 = -1 < 0$ indicates that the system is _stable_ at the stationary point $vv(x) = (0, 0)^TT$.

As such, for all points (except the trivial case of the origin), we expect the trajectories to point towards the origin as $t → ∞$. From @perko[Chapter 1.5, Case II], we know that the trajectories will be tangent to the eigenvector $vv(v)_1 = (1, 0)^TT$ as they approach the origin.

Lastly, we consider the coordinate transformation matrix $mm(P) = mat(1, 0; 0, 1/2)$ from @sec:3a_transform, which indicates that the $y$-coordinate is scaled by a factor of $1/2$ in the eigenbasis. This means that the trajectories will be compressed in the $y$-direction relative to the $x$-direction,
which will be largely inconsequential for the sketch of the phase portrait.

With this, we are able to sketch the phase portrait of the system as shown in @fig:3b.

// Note: This can easily be validated using https://choosedews.github.io/PhasePlane/
#figure(
  image("diagram_3b.jpg", height: 8cm),
  caption: [Phase portrait of the system showing trajectories approaching the stable, improper node at the origin.],
) <fig:3b>

== Invariant Sets and Subspaces <sec:3c>

We understand an invariant set to be a set such that if the system starts in the set, it remains in the set for all time $t > 0$. We recall from @perko[Chapter 1.9] that we may write the phase space as
the outer product of the stable, unstable, and center subspaces given by $E^s$, $E^u$, and $E^c$ respectively:
$
  ℝ^2 = E^s ⊕ E^u ⊕ E^c
$

Where the subspaces may be found using the spans of the eigenvectors corresponding to the eigenvalues with negative, positive, and zero real parts respectively, for $i ∈ {1, 2}$:
$
  E^s &= Span{vv(v)_i | Re(λ_i) < 0}\
  E^u &= Span{vv(v)_i | Re(λ_i) > 0}\
  E^c &= Span{vv(v)_i | Re(λ_i) = 0}\
$

As such, we find using the eigenvectors $λ_1 = λ_2 = -1 < 0$:
$
  E^s = ℝ^2 wider E^u = {0} wider E^c = {0}
$

We may subdivide the stable subspace, which in this case is the entire phase space, smaller invariant sets as well. Trivially, the stationary point at the origin is an invariant set:
$
  {(0, 0)^TT} ⊂ E^s
$

Additionally, as illustrated by the sketch in @fig:3b and by inspection of @eq:3_system, we note that
the first axis itself is an invariant set, as any point on the axis will remain on the axis for all time $t > 0$:
$
  {(x_1, 0)^TT | x_1 ∈ ℝ} ⊂ E^s
$

Trivially, the entire phase space is of course also an invariant set.

In conclusion, we find the following invariant sets:
$
  {(0, 0)^TT} ⊂ {(x_1, 0)^TT | x_1 ∈ ℝ} ⊂ ℝ^2
$

Which are all stable invariant sets.

#bibliography("refs.bib")

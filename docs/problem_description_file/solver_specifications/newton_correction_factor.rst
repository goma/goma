****************************
Newton Correction Factor
****************************

::

	Newton correction factor = <float_list>

-----------------------
Description / Usage
-----------------------

This required card indicates the damping (or relaxation) factor for the Newton updates
and offers customization of the relaxation choice based on the size of the nonlinear
residual from the Newton iteration. Definitions of the <float_list> input parameters,
from one (f\ :sub:`1`) to six (f\ :sub:`2`, ... f\ :sub:`6`) floating point numbers (one required and five optional),
are as follows:

<float1>
    f\ :sub:`1`, damping factor for the Newton updates, where ( 0.0 < f\
    :sub:`1` ≤ 1.0 ). A value of 1.0 gives the usual Newton’s method,
    otherwise, only a portion of the Newton update is applied to the solution.
    Values near 0 (e.g., 0.1) may be used effectively to aid convergence for
    sensitive problems where the initial guess is not very close to the final
    solution for the first several Newton iterations. This parameter can also
    be controlled from the command line (see **-r** option, Command-line
    Arguments, Chapter 3).
[floatn]
    These five floats [f\ :sub:`2`, ... f\ :sub:`6`] are optional but give
    a way to more finely control the amount of relaxation applied to Newton
    updates. See the description below and the example for an explanation.

------------
Examples
------------

A simple example is the following:
::

	Newton correction factor = 0.1

This tells *Goma* to take the specified number of Newton iterations (from the *Number of
Newton Iterations* card) at a fixed relaxation parameter of 0.1. This is a moderately
large amount of relaxation, but of course “moderately large” is always problem
dependent.

A more interesting example:
::

	Newton correction factor = 0.8 1.0e-6 0.4 1.0e-4 0.1 1.0e-3

causes the following relaxation scheme to be used according to the L\ :sub:`∞` norm of the
nonlinear residual:

* If L\ :sub:`∞` > 1.0e – 3, the relaxation factor is taken as 0.1.

* If 1.0e – 4 < L\ :sub:`∞` ≤ 1.0e – 3, the relaxation factor is taken as 0.4.

* If 1.0e – 6 < L\ :sub:`∞` ≤ 1.0e – 4, the relaxation factor is taken as 0.8.

* If L\ :sub:`∞` ≤ 1.0e – 6, the relaxation factor is taken as the usual Newton’s method
  relaxation of 1.0.

The default relaxation level for small residuals is 1.0.

-------------------------
Technical Discussion
-------------------------

The relaxation factor is used to intentionally shorten the solution update vector
computed by the Newton iteration. There are many factors that can cause the effective
radius of convergence of Newton’s method to be quite small or malformed:

* the underlying nonlinear problem is stiff,

* the initial solution is poor,

* non-analytic constitutive models or boundary conditions,

* poor linear solver performance, etc.

Under these kinds of circumstances, the update computed by Newton’s method may be
too large and end up not improving the overall solution. In such cases it is
recommended that one uses some relaxation (e.g., 0.9), and possibly a lot (e.g., 0.05).

What one really wishes to do is to use shortened updates when far from convergence,
and full updates as the solution converges. This is the capability that the optional five
parameters makes available. While they don’t directly measure how far the solution is
from convergence, it does use the residual as an indicator. The full set of six parameters
allows the user to specify four different residual intervals with four different relaxation
factors. The f\ :sub:`1`, f\ :sub:`3` and f\ :sub:`5` values are relaxation factors and must lie in 0.0 < f\ :sub:`i` ≤ 1.0,
while the f\ :sub:`2`, f\ :sub:`4`, and f\ :sub:`6` values are interval endpoints. The supplied interval endpoints
must be in ascending order, 0 < f\ :sub:`2` < f\ :sub:`4` < f\ :sub:`6`. Although no such restriction is put on the
relaxation factors, they should generally satisfy 0 < f\ :sub:`5` ≤ f\ :sub:`3` ≤ f\ :sub:`1` ≤ 1.0.


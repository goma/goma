***********************************
Normalized Correction Tolerance
***********************************

::

	Normalized Correction Tolerance = <float>

-----------------------
Description / Usage
-----------------------

This optional card sets the tolerance for a mixed measure of the size of the update
vector which must be satisfied for the solution to be considered converged. The input
parameter is defined as

<float>
    **rel,** a floating point value ( rel ≥ 0.0 ) used as the convergence
    tolerance for the mixed measure of the update vector (defined in the
    Technical Discussion).

When the *Normalized Correction Tolerance* card is omitted, the default value of **rel** is
1.0e+10.

------------
Examples
------------

Following is a sample card:
::

	Normalized Correction Tolerance = 1.0e-4

-------------------------
Technical Discussion
-------------------------

The mixed measure used here is:

.. math::

   \sqrt{ \sum \frac{ {\Delta x_i}^2}{1 + {x_i}^2} }

This measures the relative size of the update vector when the solution vector is large
(i.e., size of unknowns is greater than 1), and measures the absolute size of the update
vector when the solution vector is small (i.e., size of unknowns is much less than 1).

This mixed measure must be less than **rel,** in addition to the nonlinear residual
satisfying the absolute residual tolerance specified in the *Normalized Residual
Tolerance* card for a solution to be considered converged.

If rel < 1.0 (larger values are not really imposing any restrictions), mixed measure
values are output instead of the update vector norms.


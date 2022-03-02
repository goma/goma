*********************************
Normalized Residual Tolerance
*********************************

::

	Normalized Residual Tolerance = <float>

-----------------------
Description / Usage
-----------------------

This is the per matrix version of this card. This overrides Normalized Residual Tolerance for each matrix

This required card indicates the value of the L\ :sub:`2` norm of the global nonlinear residual
vector that indicates termination of Newton’s method (i.e., convergence). The input
parameter is defined as

<float>
    **tol,** a non-negative floating point number ( tol ≥ 0.0 ) specifying the
    L\ :sub:`2` convergence tolerance for the global nonlinear residual vector.

The *Normalized Residual Tolerance* card is required; there is no default.

------------
Examples
------------

Following is a sample card:
::

   MATRIX = 1
       Normalized Residual Tolerance = 1.0e-8
       Number of EQ = 3

-------------------------
Technical Discussion
-------------------------

Newton’s method is terminated when the global nonlinear residual falls below tol, or
the maximum number of iterations specified in the *Number of Newton Iterations* is
reached.


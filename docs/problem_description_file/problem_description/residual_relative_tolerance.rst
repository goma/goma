*********************************
Residual Relative Tolerance
*********************************

::

	Residual Relative Tolerance = <float>

-----------------------
Description / Usage
-----------------------

This is the per matrix version of this card. This overrides Residual Relative Tolerance for each matrix


<float>
    **tol,** a non-negative floating point number ( tol ≥ 0.0 ) specifying the
    L\ :sub:`2` ratio convergence tolerance for the global nonlinear residual vector.

The *Residual Relative Tolerance* card is not required; the default is 1e10 effectively turning it off.

------------
Examples
------------

Following is a sample card:
::

   MATRIX = 1
       Residual Relative Tolerance      = 1e-6
       Number of EQ = 3

-------------------------
Technical Discussion
-------------------------



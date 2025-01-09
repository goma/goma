*********************************
Residual Relative Tolerance
*********************************

::

	Residual Relative Tolerance = <float>

-----------------------
Description / Usage
-----------------------

Residual Relative Tolerance requires the convergence to reach a relative
tolerance compared to the initial normalized residual.

For example if your initial normalized residual is :math:`r0=1e2` and you set the relative tolerance to :math:`rt=1e-6`
then you must reach :math:`rf <= r0*rt`


<float>
    **tol,** a non-negative floating point number ( tol ≥ 0.0 ) specifying the
    L\ :sub:`2` ratio convergence tolerance for the global nonlinear residual vector.

The *Residual Relative Tolerance* card is not required; the default is 1e10 effectively turning it off.

------------
Examples
------------

Following is a sample card:
::

   Residual Relative Tolerance      = 1e-6

-------------------------
Technical Discussion
-------------------------



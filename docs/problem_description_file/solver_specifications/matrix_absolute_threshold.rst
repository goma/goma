*****************************
Matrix Absolute Threshold
*****************************

::

	Matrix Absolute Threshold = <float>

-----------------------
Description / Usage
-----------------------

This card is only available with the Trilinos library. It allows the user to specify a lower
bound for either a diagonal entry or a singular value. The exact meaning depends on the
kind of preconditioner used (scalar-based or block-based). The legal values are:

<float>
    **t,** a floating point number ( t ≥ 0.0 ) that specifies a minimum
    threshold value for diagonal or singular value.

Along with the *Matrix Relative Threshold* card, this card gives the user the ability to
modify what matrix the preconditioner operates on. See the *Matrix Relative Threshold*
card for a full description.

If this card is omitted, the default is 0.0.

------------
Examples
------------

A sample input card follows:
::

	Matrix Absolute Threshold = 1.e-4

-------------------------
**Technical Discussion**
-------------------------

Refer to the discussion for card *Matrix Relative Threshold*. The appropriate values for
the threshold can vary over many orders of magnitude depending on the situation.
Refer to Schunk, et. al., 2002 for information and for further guidance.



--------------
References
--------------

SAND2001-3512J: Iterative Solvers and Preconditioners for Fully-coupled Finite
Element Formulations of Incompressible Fluid Mechanics and Related Transport
Problems, P. R. Schunk, M. A. Heroux, R. R. Rao, T. A. Baer, S. R. Subia and A. C.
Sun, March 2002.

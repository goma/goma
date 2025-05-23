*****************************
Matrix Relative Threshold
*****************************

::

	Matrix Relative Threshold = <float>

-----------------------
Description / Usage
-----------------------

This card is only available with the Trilinos library. The effect of this card is to impose
a relative lower bound to either a diagonal value or a singular value. The legal values
for <float> are:

<float>
    **r**, a floating point number ( r ≥ 0.0 ) that specifies a relative
    threshold.

If this card is omitted, the default is 0.0.

------------
Examples
------------

A sample input card follows:
::

	Matrix Relative Threshold = 1.e-4

-------------------------
Technical Discussion
-------------------------

This card, along with the *Matrix Absolute Threshold* card, allow the user to modify the
linear system prior to calculation of the preconditioner. Note that the modification is
only to change the “initial condition” of the preconditioner--it does not actually change
the linear system.

Let **t** be the value specified with the Matrix Absolute Threshold card. For a scalar-based
preconditioner (**ilut, ilu, rilu, icc**), each value on the diagonal undergoes the following
substitution:

 .. math::

    d_{\mathrm{new}} = r * d_{\mathrm{old}} + \mathrm{sgn} \left( d_{\mathrm{old}}  \right) * t

For the **bilu** preconditioner, each singular value of the diagonal block preconditioner is
compared to:

 .. math::

    \sigma_{\mathrm{min}} = r * \sigma_1 + t

where σ1 is the largest singular value of the diagonal block under consideration. All :math:`\sigma_k`
are modified (if necessary) to be at least as large as :math:`\sigma_{\mathrm{min}}`.

The appropriate values for the threshold can vary over many orders of magnitude
depending on the situation. Refer to Schunk, et. al., 2002 for information and for
further guidance.

--------------
References
--------------

SAND2001-3512J: Iterative Solvers and Preconditioners for Fully-coupled Finite
Element Formulations of Incompressible Fluid Mechanics and Related Transport
Problems, P. R. Schunk, M. A. Heroux, R. R. Rao, T. A. Baer, S. R. Subia and A. C.
Sun, March 2002.

***************************
Matrix Polynomial Order
***************************

::

	Matrix polynomial order = <integer>

-----------------------
Description / Usage
-----------------------

This optional card allows selection of polynomial order when a polynomial
preconditioning option is selected (see the *Preconditioner* card). The input parameter is
defined as:

<integer>
    Number of steps, **k** ( ≥ 0 ), to take when using matrix polynomial based
    preconditioners (Jacobi and symmetric Gauss-Seidel, for example).

If the Matrix polynomial order card is omitted, then the default selection is **k=3.**

------------
Examples
------------

Following is a sample card:
::

	Matrix polynomial order = 4

-------------------------
Technical Discussion
-------------------------

When used, the value of this parameter should be greater than 0, and probably no more
than 10. In some, if not all, cases, a value of 0 is meaningless.

This card is not used if the preconditioner does not use matrix polynomials.




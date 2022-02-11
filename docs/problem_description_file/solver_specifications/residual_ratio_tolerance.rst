****************************
Residual Ratio Tolerance
****************************

::

	Residual Ratio Tolerance = <float>

-----------------------
Description / Usage
-----------------------

This optional card sets the convergence criterion for the iterative solution of the linear
matrix system solved at each Newton iteration. The input parameter is defined as

<float>
    **tol,** a non-negative real number ( tol ≥ 0.0 ) specifying the value of
    the convergence criterion.

The default value of **tol** is 1.0e-6.

------------
Examples
------------

Following is a sample card:
::

	Residual Ratio Tolerance = 1.0e-3

-------------------------
Technical Discussion
-------------------------

The value of **tol** is ignored when a direct factorization algorithm (such as **lu**) for 
the
linear solve is specified in the *Solution Algorithm* card. When an iterative matrix
solution technique is specified (such as **gmres**), **tol** acts as the inner iteration
termination relative tolerance. Letting r\ :sub:`0` represent the initial residual norm, 
when the
n\ :sub:`th` iteration’s linear residual norm r\ :sub:`n` satisfies 					
r\ :sub:`n` / r\ :sub :`0` ≤ tol, the iterative 
solution is
deemed acceptable and the inner iterations terminate. The number of iterations required
is reported under the LIS column of the Newton iteration output. If the maximum
number of iterations (specified in the *Maximum Linear Solve Iterations* card) is
reached, then **max** appears instead of a number. Although the standard residual is
usually used as the residual norm, the type of matrix residual norm used can be
changed through the *Matrix residual norm type* card.


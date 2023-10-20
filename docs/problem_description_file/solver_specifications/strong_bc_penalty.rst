*********************************
Strong Boundary Condition Penalty
*********************************

::

	Strong Boundary Condition Penalty = <float>

-----------------------
Description / Usage
-----------------------

<float>
    **t,** a floating point number ( t :math:`>` 0.0 ) that specifies the penalty applied to strong boundary conditions
    

Currently only for `COLLOC_INT_SURF`, `STRONG_INT_SURF`, `STRONG_INT_NEDELEC`

------------
Examples
------------

A sample input card follows:
::

	Strong Boundary Condition Penalty = 1e12

-------------------------
Technical Discussion
-------------------------

The strong boundary condition penalty is used to enforce the boundary conditions on the surface by using a penalty 
to override the residual of equation with the boundary condition. This is the default behavior, see `Strong Boundary Condition Replace Equation` for an alternative.

--------------
References
--------------

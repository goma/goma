Eigen Linear Solver tolerance
==============================

.. code-block:: none

    Eigen Linear Solver tolerance = <float>

Description/Usage
-----------------

This optional card is used to specify the convergence tolerance used by Aztec linear 
solvers when solving the shifted linear system used by the ARPACK eigensolver, and is 
applicable only when an Aztec iterative solver is specified for the steady state problem. 
The absolute tolerance passed from ARPACK to Aztec is the product of this value and 
the norm of the mass matrix, which is recalculated at each iteration. Valid input are 
defined as:

LStol>=0.0
    The base linear solver convergence criterion used within the eigensolver.

The default value of LStol is the value specified for the linear solver (see the Residual 
Ratio Tolerance card) or its default of 1.0e-6.

Examples
--------

Here is a sample card:

::

    Eigen Linear Solver tolerance = 1.0e-10

Technical Discussion
--------------------

This card value differs from the one in the Solver Specifications section because it has 
a mass matrix norm factor. It can also be used to set a smaller tolerance for eigensolves 
than for steady state solves.

This card is not applicable to the eggroll eigensolver, which does not use Aztec.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

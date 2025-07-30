Eigen Tolerance
===============

.. code-block:: none

    Eigen Tolerance = <float>

Description/Usage
-----------------

This optional card specifies the absolute residual tolerance used in the eggroll 
eigensolver. Valid input are defined as:

tol>=0.0
    The termination criteria used within the eigensolver.

The default value of tol is 1.0e-6.

Examples
--------

Here is a sample card:

::

    Eigen Tolerance = 1.0e-16

Technical Discussion
--------------------

The value of tol specifies a termination condition used within the eigensolver. At the 
time of this writing it was difficult to guarantee that a small tol would lead to very 
accurate eigenpairs with the eggroll solver. This is one of the primary reasons for 
switching to the ARPACK eigensolver. The correspondence of "internal residual" to 
"error in eigenpair" is even more difficult with the generalized eigenvalue problem 
than "residual" to "error in solution vector" with regular linear solvers.

This card is not applicable to the ARPACK eigensolver, which uses a relative tolerance 
(see the Eigen Relative tolerance card).

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

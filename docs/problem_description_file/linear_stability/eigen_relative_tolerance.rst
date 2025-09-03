Eigen Relative tolerance
========================

.. code-block:: none

    Eigen Relative tolerance = <float>

Description/Usage
-----------------

This optional card specifies the relative residual tolerance used in the ARPACK 
eigensolver. Valid input are defined as:

tol>=0.0
    The relative termination criterion used within the eigensolver.

The default value of tol is 1.0e-6.

Examples
--------

Here is a sample card:

::

    Eigen Relative tolerance = 1.0e-10

Technical Discussion
--------------------

The value of tol specifies a termination condition used within the eigensolver which is 
applied to the transformed eigensystem. The correspondence of "internal residual" to 
"error in eigenpair" is even more difficult with the generalized eigenvalue problem 
than "residual" to "error in solution vector" with regular linear solvers.

This card is not applicable to the eggroll eigensolver, whose tolerance is specified with 
the Eigen Tolerance card.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

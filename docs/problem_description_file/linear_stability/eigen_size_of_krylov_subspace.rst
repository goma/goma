Eigen Size of Krylov subspace
==============================

.. code-block:: none

    Eigen Size of Krylov subspace = <int>

Description/Usage
-----------------

This optional card specifies the dimension of the Krylov subspace used in the 
eigensolver. Valid inputs are:

m>0
    A positive integer specifying the maximum dimension of the Krylov subspace.

The default value of m is 30.

Examples
--------

Here is a sample card:

::

    Eigen Size of Krylov subspace = 120

Technical Discussion
--------------------

A restarted Lanczos-type algorithm is used internally to solve the generalized 
eigenvalue problem. Specifying a very large m can lead to wasted effort, while an m
that is too small can lead to nonconvergence, or possibly very slow convergence. The 
optimal value of m is completely problem dependent.

This card is applicable to both eggroll and ARPACK eigensolvers. When ARPACK is 
used, however, this value is also used as the maximum number of inner eigensolver 
iterations (see the Eigen Maximum Iterations card).

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

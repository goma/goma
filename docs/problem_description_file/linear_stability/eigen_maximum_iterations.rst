Eigen Maximum Iterations
========================

.. code-block:: none

    Eigen Maximum Iterations = <int>

Description/Usage
-----------------

This optional card specifies the maximum number of eggroll eigensolver iterations. Its 
valid argument is given by:

n>0
    Maximum number of eigensolver iterations.

The default value of n is 100.

Examples
--------

Here is a sample card:

::

    Eigen Maximum Iterations = 32

Technical Discussion
--------------------

Within the eigensolver is a restarted Lanczos-type method. The maximum number of 
restarts is n. Also, see the Eigen Size of Krylov subspace card.

This card is not applicable to the ARPACK eigensolver, for which the maximum 
number of iterations is set to the Krylov subspace size.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

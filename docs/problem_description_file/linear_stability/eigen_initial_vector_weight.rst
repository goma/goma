Eigen Initial Vector Weight
============================

.. code-block:: none

    Eigen Initial Vector Weight = <float>

Description/Usage
-----------------

This optional card allows the user to modify the initial vector used within the eggroll 
eigensolver to contruct a Krylov subspace. Legal inputs are:

w
    A scale factor for the amount of uniform randomness to apply to each component of the initial vector.

The default value of w is 0.5.

Examples
--------

Here is a sample card:

::

    Eigen Initial Vector Weight = 0.1

Technical Discussion
--------------------

Because of the way a Krylov space is generated, it is very important that the initial 
vector contain some component in the direction of the eigenvector. By adding a 
random vector to the initial vector this can be guaranteed. A random vector r is 
constructed to have components from U[0,1] (uniform distribution). It is then 
normalized to have length 1. The initial vector is then modifed as x = w * r + (1-w) * x. 
The effect this value has should be minimal for w not equal to 0 or 1.

This card is not applicable to the ARPACK eigensolver.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

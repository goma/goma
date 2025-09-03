Eigen Initial Shifts
====================

.. code-block:: none

    Eigen Initial Shifts = <float1> <float2> <float3> <float4>

Description/Usage
-----------------

This optional card allows the user to supply a set of initial shifts for the shift-and-invert 
algorithm used in the eggroll eigensolver. Valid syntax is:

f1, f2, f3, f4
    All real values that should be close to the eigenvalue the user is trying to compute.

The default value(s) are all -1.0.

Examples
--------

Here is a sample card:

::

    Eigen Initial Shifts = -0.001 -0.001 -0.001 -0.001

Technical Discussion
--------------------

The eigensolver algorithm beings by trying to find eigenvalues near f1. If it fails to do 
so, it will try for eigenvalues near f2, and so on. The only time when shifts other than 
f1 are used is when f1 is a truly bad approximation to an eigenvalue. This happens 
rarely (the locations of the eigenvalues need to be clustered in a particularly bad way).

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

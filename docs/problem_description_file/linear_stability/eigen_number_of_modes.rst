Eigen Number of modes
=====================

.. code-block:: none

    Eigen Number of modes = <int>

Description/Usage
-----------------

This optional card specifies how many modes (eigenvalue/eigenvector pairs) are to be 
computed during linear stability analysis. The input values are:

N
    Number of eigenvalue/eigenvector pairs (modes) to compute.

The default value for N is 10.

Examples
--------

Here is a sample card:

::

    Eigen Number of modes = 10

Technical Discussion
--------------------

When the Linear Stability card indicates stability analysis should be performed, this 
card is used when computing the spectrum. Generally, the user is indicating that they 
would like to compute the leading N modes. Judicious selection of other pameters is 
required for this to happen. See the Eigen Initial Shifts card, especially.

This card is applicable to both eggroll and ARPACK eigensolvers.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

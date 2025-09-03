Eigen Cayley Sigma
==================

.. code-block:: none

    Eigen Cayley Sigma = <float>

Description/Usage
-----------------

This optional card is used by the ARPACK eigensolver as a shift parameter. When the 
Cayley transformation is selected (see the Eigen Algorithm card), this value and Eigen 
Cayley Mu are the two shift parameters. When the shift and invert (default) algorithm is 
used, this is the single shift parameter.

The default value is 100.

Examples
--------

Here is a sample card:

::

    Eigen Cayley Sigma = 10.0

Technical Discussion
--------------------

When the Cayley algorithm is selected, this value and Eigen Cayley Mu also determine 
which of the two Cayley transformation methods (A or B) is used by ARPACK. For 
either of these methods, sigma must be to the right of (greater than) the real parts of all 
eigenvalues -- any eigenvalues greater than sigma will not be found and will generate 
an error message. If this happens, just increase sigma and try again.

This card is not applicable to the eggroll eigensolver, which takes shift values from the 
Eigen Initial Shifts card.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

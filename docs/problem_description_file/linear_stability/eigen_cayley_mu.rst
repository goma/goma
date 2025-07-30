Eigen Cayley Mu
===============

.. code-block:: none

    Eigen Cayley Mu = <float>

Description/Usage
-----------------

This optional card is used by the ARPACK eigensolver as the second Cayley shift 
parameter. When the Cayley transformation is selected (see the Eigen Algorithm card), 
this value and Eigen Cayley Sigma are the two shift parameters. When the shift and 
invert (default) algorithm is used, this card is not applicable.

The default value is 1000.

Examples
--------

Here is a sample card:

::

    Eigen Cayley Mu = 10.0

Technical Discussion
--------------------

When the Cayley algorithm is selected, this value and Eigen Cayley Sigma also 
determine which of the two Cayley transformation methods (A or B) is used by 
ARPACK.

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

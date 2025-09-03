Eigen Algorithm
===============

.. code-block:: none

    Eigen Algorithm = {cayley | si}

Description/Usage
-----------------

This optional card is used by the ARPACK eigensolver to select between the Cayley 
transformation algorithm and the shift and invert algorithm. When cayley is selected, 
Eigen Cayley Sigma card and Eigen Cayley Mu are the two shift parameters. When the 
shift and invert (default) algorithm is used, Eigen Cayley Sigma is the single shift 
parameter. The valid options are:

cayley
    Cayley transformation algorithm.

si
    Shift and invert algorithm.

The default value is si.

Examples
--------

Here is a sample card:

::

    Eigen Algorithm = cayley

Technical Discussion
--------------------

When the Cayley algorithm is selected, the relative values of Eigen Cayley Sigma and 
Eigen Cayley Mu also determine which of the two Cayley transformation methods (A 
or B) is used by ARPACK.

Another use of this card is to determine whether eggroll or ARPACK will be used for 
non-continuation problems: If this card is present with a valid option, ARPACK will 
be used; otherwise eggroll will be used. Continuation runs with LOCA will always use 
ARPACK, and non-LOCA continuation runs cannot be combined with stability 
analysis.

This card is not applicable to the eggroll eigensolver, which uses only the shift and 
invert algorithm.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

Eigen Recycle
=============

.. code-block:: none

    Eigen Recycle = { yes | no }

Description/Usage
-----------------

This optional card determines whether or not the eigenpairs will be recycled within the 
eggroll eigensolver algorithm. Valid input is defined by:

yes
    Do recycle eigenpairs.

no
    Do not recycle eigenpairs.

The default value is no.

Examples
--------

Here is a sample card:

::

    Eigen Recycle = yes

Technical Discussion
--------------------

In practice, this card seems to have little effect.

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

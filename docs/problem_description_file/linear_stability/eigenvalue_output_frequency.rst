Eigenvalue output frequency
===========================

.. code-block:: none

    Eigenvalue output frequency = <int>

Description/Usage
-----------------

This optional card is used by LOCA to determine how often to call the ARPACK 
eigensolver during a continuation problem. When this value is selected to n, ARPACK 
will be called on the first continuation step and every nth step thereafter. n can also be 
set to -1, in which case ARPACK is called only after the last continuation step.

The default value is -1, in which case ARPACK is called on the last continuation step.

Examples
--------

Here is a sample card:

::

    Eigenvalue output frequency = 2

Technical Discussion
--------------------

This card is useful for continuation problems when it is not desired to calculate 
eigenvalues at every continuation step, enabling time to be saved.

This card is not applicable to the eggroll eigensolver, which cannot be used during 
continuation.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

Eigen Matrix Output
===================

.. code-block:: none

    Eigen Matrix Output = { yes | no }

Description/Usage
-----------------

This optional card is used to indicate that the Jacobian and mass matrices which are 
calculated for linear stability analysis should be written to output files.

Examples
--------

To write the Jacobian and mass matrices when doing LSA, use:

::

    Eigen Matrix Output = yes

Technical Discussion
--------------------

These files can be quite large, so it is advisable to output them only if you need them. 
The intent of matrices availble is so that the experienced human eigensolver can import 
the matrices into other tools (i.e., Matlab) to perform their own investigations. Note 
that when the Linear Stability choice is "file" or "3Dfile", this card must be set to yes in 
order to create and write the desired stability output files.

This card is applicable to both eggroll and ARPACK. The same matrix output function 
is called in either case.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

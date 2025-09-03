Eigenvector output frequency
============================

.. code-block:: none

    Eigenvector output frequency = <integer>

Description/Usage
-----------------

This optional card is used by LOCA to determine how often to write out a set of 
eigenvectors from the ARPACK eigensolver to EXODUS II files during a continuation 
problem. When this value is selected to n, ARPACK eigenvectors will be written on 
the first continuation step and every nth step thereafter. n can also be set to -1, in which 
case eigenvectors will be output only on the last continuation step.

The default value is -1, in which case eigenvectors are written on every continuation 
step.

Examples
--------

Here is a sample card:

::

    Eigenvector output frequency = 2

Technical Discussion
--------------------

This card is useful for continuation problems when is is not desired to save 
eigenvectors at every continuation step, thus enabling savings of time and disk space.

This value is independent of the Eigenvalue output frequency, so if a frequency m is 
specified for eigenvalues and n is not a multiple of m, then eigenvectors will be written 
on step numbers which are common multiples of m and n.

When requesting linear stability analysis during a LOCA continuation run, the 
eigenvectors from all continuation steps corresponding to a given mode number (or 
mode number/wave number combination for 3D of 2D) will be written sequentially to 
the same EXODUS II file. The time stamp for each step will be the eigenvalue, rather 
than the continuation parameter value. Complex eigenvectors are written as two steps 
in the same file: first the real part (time-stamped with the real part of the eigenvalue), 
then the imaginary part (time-stamped with the imaginary part of the eigenvalue).

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

********************
**Von Mises Strain**
********************

::

   Von Mises Strain = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This option allows you to plot the Von Mises strain invarients of the strain tensor, for
use with the FAUX_PLASTICITY model of the modulus. These quantities are written
to the *Output EXODUS II* file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the von Mises strain invariants and 
         output as a nodal variable in the Output 
         EXODUS II file.
**no**   Do not calculate the invariants (default).
======== ===============================================

------------
**Examples**
------------

The following sample card requests **LAMBDA** be written to the EXODUS II file:
::

   Von Mises Strain = yes

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
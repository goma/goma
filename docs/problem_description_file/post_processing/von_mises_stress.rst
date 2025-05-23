********************
**Von Mises Stress**
********************

::

   Von Mises Stress = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This option allows you to plot the Von Mises stress tensor invarients, for use with the
FAUX_PLASTICITY model of the modulus. These quantities are written to the
*Output EXODUS II* file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the von Mises stress invariants and 
         output as a nodal variable in the Output 
         EXODUS II file.
**no**   Do not calculate the invariants (default).
======== ===============================================

------------
**Examples**
------------

The following sample card requests **LAMBDA** be written to the EXODUS II file:
::

   Von Mises Stress = yes

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




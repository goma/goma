***************
**Lame LAMBDA**
***************

::

   Lame LAMBDA = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This option allows you to plot the Lame LAMDA mechanical property, which is
written to the *Output EXODUS II* file as the variable **LAMBDA**. This is a useful
feature for temperature dependent mechanical properties and the like. Contouring this
variable **LAMBDA** over the domain can be useful in explaining some physical
phenomena.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the Lame LAMBDA and output as a nodal
         variable in the Output EXODUS II file.
**no**   Do not calculate the the coefficient (default).
======== ===============================================

------------
**Examples**
------------

The following sample card requests **LAMBDA** be written to the EXODUS II file:
::

   Lame LAMBDA = yes

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
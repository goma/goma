***********
**Density**
***********

::

   Density = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This card is used to trigger the thermophysical or mechanical property of density (see
*Density* card) to be computed and output as an EXODUS II nodal variable in the Ouput
EXODUS II file with the variable name **RHO**.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the density and store it as a nodal 
         variable in the output EXODUS II file.
**no**   Do not calculate density.
======== ===============================================

------------
**Examples**
------------

This is an example of the input to request density be written to the EXODUS II file.
::

   Density = yes

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
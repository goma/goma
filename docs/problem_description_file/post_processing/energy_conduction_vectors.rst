*****************************
**Energy Conduction Vectors**
*****************************

::

   Energy Conduction Vectors = {yes | no}

-----------------------
**Description / Usage**
-----------------------

Activation of this option can be used to visualize the energy conduction paths in a
solution. The resulting nodal variables are called **TCOND0** (conduction in x direction),
**TCOND1** (conduction in y direction), and **TCOND2** (conduction in z direction) in the
output EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the energy conduction vectors and 
         store in the output EXODUS II file.
**no**   Do not calculate the energy conduction vectors.
======== ===============================================

------------
**Examples**
------------

This example card requests that energy conduction vectors be written to the EXODUS
II file:
::

   Energy Conduction Vectors = yes

-------------------------
**Technical Discussion**
-------------------------

These quantities can be employed in a hedge-hog or vector plot to visualize the energy
conduction pathways across a domain (cf. the vector option in BLOT, or the hedge-hog
option in Mustafa, for example).



--------------
**References**
--------------

No References.
*************************
**Lagrangian Convection**
*************************

::

   Lagrangian Convection = {yes | no}

-----------------------
**Description / Usage**
-----------------------

In deformable solids with a Lagrangian mesh, convection in the stress-free state can be
mapped to the deformed configuration; this variable stores the velocity vectors of this
solid motion (see *Convective Lagrangian Velocity* card). This variable is called **VL1, VL2, VL3** in the output EXODUS II file; the integer values **1, 2** and **3** denote
coordinate directions.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the Lagrangian convection and store 
         as nodal variable velocity fields in the output 
         EXODUS II file.
**no**   Do not calculate the Lagrangian convection.
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   Lagrangian Convection = no

-------------------------
**Technical Discussion**
-------------------------

This option only applies to *Mesh Motion* type of *LAGRANGIAN*.



--------------
**References**
--------------

No References.
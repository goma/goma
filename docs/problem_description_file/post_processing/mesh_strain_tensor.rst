**********************
**Mesh Strain Tensor**
**********************

::

   Mesh Strain Tensor = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The mesh strain tensor is associated with the equations of elasticity. The strain tensor
has six entries (in three dimensions, because it is symmetric) called **E11, E22, E33, E12, E13,** and **E23** in the output EXODUS II file, corresponding to the six independent
components (the numbers **1, 2,** and **3** indicate the basis direction, e.g. 1 means xdirection
for a Cartesian system).

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the mesh strain tensor and write the 
         components as nodal variables to the output 
         EXODUS II file.
**no**   Do not calculate the mesh strain tensor.
======== ===============================================

------------
**Examples**
------------

The following example input card does not request output of the stain tensor:
::

   Mesh Strain Tensor = no

-------------------------
**Technical Discussion**
-------------------------

Definitions of the strain tensor depend on the solid constitutive equation type (see
description for *Solid Constitutive Equation* card).



--------------
**References**
--------------

No References.
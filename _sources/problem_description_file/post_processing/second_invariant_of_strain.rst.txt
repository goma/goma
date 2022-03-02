******************************
**Second Invariant of Strain**
******************************

::

   Second Invariant of Strain = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The strain tensor is associated with the deformation of the mesh. Its second invariant
indicates the level of shear strain of the mesh. This variable is called **IIE** in the output EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the second invariant.
**no**   Do not calculate the second invariant.
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   Second Invariant of strain = yes

-------------------------
**Technical Discussion**
-------------------------

The second invariant is computed in *Goma* as

.. figure:: /figures/323_goma_physics.png
	:align: center
	:width: 90%

Here Einstein’s summation convention applies, viz.

.. figure:: /figures/324_goma_physics.png
	:align: center
	:width: 90%




.. 
	TODO - Lines 38 and 44 are photos that need to be swapped with the equations.
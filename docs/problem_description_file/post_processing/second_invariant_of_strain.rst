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

.. math::

   II_E=\frac12(E_{ij}E_{ij}-E_{ii}E_{jj})

Here Einstein’s summation convention applies, viz.

.. math::

   E_{ij}E_{ij}=\sum_i\sum_j E_{ij}E_{ij}

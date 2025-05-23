*****************************
**First Invariant of Strain**
*****************************

::

	First Invariant of Strain = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The strain tensor is associated with the deformation of the mesh. Its first invariant is its
trace and represents the volume change in the small strain limit. This variable is called
**IE** in the output EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the first invariant.
**no**   Do not calculate the first invariant.
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   First Invariant of strain = yes

-------------------------
**Technical Discussion**
-------------------------

Computation of the strain tensor :math:`\underline{E}` in *Goma* is discussed on the *Solid Constitutive Equation* card. The trace is related to the divergence of the tensor, and hence related to a measure of volume change in a material.

It should be noted that the mesh strain is equivalent to the material strain for
*LAGRANGIAN* mesh motion types. For *ARBITRARY* or *TOTAL_ALE* mesh motion
types (see *Mesh Motion* card), the strain is strictly related to mesh and not the material.



--------------
**References**
--------------

No References.
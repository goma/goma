*******************
**Stream Function**
*******************

::

	Stream Function = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The stream function provides a visual representation of the flow field in incompressible
fluids and is derived from the fluid velocity components identified in the *Output
Exodus II file* card.

This auxiliary field triggered by “yes” on this card results in a nodal variable that is
called **STREAM** in the output EXODUS II file.

The permissible values for this postprocessing option are:

======== =====================================
**yes**  Calculate the stream function.
**no**   Do not calculate the stream function.
======== =====================================

------------
**Examples**
------------

Following is a sample card:
::

   Stream Function = no

-------------------------
**Technical Discussion**
-------------------------

This function is computed with an element-by-element volumetric flow calculation
routine. Poor element quality can result in “kinks” in the stream function field when
contoured.

It is important to construct a mesh whose elements are contiguously ordered in such a
way that there are no isolated clusters as the elements are swept, i.e., element n+1 must
be in contact with one of the previous n elements. NOTE: as of 4/2001 an automatic
element reordering scheme based on Reverse Cuthill-McKee algorithm has been
implemented in *Goma*. Automatic ordering can be assured by issuing the OPtimize
command to the FASTQ meshing module (cf. Blacker 1988). Most other mesh
generators do not provide this service, viz. they do not put out an element order-map
field in the EXODUS II file.

NOTE: THIS FUNCTION IS NOT AVAILABLE IN THREE DIMENSIONS, but
pathlines, which are equivalent to streamlines for steady flows can be computed in
many graphics packages, like Mustafa (Glass, 1995).



--------------
**References**
--------------

SAND88-1326: FASTQ Users Manual: Version 1.2, Sandia Technical Report, Blacker,
T. D. 1988.

Mustafa, Glass, M. W., Personal Communication, 1995
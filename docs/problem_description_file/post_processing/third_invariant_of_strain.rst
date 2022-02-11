*****************************
**Third Invariant of Strain**
*****************************

::

   Third Invariant of strain = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The strain tensor is associated with the deformation of the mesh. Its third invariant
indicates the volume change from the stress-free state (**IIIE** = 1.0 indicates no volume
change). This variable is called **IIIE** in the output EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the third invariant.
**no**   Do not calculate the third invariant.
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   Third Invariant of strain = yes

-------------------------
**Technical Discussion**
-------------------------

The mathematical definition of the third invariant is related to the determinant of the
strain tensor, which is defined for the various constitutive equations in the manual entry
for the *Solid Constitutive Equation* card.



--------------
**References**
--------------

No References.
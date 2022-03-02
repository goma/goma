*******************
**Stress Contours**
*******************

::

	Stress contours = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This card allows the user to invoke the components of the stress tensor for all
viscoelastic modes be included as nodal post-processing variables. Often times this is
not desirable on long time-dependent runs because of the voluminous data that will
appear in the output EXODUS II file. This variable is called *csij_mode* in the output
EXODUS II file, where *i* and *j* indicate components of the stress tensor and *mode*
indicates the desired viscoelastic mode; for example, **cs23_4** represents the stress
contour for the fifth mode of polymer stress component *yz*.

The permissible values for this postprocessing option are:

======== =======================================================================
**yes**  Calculate and include the stress-tensor components for all
         modes of viscoelasticity.
**no**   Do not calculate and include the stress-tensor components.
======== =======================================================================

These stresses become dependent variables if the *Polymer Constitutive Equation* card
is given any model but the *NOPOLYMER* model.

------------
**Examples**
------------

An example card requesting viscoelastic stress components be written:
::

   Stress contours = yes

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

GT-014.1: Tutorial for Running Viscoelastic Flow Problems with GOMA, June 21,
2000, R. R. Rao
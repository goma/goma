*****************
**Fill Contours**
*****************

::

	Fill contours = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This card triggers the inclusion of the level set or VOF fill function as a nodal variable
in the output EXODUS II file.

The nodal variable appears as **FILL** in the output EXODUS II file. This function is
computed with the FILL equation (see EQ card).

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the fill contours.
**no**   Do not calculate the fill contours.
======== ===============================================

------------
**Examples**
------------

An example card requesting **FILL** contours be written to the EXODUS II file is:
::

   Fill contours = yes

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer
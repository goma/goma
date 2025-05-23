*********************
**Porous Saturation**
*********************

::

   Porous Saturation = {yes | no}

-----------------------
**Description / Usage**
-----------------------

In partially saturated porous media, the saturation represents the volume fraction of the
pore space that is filled with liquid. If this option is selected, then the saturation level
(an auxiliary variable) is included in the output EXODUS II file. This variable is called
**SAT** in the output EXODUS II file.

The permissible values for this postprocessing option are:

============= ================================================================
**yes**       Calculate the porous saturation and write to the output
              EXODUS II file.
**no**        Do not calculate the porous saturation.
============= ================================================================

------------
**Examples**
------------

This sample input card turns off writing of saturation to the EXODUS II file:
::

   Porous Saturation = no

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

SAND96-2149: Drying in Deformable Partially-Saturated Porous Media: Sol-Gel
Coatings, Cairncross, R. A., P. R. Schunk, K. S. Chen, S. S. Prakash, J. Samuel, A. J.
Hurd and C. Brinker (September 1996)
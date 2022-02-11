******************
**Mass Fluxlines**
******************

::

   Mass Fluxlines = {yes | no}

-----------------------
**Description / Usage**
-----------------------

With this post-processing option mass-diffusion pathlines are calculated and stored as
post-processing nodal variables in the output EXODUS II file. This variables are called
**Y0FLUX**, **Y1FLUX**, . . .(by species number) in the file and can be contoured in the
visualization program. These flux lines are analogous to the stream function, viz.
contours of the flux function represent pathlines for each species in solution.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the mass fluxlines and include in the
         output EXODUS II file.
**no**   Do not calculate the mass fluxlines.
======== ===============================================

------------
**Examples**
------------

The following sample card requests that mass-diffusion pathlines be written to the
EXODUS II file:
::

   Mass Fluxlines = yes

-------------------------
**Technical Discussion**
-------------------------

Currently this option is available only for *FICKIAN* and *HYDRODYNAMIC* mass flux
types (see *Diffusion Constitutive Equation* card). In the *FICKIAN* case, the flux is
computed with the base, constant diffusivity. Also, the *Mass Diffusion Vectors* post
processing option must also be activated for this option to work.



--------------
**References**
--------------

No References.
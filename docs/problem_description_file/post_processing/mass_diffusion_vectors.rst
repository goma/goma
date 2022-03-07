**************************
**Mass Diffusion Vectors**
**************************

::

   Mass Diffusion Vectors = {yes | no}

-----------------------
**Description / Usage**
-----------------------

Activating this post-processing option allows the user to visualize the diffusive mass
flux directions of all species components in a problem. Species components result from
the *EQ = species_bulk* equation card. With this option selected, the output EXODUS II
file will contain nodal variables called **Y0dif0** (diffusion of first species in x direction),
**Y0dif1** (diffusion of first species in y direction), **Y0dif2** (diffusion of first species in z
direction), **Y1dif0** (diffusion of second species in x direction), **Y2dif1** (diffusion of
second species in y direction), . . . and so on, depending on the number of species
components in the problem.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the mass diffusion vectors and 
         include in the output EXODUS II file.
**no**   Do not calculate the mass diffusion vectors.
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   Mass Diffusion Vectors = yes

-------------------------
**Technical Discussion**
-------------------------

Currently this option is available only for *FICKIAN* and *HYDRODYNAMIC* mass flux
types (see *Diffusion Constitutive Equation* card). In the *FICKIAN* case, the flux is
computed with the base, constant diffusivity.



--------------
**References**
--------------

No References.
**************************************
**Grid Peclet Number in Porous Media**
**************************************

::

   Grid Peclet Number in porous media = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This option triggers the computation and output of the so-called grid-level Peclet
number as a nodal variable in the output EXODUS II file. It appears as a nodal variable
called **Por_Grid_Peclet**. This quantity gives the user a measure of advective transport
relative to diffusive transport in a porous medium, and is strongly correlated to the
steepness of a saturation front. This quantity is actually used to scale the formulation
which employs the streamline upwind Petrov-Galerkin method for stabilizing the
equations for partially saturated flow. This option only applies for unsaturated media
and only for the *SUPG* option on the *Porous Weight Function* card.

The permissible values for this postprocessing option are:

============= ================================================================
**yes**       Compute the grid-level Peclet Number and write to output
              EXODUS II file.
**no**        Do not calculate the grid-level Peclet Number.
============= ================================================================

------------
**Examples**
------------

This is a sample input card to activate calculation of the Peclet Number:
::

   Grid Peclet Number in porous media = yes

-------------------------
**Technical Discussion**
-------------------------

See discussion for the *Porous Weight Function* card.



--------------
**References**
--------------

GTM-029.0: SUPG Formulation for the Porous Flow Equations in Goma, H. K.
Moffat, August 2001 (DRAFT).
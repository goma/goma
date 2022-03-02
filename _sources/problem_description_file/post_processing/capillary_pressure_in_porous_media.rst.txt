**************************************
**Capillary Pressure in Porous Media**
**************************************

::

   Capillary pressure in porous media = {yes | no}

-----------------------
**Description / Usage**
-----------------------

In partially saturated porous media, the capillary pressure is the difference between the
gas and liquid pressures. This option only takes affect for *POROUS_TWO_PHASE* and
*POROUS_UNSATURATED* media types (see *Media Type* card). This variable is called
**PC** in the output EXODUS II file.

The permissible values for this postprocessing option are:

============= ================================================================
**yes**       Calculate the capillary pressure and write to output
              EXODUS II file.
**no**        Do not calculate the capillary pressure.
============= ================================================================

------------
**Examples**
------------

This is a sample input card to activate calculation of capillary pressure:
::

   Capillary pressure in porous media = yes

-------------------------
**Technical Discussion**
-------------------------

The capillary pressure is a critical variable for partially saturated porous media, and is
in fact the dependent variable for unsaturated (not two-phase) flows for which the gasphase
pressure is taken as constant. It is simply defined as

.. figure:: /figures/335_goma_physics.png
	:align: center
	:width: 90%

As such, positive capillary pressures imply liquid phase pressure being greater than gas
phase pressure. Because liquid phase saturation strongly correlates to capillary
pressure, this current quantity is a good indicator of the level of liquid inventory in
smaller pores in the skeleton relative to large pores. Contouring this quantity can give
some indication of the level of suction exerted on the porous-skeleton, which is
relevant when the skeleton is taken as deformable.



--------------
**References**
--------------

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

.. 
	TODO - Line 43 is a photo that needs to be swapped with the equation.
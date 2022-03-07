*******************
**POROUS_PRESSURE**
*******************

::

	BC = POROUS_PRESSURE SS <bc_id> <integer1> <integer2>

-----------------------
**Description / Usage**
-----------------------

**(PCC/POR_LIQ_PRES)**

This condition enforces a continuous fluid-phase pressure between material types, and
is applied to a side set between two materials, one of type *POROUS_SATURATED,
POROUS_UNSATURATED*, or *POROUS_TWO_PHASE*, and the other of type
*CONTINUOUS* (see material card *Media Type*). Basically it sets the continuity of
hydrodynamic pressure in the continuous fluid to the liquid Darcy pressure in the
porous medium, at the interface. The input data is as follows:

=================== ===========================================================
**POROUS_PRESSURE** Name of the boundary condition (<bc_name>).
**SS**              Type of boundary condition (<bc_type>), where **NS**
                    denotes side set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (node set in
                    EXODUS II) in the problem domain.
<integer1>          Element block ID of the porous phase medium.
<integer2>          Element block ID of the continuous fluid phase medium.
=================== ===========================================================

------------
**Examples**
------------

An example input card for this boundary condition follows:
::

   BC = POROUS_PRESSURE NS 101   1   2

This card sets the Darcy liquid phase pressure (*p_liq* in the output EXODUS II file) in
element block 1 equal to the continuous phase hydrodynamic pressure (*P* in the output
EXODUS II file) in element block 2.

-------------------------
**Technical Discussion**
-------------------------

The mathematical form for this boundary condition is as follows

.. figure:: /figures/176_goma_physics.png
	:align: center
	:width: 90%

but its implementation is not; a memo describing the details of this boundary condition
and how it is applied is cited below. This continuity of pressure is critical for the
sensitivity of pressurizing the continuos phase to the penetration rate of the porous
phase. Interestingly, it forces one to set the pore-phase pressure datum to the same
datum in the continuous phase, and that effects the level of the Saturation versus
capillary pressure curve (see *Saturation* material card).



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

.. 
	TODO - Line 52 has a photo in which the equation that needs to be written.
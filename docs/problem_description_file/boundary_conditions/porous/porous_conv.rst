***************
**POROUS_CONV**
***************

::

	BC = POROUS_CONV SS <bc_id> <integer>

-----------------------
**Description / Usage**
-----------------------

**(WIC/POR_LIQ_PRES)**

This boundary condition is used to set the total flux of the liquid phase solvent (in both the gas and liquid phase) at the surface of a *POROUS_UNSATURATED* or
*POROUS_TWO_PHASE* medium to the net convection of solvent due to a
superimposed convective Lagrangian velocity (see *Media Type* card and *Convective
Lagrangian Velocity* card). The only input is an integer indicating which component of
the liquid phase solvent is to be set (*as of 11/2/01 this component selectability option is not available and as indicated below should be set to zero; this card has not been tested*).

Definitions of the input parameters are as follows:

=============== ==========================================================
**POROUS_CONV** Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS**
                denotes side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set
                in EXODUS II) in the problem domain.
<integer>       Species number of transported species. (currently only
                used for multicomponent species in the phases, which as
                of 11/2/01 is not active, so set to zero).
=============== ==========================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = POROUS_CONV SS 12 0

that applies a convective flux to side set 12 for porous liquid phase species 0. This
species number is currently not used and ignored.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition has the following form

.. figure:: /figures/174_goma_physics.png
	:align: center
	:width: 90%

where the left hand side is the total flux of the solvent i in the medium, which includes,
in order, the flux due to Darcy flow of gas vapor, the Darcy flow of liquid solvent, the
diffusive flux of gas vapor in the pore space and the diffusive flux of liquid solvent in
the liquid phase. :math:`v_s` is the user supplied convection velocity of the stress-free state as
defined on the *Convective Lagrangian Velocity* card. As of now (11/2/01), this
condition is used for a single component liquid solvent and has not been furbished for a
single component of that solvent. Also, as of 11/02/01 the condition has not been
tested.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

..
	TODO _ Line 53 has a photo that needs to be changed into an equation.
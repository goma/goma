*************************
**POROUS_LIQ_FLUX_CONST**
*************************

::

	BC = POROUS_LIQ_FLUX_CONST SS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(WIC/POR_LIQ_PRES)**

This boundary condition sets the flux of liquid-phase solvent to a constant value in the
Galerkin finite element weak sense. Specifically, this flux is applied to a side set as a
weak-integrated constant and will set the net flux of liquid phase solvent component (in
both gas and liquid phases) to a specified value. It can be applied to material regions of
*Media Type POROUS_SATURATED, POROUS_UNSATURATED*, and
*POROUS_TWO_PHASE* (see Technical Discussion below).

Definitions of the input parameters are as follows:

========================== ======================================================
**POROUS_LIQ_FLUX_CONST**  Name of boundary condition (<bc_name>).
**SS**                     Type of boundary condition (<bc_type>), where **SS**
                           denotes side set in the EXODUS II database.
<bc_id>                    The boundary flag identifier, an integer associated with
                           <bc_type> that identifies the boundary location (side set
                           in EXODUS II) in the problem domain.
<float1>                   Value of the liquid-solvent total flux, in 
                           M/ :math:`L^2` -t.
[float2]                   This optional parameter is not applicable to this
                           boundary condition type, even though it is parsed if
                           present. This parameter is used for boundary conditions
                           of the Dirichlet type.
========================== ======================================================

------------
**Examples**
------------

The input card
::

   BC = POROUS_LIQ_FLUX_CONST SS 102 200.0

sets the total liquid-solvent mass flux, in both gas and liquid phases, to 200.0 along the side set 102.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is of the mathematical form:

.. figure:: /figures/171_goma_physics.png
	:align: center
	:width: 90%

where :math:`v_s` is the user supplied convection velocity of the stress-free state as defined on
the *Convective Lagrangian Velocity* card (this is usually zero except in advanced
cases), :math:`p_l^T` is the total bulk density of liquid phase solvent (in both gas and liquid
phase, and hence depends on the local saturation), :math:`p_l` is the pure liquid density, :math:`\phi` is the
porosity, :math:`p_l` is the liquid phase pressure, and the other quantities on the second term
help define the Darcy velocity. The *const* quantity is the input parameter identified
above (<float1>). Note that this sets the flux relative to the boundary motion to the
*const* value, but by virtue of the Galerkin weak form this condition is automatically
applied with *const* =0 if no boundary condition is applied at the boundary. In a saturated
case, viz. *POROUS_SATURATED* media type, this condition is applied as

.. figure:: /figures/172_goma_physics.png
	:align: center
	:width: 90%



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

.. TODO - Lines 56 and 71 have photos that need to be replaced with the proper equations.

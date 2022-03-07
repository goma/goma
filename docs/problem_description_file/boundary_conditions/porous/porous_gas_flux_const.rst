*************************
**POROUS_GAS_FLUX_CONST**
*************************

::

	BC = POROUS_GAS_FLUX_CONST SS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(WIC/POR_GAS_PRES)**

This boundary condition card is used to set the flux of gas-phase solvent to a constant
value in the Galerkin finite element weak sense. Specifically, this flux is applied to a
side set as a weak-integrated constant and will set the net flux of gas phase solvent
component (in both gas and liquid phases, but because the gas solvent is assumed
insoluble in the liquid phase, the liquid phase portion vanishes) to a specified value.
This boundary condition can be applied to material regions of *Media Type
POROUS_TWO_PHASE* only, as only this type contains a field of gas-phase solvent
flux. (See technical discussion below).

Definitions of the input parameters are as follows:

========================= =======================================================
**POROUS_GAS_FLUX_CONST** Name of the boundary condition (<bc_name>).
**SS**                    Type of boundary condition (<bc_type>), where **SS**
                          denotes side set in the EXODUS II database.
<bc_id>                   The boundary flag identifier, an integer associated with
                          <bc_type> that identifies the boundary location (side set
                          in EXODUS II) in the problem domain.
<float1>                  Value of the gas-solvent total flux, in 
                          M/ :math:`L^2` -t.
[float2]                  This optional parameter is not applicable to this
                          boundary condition type, even though it is parsed if
                          present. This parameter is used for boundary conditions
                          of the Dirichlet type.
========================= =======================================================

------------
**Examples**
------------

The input card
::

   BC = POROUS_LIQ_FLUX_CONST SS 102 200.0

sets the total gas-solvent mass flux, in the gas phase only, to 200.0 along the side set
102.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is of the mathematical form:

.. figure:: /figures/173_goma_physics.png
	:align: center
	:width: 90%

where :math:`v_s` is the user supplied convection velocity of the stress-free state as defined on
the *Convective Lagrangian Velocity* card (this is usually zero except in advanced
cases), :math:`p_g^T` is the total bulk density of gas phase solvent (in both gas and liquid phase,
and hence depends on the local saturation), :math:`p_g` is the pure gas density, 
:math:`\phi` is the porosity, :math:`p_g`
is the gas-phase pressure, and the other quantities on the second term help define the
Darcy velocity. The *const* quantity is the input parameter described above (<float1>).
Note that this sets the flux relative to the boundary motion to the *const* value, but by virtue of the Galerkin weak form this condition is automatically applied with *const = 0* if no boundary condition is applied at the boundary.




..
	TODO - Line 59 has a photo that needs to be changed into an equation.
***************
**POROUS_FLUX**
***************

::

	BC = POROUS_FLUX SS <bc_id> <integer> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/POR_LIQ_PRES)**

This boundary condition is used to set the total flux of the liquid phase solvent (in both the gas and liquid phase) at the surface of a *POROUS_UNSATURATED* or
*POROUS_TWO_PHASE* medium to mass transfer coefficient times driving force (see
*Media Type* card). The flux quantity is specified on a per mass basis so the mass
transfer coefficient is in units of L/t, and the sink density is in units of M/
:math:`L^3`.

The <float_list> for this boundary condition has four values; the definitions of the input parameters are as follows:

================ ==================================================================
**POROUS_FLUX**  Name of the boundary condition (<bc_name>).
**SS**           Type of boundary condition (<bc_type>), where **SS**
                 denotes side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain.
<integer>        Species number of transported species (currently only
                 used for multicomponent species in the phases, which as
                 of 11/2/01 is not active; so set to zero).
<float1>         Value of mass transfer coefficient, :math:`h_1` in units 
                 of L/t, consistent with gas phase concentration driving force.
<float2>         Driving force concentration in external phase, i.e., sink
                 density, :math:`p_{gt}^0` in units of M/ :math:`L^3`.
<float3>         Value of pressure-driven mass transfer coefficient, 
                 :math:`h_2` in
                 units of 1/L, for a liquid exiting a partially saturated
                 domain.
<float4>         Driving force concentration in external phase, i.e., sink
                 pressure for liquid extraction, :math:`p_{liq}^0` in units 
                 of M/L/t.
================ ==================================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = POROUS_FLUX SS 12 0 0.03 0. 0. 0. 0.

This card applies the lumped mass transfer model for the liquid phase solvent with a
mass transfer coefficient of 0.03 and a sink density of 0.0 for the total flux. The
boundary condition is applied to side set 12 and to the species number 0. This species
number is currently not used and ignored.

-------------------------
**Technical Discussion**
-------------------------

The mathematical form for this boundary condition is as follows

.. figure:: /figures/175_goma_physics.png
	:align: center
	:width: 90%

where the left hand side is the total flux of the liquid solvent i in the medium, which
includes, in order, the flux due to Darcy flow of gas vapor, the Darcy flow of liquid
solvent, the diffusive flux of gas vapor in the pore space and the diffusive flux of liquid
solvent in the liquid phase. The parameters are 
:math:`h_1`, :math:`p_g^0` , :math:`h_2` , :math:`p{liq}^0` and as defined on the
input card. :math:`v_s` is the user supplied convection velocity of the stress-free state as defined
on the *Convective Lagrangian Velocity* card.

At the present time (11/2/01), this condition is only used for single component liquid
phases and has not been furbished for multicomponent capability yet. Note that usually
the second term on the right is turned off, as in the example above, unless the liquid
pressure at the surface of the sample is greater than the external pressure. This term was
added for applications in which liquid is being squeezed out of a medium and then
drips off or disappears, as liquid is not allowed to be sucked back in (Heaviside
function, *H*), although the condition could be furbished for this.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk
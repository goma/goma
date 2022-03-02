***************
**CAP_REPULSE**
***************

::

	BC = CAP_REPULSE SS <bc_id> <float_list> [mat_id]

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This card functions much in the same way as the *CAPILLARY* card. It imposes surface
tension forces on free-surfaces. The addition to this card, however, is that in the vicinity
of a specified planar boundary, an additional repulsive force is added to the surface
tension force. This force is directed away from the planar surface and increases in
proportion to 1/ :math:`r^2` as the free-surface approaches the planar surface. This condition can
be used to contend with the difficult problem of fluid/solid contact in an approximate
way. This boundary condition is only applicable to two-dimensional problems; trying
to apply it in a three-dimensional problem will cause an error.

There are seven values in the <float_list>; definitions of the input parameters are as
follows:

================ ============================================================
**CAP_REPULSE**  Name of the boundary condition.
**SS**           Type of boundary condition (<bc_type>), where **SS**
                 denotes side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain.
<float1>         :math:`\sigma`, the surface tension value or a multiplier 
                 for the capillary effect. See the important caveat under 
                 the *CAPILLARY* card regarding the use of this parameter
                 when a surface tension value is supplied in the material
                 file.
<float2>         :math:`P_{ex}`, the applied external pressure field on the 
                 free surface.
<float3>         :math:`P_{rep}`, the coefficient on the surface repulsion term.
                 This parameter should have units of (ML/ :math:`T^2`). See 
                 below for an exact description of the surface repulsion term.
<float4>         *a*, the sensitivity with respect to x-coordinate (the *a*
                 coefficient) of the plane surface that is repelling the free
                 surface sideset.
<float5>         *b*, the sensitivity with respect to y-coordinate (the *b*
                 coefficient) of the plane surface that is repelling the free
                 surface sideset.
<float6>         *c*, the sensitivity with respect to z-coordinate (the *c*
                 coefficient) of the plane surface that is repelling the free
                 surface sideset.
<float7>         *d*, the constant *d* coefficient of the plane surface
                 equation that is repelling the free surface sideset.
[mat_id]         In the case of a surface node shared by more than one
                 material, this optional integer parameter allows the user
                 to specify which material the condition will be applied
                 in. This is rarely used.
================ ============================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = CAP_REPULSE SS 24 1.0 0.0 0.1 1. 1. 0. 2.

applies a standard capillary surface tension pressure jump condition to side set 24,
except as the free surface approaches the plane surface defined by a solution to the
equation *x + y = -2.0*.

-------------------------
**Technical Discussion**
-------------------------

* This boundary condition applies the following force term to the fluid momentum 
  equation:

.. figure:: /figures/102_goma_physics.png
	:align: center
	:width: 90%

which is almost identical to the force applied by the *CAPILLARY* card. The only
difference is the last term on the right in which *d* is the normal distance from a
given point on the free-surface side set and the planar surface defined by the
equation:

.. figure:: /figures/103_goma_physics.png
	:align: center
	:width: 90%




.. TODO - Lines 81 and 90 have photos that needs to be replaced with the real equation.
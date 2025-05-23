*************
**CAPILLARY**
*************

::

	BC = CAPILLARY SS <bc_id> <float_list> [integer]

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition card is used to apply capillary forces (surface tension) to the
momentum equation on a free-surface.

Definitions of the input parameters are as follows:

============= ==================================================================
**CAPILLARY** Name of the boundary condition.
**SS**        Type of boundary condition (<bc_type>), where **SS** denotes
              side set in the EXODUS II database.
<bc_id>       The boundary flag identifier, an integer associated with
              <bc_type> that identifies the boundary location (side set in
              EXODUS II) in the problem domain.
<float1>      :math:`\sigma`, surface tension or capillary term multiplier. 
              **IMPORTANT NOTE: if no Surface Tension card appears in the material
              file, this parameter is the surface tension value used here. If
              Surface Tension is set in the material file, however, this float
              value will multiply the surface tension value from the
              material file prior to it’s application in this boundary
              condition. Best practice is to set this parameter to 1.0 and
              set the surface tension value in the material file.**
<float2>      :math:`P_{ex}`, the external applied isotropic pressure on the free
              surface.
<float3>      :math:`P_r`, deprecated capability. Always set to zero.
[integer]     Optional integer value indicating the element block id from
              which to apply the boundary condition. This is used to force
              the capillary stresses to be applied from within a phase
              where the momentum equations are defined.
============= ==================================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = CAPILLARY SS 12   1.0 10.0 0.0

This card specifies that capillary forces be applied to the free surface on side set 12. If a
surface tension material parameter value or model is supplied, this is the surface
tension value used. If not, the surface tension value used is 1.0. An external isotropic
pressure of 10.0 is applied from the surrounding environment.

-------------------------
**Technical Discussion**
-------------------------

* One of the primary characteristics of a free-surface is the presence of surface
  tension-related forces. This boundary condition permits application of such forces.
  The forces on the fluid at the free-surface are set via the following relation:

.. figure:: /figures/101_goma_physics.png
	:align: center
	:width: 90%

where *n* is the outward normal to the surface, *T* is the fluid stress tensor, 
:math:`P_{ex}` is the
external applied pressure described above, *H* is the surface curvature defined as,
H = –:math:`\Delta_s` ⋅ n ⁄ 2, :math:`\sigma` is the surface tension, and :math:`\Delta_s`
is the surface divergence operator defined as :math:`\Delta_s^f` = (I – nn) ⋅ 
:math:`\Delta` f .

* Typical usage of this boundary condition is in conjunction with a *KINEMATIC*
  boundary condition. The latter enforces no penetration of fluid through a free
  surface by deforming the mesh and this boundary condition acts on the fluid
  momentum equation to enforce the capillary jump condition given above.

* No end of confusion results from use of this card primarily because of overloading
  the surface tension parameter. To reiterate, the value for surface tension that
  appears on this card is the actual (constant) value of surface tension that is used *if a
  surface tension model has NOT been specified explicitly in the material file*. If
  such a model has been identified, the surface tension parameter in the *CAPILLARY*
  card is a *multiplier* to the surface tension. The best practice is to simply always use
  1.0 for this parameter and set the surface tension in the material file.

* The optional (integer) element block ID corresponds to the material numbers given
  in the *Problem Description* section of the input file.




.. TODO - Line 66 has a photo that needs to be replaced with the real equation.
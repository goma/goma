*********************
**VELO_TANGENT_EDGE**
*********************

::

	BC = VELO_TANGENT_EDGE SS <bc_id1> <bc_id2> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(PCC-EDGE/ROTATED MOMENTUM)**

This boundary condition card is used to make the velocity component tangent to the
contact line in the plane of the web equal to the component of web velocity
(:math:`W_x`, :math:`W_y`, :math:`W_z`) along the contact line. This constraint replaces the tangential component
of the *MOMENTUM* equation along the contact line. It is used with the
*VELO_NORMAL_EDGE* condition to impose a wetting line model onto dynamic
contact lines in three-dimensions. The constraint is a rotated collocated condition.

Definitions of the input parameters are as follows:

===================== ================================================================
**VELO_TANGENT_EDGE** Name of the boundary condition.
**SS**                Type of boundary condition (<bc_type>), where **SS**
                      denotes side set in the EXODUS II database.
<bc_id1>              The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location (side set
                      in EXODUS II) of the primary side set defining the edge
                      geometry in the problem domain. When applied to
                      dynamic contact lines, this side set should correspond to
                      the free surface.
<bc_id2>              The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location (side set
                      in EXODUS II) of the secondary side set defining the
                      edge geometry in the problem domain. The boundary
                      condition is applied to the curve defined as the
                      intersection of this side set with the primary side set
                      When applied to dynamic contact lines, this side set
                      should correspond to the substrate.
<float1>              :math:`W_x`, x-component of the substrate (or web) velocity.
<float2>              :math:`W_y`, y-component of the substrate (or web) velocity.
<float3>              :math:`W_z`, z-component of the substrate (or web) velocity.
===================== ================================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = VELO_TANGENT_EDGE SS 5 4 -1.0 0.0 0.0

This card imposes a tangent velocity component along the curve formed by the
intersection of sidesets 5 and 4. The value of the component is the projection of the
substrate velocity (-1.0, 0. ,0.) into the tangent direction. The tangent direction is along
the curve itself.

-------------------------
**Technical Discussion**
-------------------------

* This equation imposes the following constraint as a point collocated condition at
  the integration points of the elements along the curve:

.. figure:: /figures/086_goma_physics.png
	:align: center
	:width: 90%

where :math:`t_{cl}` is a vector tangent to the curve, *v* is the fluid velocity, and *W* is the
(constant) velocity of the moving substrate. The reader is referred to the sketch
appearing with the *VELO_NORMAL_EDGE* card for a depiction of these vectors.
It is applied as a point collocated condition at the integration points of the line
elements along the curve.

* As noted above this boundary condition is used in concert with the
  *VELO_NORMAL_EDGE* condition to impose a model of wetting line physics
  along a dynamic contact line in three dimensions. The reader is referred to the
  discussion section of this latter boundary condition for a thorough exposition of
  this model. Suffice it to say that this boundary condition enforces no-slip between
  substrate and fluid in the tangent direction to the contact line. This is an essential
  part of the wetting line model because it implies that the wetting line forces related
  to surface tension etc. do not act tangential to the wetting line. Therefore, there is
  no agent in this direction which could account for departures from a strictly no-slip
  boundary condition.

* The astute user might note that the mesh velocity doesn’t appear in this expression
  whereas it does in the expression for *VELO_NORMAL_EDGE*. In the latter
  expression, the normal motion of the mesh represents the wetting velocity of the
  contact line normal to itself. It has a physical significance and so it make senses to
  connect it to the fluid velocity at that point. In the case of the tangential mesh
  motion velocity, it cannot be attached to any obvious physical part of the wetting
  model. It makes no sense that the tangential motion of nodes along the contact line
  should induce velocity in the fluid and vice versa. As a result, mesh motion is left
  out of the preceding relation.




.. TODO - In line 68, the photo needs to be repalced by the proper equation.
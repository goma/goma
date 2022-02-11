*****************
**CA_EDGE_CURVE**
*****************

::

	BC = CA_EDGE_CURVE SS <bc_id1> <bc_id2> <float1>

-----------------------
**Description / Usage**
-----------------------

**(PCC-EDGE/ROTATED MESH)**

This boundary condition allows the user to specify a constant contact angle along an
edge in three-dimensions. It is similar in function to the *CA_EDGE* boundary condition
in which the contact angle is enforced with respect to a fixed vector. However, for this
boundary condition, the contact angle is enforced with respect to the normal of the
*secondary* side set thereby permitting a contact angle constraint to be applied on a
curving surface. The boundary condition is applied to the edge curve defined by the
intersection of the primary and secondary side sets.

Definitions of the input parameters are as follows:

================== ================================================================
**CA_EDGE_CURVE**  Name of the boundary condition.
**SS**             Type of boundary condition (<bc_type>), where **SS** denotes
                   side set in the EXODUS II database.
<bc_id1>           The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (side set in
                   EXODUS II) in the problem domain. This identifies the
                   *primary side set*; in almost all cases it should also be a free
                   surface.
<bc_id2>           The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (side set in 
                   EXODUS II) in the problem domain. This identifies the
                   *secondary side set*. The outwards-pointing normal vector to
                   this side set is used as the *substrate* vector when enforcing
                   the contact angle constraint.
<float1>           the enforced contact angle, in degrees. Its value should lie in
                   the range 0 ≤ angle ≤ 180.
================== ================================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = CA_EDGE_CURVE SS 40 50 135.0

This boundary condition will enforce a 135 degree angle between the normal to the free
surface on side set 40 and the outward-pointing normal to side set 50 at all points along
the edge defined by side set 40 and 50. There is no restriction on whether side set 50’s
normal vectors must be constant.

-------------------------
**Technical Discussion**
-------------------------

* Although this boundary condition deals with vector quantities it is a scalar
  constraint. The actual requirement that is imposed is:

  .. figure:: /figures/059_goma_physics.png
	:align: center
	:width: 90%

  where n\ :sub:`f` is the outward-pointing normal to the
  primary side set, n\ :sub:`s` is the outward-pointing normal
  to the secondary side set, and θ is the angle
  supplied on the card. There is always some
  confusion regarding the sense of the angle; use the
  figure to the right for guidance. Note that the sense
  depicted here is at odds with the usual contact
  angle convention. Keep this in mind when using
  this card.

    .. figure:: /figures/060_goma_physics.png
	:align: center
	:width: 90%

* As in the case of the *CA_EDGE* condition, this condition is also a strongly
  enforced point collocated condition.

* Related boundary conditions: *CA_EDGE, CA_EDGE_INT,
  CA_EDGE_CURVE_INT, VAR_CA_EDGE, VAR_CA_USER*.



--------------
**References**
--------------

No References.
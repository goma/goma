***********
**CA_EDGE**
***********

::

	BC = CA_EDGE SS <bc_id1> <bc_id2> <float_list>

-----------------------
**Description / Usage**
-----------------------

.. warning::

  The angle in this BC is specified in degrees, not radians. This is 
  different than the 2D equivalent.

**(PCC-EDGE/ROTATED MESH)**

This boundary condition card specifies a constant contact angle on the edge defined by
the intersection of the primary and secondary side sets. This card is used most often to
enforce contact angle conditions on three-dimensional static contact lines. *It should not
be used in two-dimensional problems*, where the *CA* boundary condition is the
appropriate choice.

The contact angle supplied on the card will be enforced so that it is the angle between
the outward-pointing normal of the primary side set and the unit vector supplied on the
card. It is important to note that this outward-pointing normal should be variable, that is
to say, the primary side set is most likely a free-surface.

Definitions of the input parameters are as follows:

============== =================================================================
**CA_EDGE**    Name of the boundary condition.
**SS**         Type of boundary condition (<bc_type>), where **SS** denotes
               side set in the EXODUS II database.
<bc_id1>       The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set in
               EXODUS II) in the problem domain. This identifies the
               *primary side set*; in almost all cases it should also be a free
               surface.
<bc_id2>       The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set in
               EXODUS II) in the problem domain. This identifies the
               *secondary side set*, which plays no other role in this
               boundary condition than to provide a means of defining the 
               appropriate edge geometry in conjunction with the primary
               side set. Thus, the secondary side set will often represent a
               solid boundary.
<float1>       angle, value specifying the enforced angle, in degrees; it
               should lie in the range 0 ≤ angle ≤ 180.
<float2>       n\ :sub:`x`, the x-component of the fixed unit vector.
<float3>       n\ :sub:`y`, the y-component of the fixed unit vector.
<float4>       n\ :sub:`z`, the z-component of the fixed unit vector.
============== =================================================================

This boundary condition is a point collocated condition so it will be enforced exactly at
every node that lies on the edge (subject to overriding *ROT* cards or Dirichlet
conditions).

------------
**Examples**
------------

The following is a sample input card:
::

     BC = CA_EDGE SS 40 50 33.0 0. 1. 0.

This card will result in an angle of 33 degrees between the outward-pointing normal to
side set 40 and the vector (0,1,0) at all points on the edge defined by the intersection of
side set 40 and side set 50.

-------------------------
**Technical Discussion**
-------------------------

* Although this constraint deals with vector quantities, it is a scalar constraint. The
  actual requirement that is imposed is:

  .. math::
  
    n_f \cdot n = cos\ \left(\theta \right)
  

  where n\ :sub:`f` is the outward-pointing normal to the primary side set, *n* is the vector
  supplied on the card, and θ is the angle supplied on the card. It should be
  recognized that there are usually two orientations for n\ :sub:`f` which would satisfy this
  constraint. Most often the surrounding physics will choose the correct one, but
  there is nothing to guarantee this in special situations, for example, values for θ
  near zero or near 180.

* This boundary condition is a point collocated condition so the preceding
  constraint, will be enforce exactly and strongly for each node on the edge. The
  actual free surface normal is an average of vectors supplied by adjacent elements
  sharing a given node.

* As noted above, this boundary condition is most often used in three-dimensional
  free surface problems to enforce static contact angle conditions at the junction of a
  free, capillary surface and a solid boundary. The normal vector supplied on the
  card would be the normal to this solid boundary. Since this vector is a constant,
  there is the restriction that in this application this boundary condition can only be
  used to specify a contact angle with respect to a *planar* solid boundary. A different
  boundary condition, *CA_EDGE_CURVE*, should be used if the solid boundary is
  not planar.

* Related boundary conditions: *CA_EDGE_INT, CA_EDGE_CURVE,
  CA_EDGE_CURVE_INT, VAR_CA_EDGE, VAR_CA_USER*.



--------------
**References**
--------------

No References.
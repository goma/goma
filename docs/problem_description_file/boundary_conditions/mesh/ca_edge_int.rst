***************
**CA_EDGE_INT**
***************

::

	BC = CA_EDGE_INT SS <bc_id1> <bc_id2> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC-EDGE/ROTATED MESH)**

This boundary condition card specifies a constant contact angle on the edge defined by
the intersection of the primary and secondary side sets. It is identical in format and
function as the *CA_EDGE* boundary condition. The only difference is that this
boundary condition is a strong integrated constraint.

Definitions of the input parameters are as follows:

================ =================================================================
**CA_EDGE_INT**  Name of the boundary condition.
**SS**           Type of boundary condition (<bc_type>), where **SS** denotes
                 side set in the EXODUS II database.
<bc_id1>         The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set in
                 EXODUS II) in the problem domain. This identifies the
                 *primary side set*; in almost all cases it should also be a free
                 surface.
<bc_id2>         The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set in
                 EXODUS II) in the problem domain. This identifies the
                 *secondary side set*, which plays no other role in this
                 boundary condition than to provide a means of defining the
                 appropriate edge geometry in conjunction with the primary
                 side set. Thus, the secondary side set will often represent a
                 solid boundary.
<float1>         angle, value specifying the enforced angle, in degrees; it
                 should lie in the range 0 ≤ angle ≤ 180.
<float2>         n\ :sub:`x`, the x-component of the fixed unit vector.
<float3>         n\ :sub:`y`, the y-component of the fixed unit vector.
<float4>         n\ :sub:`z`, the z-component of the fixed unit vector.
================ =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = CA_EDGE_INT SS 40 50 33.0 0. 1. 0.

This card will result in an angle of 33 degrees between the outward-pointing normal to
side set 40 and the vector (0,1,0) at all points on the edge defined by the intersection of
side set 40 and side set 50.

-------------------------
**Technical Discussion**
-------------------------

* As noted above, this boundary condition is identical in function to the *CA_EDGE*
  condition. It differs only in the manner of its application. Whereas, the former was
  a point collocated constraint, this boundary condition strongly enforces the
  following *integrated* constraint at a node i**:

.. math::

  \int_{\Gamma} \phi_i \left(n_f \cdot n - cos\ (\theta) \right) d \Gamma = 0

  

|

  where φ\ :sub:`i` is the finite element trial function for node i, Γ is the edge space curve, n\ :sub:
  `f`
  is the outward-pointing normal to the primary sideset, n is the vector supplied on
  the card, and θ is the angle supplied on the card. Because it is an integrated
  constraint, evaluation of the free-surface normal vector is done at integration
  points between nodes on the edge. Therefore, there is no averaging of normal
  vectors. This is sometimes advantageous when there are discontinuities in the
  slope of the edge curve.

* Related boundary conditions: *CA_EDGE, CA_EDGE_CURVE,
  CA_EDGE_CURVE_INT, VAR_CA_EDGE, VAR_CA_USER*.



--------------
**References**
--------------

No References.
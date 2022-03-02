******************
**CA_EDGE_OR_FIX**
******************

::

	BC = CA_EDGE_OR_FIX SS <bc_id1> <bc_id2> <type_string> {float_list}

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MESH)**

In analogy to the two-dimensional condition, *CA_OR_FIX*, boundary condition, this
boundary condition imposes a contact angle on an edge feature in a three-dimensional
mesh. However, this condition also permits the user to specify a closed curve on the
substrate plane on which the contact line will attach and not move past. This permits
modeling of geometric features in which the substrate slope is discontinuous. When
contact lines encounter such sharp features, usually they arrest. The boundary
condition also permits the contact line to release from the curve if the overall fluid
mechanics would promote a recession of the contact line.

Description of the card parameters is as follows:

+-------------------+------------------------------------------------------------------+
|**CA_EDGE_OR_FIX** | Name of boundary condition.                                      |
+-------------------+------------------------------------------------------------------+
|**SS**             | Type of boundary condition (<bc_type>), where **SS** denotes     |
|                   | side set in the EXODUS II database.                              |
+-------------------+------------------------------------------------------------------+
|<bc_id1>           | The boundary flag identifier, an integer associated with         |
|                   | <bc_type> that identifies the boundary location (side set in     |
|                   | EXODUS II) in the problem domain. This identifies the            |
|                   | *primary side set* defining the edge curve on which this         |
|                   | condition applies.                                               |
+-------------------+------------------------------------------------------------------+
|<bc_id2>           | The boundary flag identifier, an integer associated with         |
|                   | <bc_type> that identifies the boundary location (side set in     |
|                   | EXODUS II) in the problem domain. This identifies the            |
|                   | *secondary side set* defining the edge curve on which this       |
|                   | condition applies. Taken together, the edge curve is the         |
|                   | intersection of the primary and secondary sidesets.              |
+-------------------+------------------------------------------------------------------+
|<type_string>      | A string identifying the type of *feature curve* being defined;  |
|                   | currently, there are only two choices: **CIRCLE** and **USER**.  |
|                   | The **CIRCLE** options indicates that the surface feature on     |
|                   | which a Gibb’s criterion is applied is a circle in the substrate |
|                   | plane. The **USER** option indicates that the user will have to  |
|                   | provide a geometric definition in the user subroutine            |
|                   | user_gibbs_criterion in the file user_bc.c.                      |
+-------------------+------------------------------------------------------------------+
|{float_list}       | A list of float parameters to be used in defining the contact    | 
|                   | angle, the normal to the substrate, and other geometric          |
|                   | parameters used to define the feature curve. For each            |
|                   | <type_string> choice there is a different set of float           |
|                   | parameters:                                                      |
+-------------------+------------------------------------------------------------------+
|                   | **CIRCLE**   <float_list>                                        |
+-------------------+------------+-----------------------------------------------------+
|                   | <float1>   | θ\ :sub:`dcl`, contact angle at dynamic contact     |
|                   |            | line, in radians                                    |
+-------------------+------------+-----------------------------------------------------+
|                   | <float2>   | n\ :sub:`x`, x-component of outward substrate normal|
+-------------------+------------+-----------------------------------------------------+
|                   | <float3>   | n\ :sub:`y`, y-component of outward substrate normal|
+-------------------+------------+-----------------------------------------------------+
|                   | <float4>   | n\ :sub:`z`, z-component of outward substrate normal|
+-------------------+------------+-----------------------------------------------------+
|                   | <float5>   | c\ :sub:`x`, x coordinate of circle center          |
+-------------------+------------+-----------------------------------------------------+
|                   | <float6>   | c\ :sub:`y`, y-coordinate of circle center          |
+-------------------+------------+-----------------------------------------------------+
|                   | <float7>   | c\ :sub:`z`, z-coordinate of circle center          |
+-------------------+------------+-----------------------------------------------------+
|                   | <float8>   | r, radius of circle                                 |
+-------------------+------------+-----------------------------------------------------+
|                   | The sign of this last parameter is important. If negative, the   |
|                   | implication is that the starting location of the contact line is |
|                   | outside of the circle. If positive, the original location is     |
|                   | assumed to be completely inside the circle.                      |
+-------------------+------------------------------------------------------------------+
|                   | **USER**   <float_list>                                          |
+-------------------+------------+-----------------------------------------------------+
|                   | <floati>   | a list of float values that are passed to the       |
|                   |            | function user_gibbs_criterion in                    |
|                   |            | the one-dimensional array *p* in the order in       |
|                   |            | which they appear on the card from left to          |
|                   |            | right. The user must be certain that the            |
|                   |            | parameters appearing here are sufficient            |
|                   |            | for applying the Gibbs criterion as well as         |
|                   |            | imposing the appropriate contact angle.             |
+-------------------+------------+-----------------------------------------------------+

------------
**Examples**
------------

An example making use of the **CIRCLE** feature curve option is as follows:
::

     BC = CA_EDGE_OR_FIX SS 10 20 CIRCLE 1.3   0. -1. 0.   0. 0. 0.   1.0

This card applies to the intersection between side sets 10 and 20. The constant contact
angle applied is 1.3 radians. The substrate outward normal is (0, -1, 0). The feature is a
circle of radius 1.0 centered at (0.0, 0.0, 0.0). The original location for the contact line
must be completely inside of the feature circle. Note also that the circle center should
lie in the substrate plane.

-------------------------
**Technical Discussion**
-------------------------

* See the Technical Discussion under the boundary condition *CA_OR_FIX* for a
  detailed discussion of the nature of the Gibb’s criterion as it applies to contact
  lines. In a nutshell, however, the basic notion is that the contact line is free to
  advance over the substrate with an imposed contact angle, constant or dependent
  on the local conditions. When the contact angle encounters the geometric feature
  defined in the function user_gibbs_criterion, it is captured at that point
  and no longer advances. The contact angle is allowed to vary as long as it is held at
  the feature. The boundary condition also permits the contact line to release from
  the feature curve and recede the way it came if the contact angle ever becomes
  larger than its mobile value.

* So the phenomena that can be modeled with this boundary condition are those in
  which a contact line moves to, for example, the edge of cylinder. At the edge, the
  very small curvature of this feature effectively presents a barrier to further advance
  of the contact line provided the deformation of the free surface beyond the vertical
  boundaries of the cylinder is not too large. In the fullness of time, it might also be
  the case that the free surface is drawn backwards in the direction of the cylinder
  axis. The contact line should also recede and this boundary condition permits this
  once the contact angle it makes with the cylinder top exceeds the mobile contact
  angle by a small amount.



--------------
**References**
--------------

No References.
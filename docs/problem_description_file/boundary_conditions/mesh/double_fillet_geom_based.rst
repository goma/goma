***************
DOUBLE_FILLET
***************

::

	BC = DOUBLE_FILLET_GEOM_BASED SS <bc_id> xpt1 ypt1 theta1 r1 xpt2 ypt2 theta2 r2 curv_mid

-----------------------
**Description / Usage**
-----------------------

This is similar to the `DOUBLE_FILLET` condition but is implemented based on
geometric partitions of the die rather than angles.

**(PCC/ROTATED MESH)**

Definitions of the input parameters are as follows:

DOUBLE_FILLET   
   Name of the boundary condition <bc_name>).
SS           
   Type of boundary condition (<bc_type>), where **SS** denotes
   side set in the EXODUS II database.
<bc_id>
   The boundary flag identifier, an integer associated with
   <bc_type> that identifies the boundary location (side set in
   EXODUS II) in the problem domain.
xpt1, ypt1, theta1, r1
   - xpt1, ypt1: Coordinates of the first point
   - theta1: Angle at the first point
   - r1: Radius at the first point
xpt2, ypt2, theta2, r2
   - xpt2, ypt2: Coordinates of the second point
   - theta2: Angle at the second point
   - r2: Radius at the second point
curv_mid
   Middle surface curvature value

------------
**Examples**
------------

The following sample input card:
::

     BC   = DOUBLE_FILLET_GEOM_BASED SS 56 {xpt1} {ypt1} {theta1} {r1} {xpt2} {ypt2} {theta2} {r2} {curv_mid}



-------------------------
**Technical Discussion**
-------------------------

This condition, like *DISTNG, PLANE*, and others that can be applied to geometry, is
applied to the normal component of the mesh motion equations along a boundary in
two dimensions; in three dimensions application needs to be further directed with the
*ROT* conditions. Examples of typical distinguishing conditions can be found in
user_bc.c in the fnc routine and companion derivative routines.

This is strictly in the x-y plane.

.. figure:: /figures/double_fillet_diagram.png
	:align: center
	:width: 90%

	Diagram of the double fillet boundary condition, showing the coordinates and
	parameters defining the two fillet points and their potential applications
	in a slot coater.  In this case there would be a `DOUBLE_FILLET` for both
	upstream and downstream die lips.





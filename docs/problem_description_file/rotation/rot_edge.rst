************
**ROT EDGE**
************

::

	ROT = {MESH | MOM} EDGE <bc_id1> <bc_id2> <string_x> <int_x> <string_y>
    <int_y> <string_z> <int_z> {seed_method} <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

This rotation specification card deals with rotation specification along edges. In this context, an edge is the intersection of two side sets. It identifies the boundary conditions that will be applied at nodes on that edge and which equation components are to be associated with them. It also identifies which components of the rotated equations will be used. Currently, only rotation of mesh and momentum equations is allowed on edges. This card can also be used for specifying a seed vector if needed.

Definitions of the input parameters are as follows:

+----------------+--------------------------------------------------+
|**{MESH | MOM}**|Equation type to which this rotation condition    |
|                |applies:                                          |
|                |                                                  |
|                | * **MESH** - Applies to mesh equations.          |
|                | * **MOM** - Applies to fluid momentum equations. |
+----------------+--------------------------------------------------+
|**EDGE**        |Type of rotation specification.                   |
+----------------+--------------------------------------------------+
|<bc_id1>        |Side set ID number of the primary side set.       |
+----------------+--------------------------------------------------+
|<bc_id2>        |Side set ID number of the secondary side set.     |
+----------------+--------------------------------------------------+

The edge is defined as the intersection of the primary and secondary side sets.

The next six parameters dictate how the x, y, and z components of the vector equation are replaced by boundary conditions or rotated components using pairs of specifiers, e.g., <string_x> and <int_x> for the x-component of the equation.

+----------------+-------------------------------------------------------+
|<string_x>      |A character string that specifies what will replace the|
|                |x-component of the vector equation (MESH or            |
|                |MOMENTUM). This string may be the name of a            |
|                |boundary condition already specified in the boundary   |
|                |condition specification section or one of the rotation |
|                |strings listed in Table 3 (Valid EDGE Tangent Equation |
|                |Rotation Strings).                                     |
+----------------+-------------------------------------------------------+
|<int_x>         |This is an integer parameter specified as follows:     |
|                |                                                       |
|                | * If <string_x> is a boundary condition name, then    |
|                |   <int_x> is the side set or node set designation     |
|                |   to which the appropriate boundary condition applies.|
|                |   This provides a means of distinguishing between     |
|                |   boundary conditions possessing the same string name |
|                |   but applied to different side sets or node sets.    |
|                | * If <string_x> is a rotation string from Table 3,    |
|                |   <int_x> should be specified as 0.                   |
+----------------+-------------------------------------------------------+
|<string_y>      |A character string that specifies what will replace the|
|                |y-component of the vector equation (MESH or MOMENTUM). |
|                |This string may be the name of a boundary condition    |
|                |already specified in the boundary condition            |
|                |specification section or one of the rotation strings   |
|                |listed in Table 3.                                     |
+----------------+-------------------------------------------------------+
|<int_y>         |This is an integer parameter specified as follows:     |
|                |                                                       |
|                | * If <string_y> is a boundary condition name, then    |
|                |   <int_y> is the side set or node set designation to  |
|                |   which the appropriate boundary condition applies.   |
|                |   This provides a means of distinguishing between     |
|                |   boundary conditions possessing the same string name |
|                |   but applied to different side sets or node sets.    |
|                | * If <string_y> is a rotation string from Table 3,    |
|                |   <int_y> should be specified as 0.                   |
+----------------+-------------------------------------------------------+
|<string_z>      |A character string that specifies what will replace the|
|                |z-component of the vector equation (MESH or MOMENTUM). |
|                |This string may be the name of a boundary condition    |
|                |already specified in the boundary condition            |
|                |specification section or one of the rotation strings   |
|                |listed in Table 3.                                     |
+----------------+-------------------------------------------------------+
|<int_z>         |This is an integer parameter specified as follows:     |
|                |                                                       |
|                | * If <string_z> is a boundary condition name, then    |
|                |   <int_z> is the side set or node set designation to  |
|                |   which the appropriate boundary condition applies.   |
|                |   This provides a means of distinguishing between     |
|                |   boundary conditions possessing the same string name |
|                |   but applied to different side sets or node sets.    |
|                | * If <string_z> is a rotation string from Table 3,    | 
|                |   <int_z> should be specified as 0.                   |
+----------------+-------------------------------------------------------+

**Table 3. Valid EDGE Tangent Equation Rotation Strings**

======================= =====================================================
**{string_x|y|z}**      **Description of Equation Rotation Selections**
======================= =====================================================
**NONE, NA,** or **NO** No rotation is performed for this equation component.
**N**                   This equation component is replaced by the normal
                        component of the residual: n • R where n is the outwardpointing normal to <bc_id1>
**T**                   This equation component is replaced by the tangential
                        component of the residual: T • R where the tangent is a line tangent along the edge defined by <bc_id1> and
                        <bc_id2>.
**B**                   This equation component is replaced by the 
                        outwardpointing binormal component of the residual:
                        B • R where the binormal is perpendicular to both the line tangent T and the outward-pointing normal to <bc_id1>.
**S**                   The equation component is replaced by the projection 
                        of the equations in the direction of the seed vector:
                        S • R
**X**                   This equation is replaced by the x-component of the
                        residual.
**Y**                   This equation is replaced by the y-component of the
                        residual.
**Z**                   This equation is replaced by the z-component of the
                        residual.
======================= =====================================================

In most cases, seeding of the tangent vectors is not needed along edges, although it is possible to specify a seed method as defined in the *ROT SURFACE* card via the parameters {seed_method}, <float1>, <float2>, and <float3>. Note also that a seed vector must be specified to use the S rotation option.

------------
**Examples**
------------

The following is an example of an edge rotation specification:
::

   ROT = MESH   EDGE 4 5   PLANE   4   PLANE 5   T   0   NONE

This card specifies rotation of the mesh equations along the edge of intersection of side sets 4 and 5. The x and y mesh equations are replaced by PLANE conditions on side sets 4 and 5, respectively. The z mesh equation is replaced by the mesh residuals rotated into the direction of the line tangent along the edge. This enables the mesh to slide freely (i.e., stress-free) along the edge.

-------------------------
**Technical Discussion**
-------------------------

* The direction of the line tangent is chosen such that the binormal 
  ( b = n × t ) with *n*, the outward-pointing normal to the primary surface 
  <bc_id1>, is outwardpointing from the edge.

* Along edges, two of the equations are normally replaced by boundary
  conditions and one equation is replaced by this tangential component. However several options are available for replacing the mesh equations by other forms of the rotated equations as listed in Table 3. (Valid EDGE Tangent Equation Rotation Strings) above.

* It is very rare to require a seed vector be specified on an edge. The 
  SEED vector choice is almost always NONE.

* A precedence rule has been established for the case when more than one 
  *Rotation Specification* could be applied at a point. The rule is as follows:

  The Rotation condition that will be applied is:

       *A>The first VERTEX condition in the input deck that could
       apply. If there is no contravening VERTEX condition then,*

       *B>The first EDGE condition in the input deck that could
       apply. If there is no contravening EDGE condition then,*

       *C>The first SURFACE condition that could apply.*

* A very important restriction exists for EDGE and VERTEX rotation conditions. It is a necessary requirement that all elements that are present on an edge have only a single segment present on the edge curve. An element may therefore never contribute more than two corner vertex nodes to the set of nodes on an edge curve. If there are more than two such nodes for a given element, *Goma* will terminate with a *“Side not connected to edge”* error. If such a situation exists, the only solution is to remesh the geometry to eliminate such elements.



--------------
**References**
--------------

GT-007.2: Tutorial on droplet on incline problem, July 30, 1999, T. A. Baer

GT-012.0: 3D Roll coating template and tutorial for GOMA, February 21, 2000, P.R. Schunk

GT-018.1: ROT card tutorial, January 22, 2001, T. A. Baer
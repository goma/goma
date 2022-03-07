**************
**ROT VERTEX**
**************

::

	ROT = {MESH | MOM} VERTEX <bc_id1> <bc_id2> <bc_id3> <string_x>
    <int_x> <string_y> <int_y> <string_z> <int_z> {seed_method} <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

This rotation specification card deals with rotation specification at vertices. In this context, a vertex is the intersection point of three side sets. It identifies the boundary conditions that will be applied at the vertex node and which equation components are to be associated with them. It also identifies which components of the rotated equations will be used. Currently, only rotation of mesh and momentum equations is allowed at a vertex. This card can also be used for specifying a seed vector if needed.

Definitions of the input parameters are as follows:

+----------------+--------------------------------------------------+
|**{MESH | MOM}**|Type of equation to which this specification      | 
|                |applies,where                                     |
|                |                                                  |
|                | * **MESH** - Applies to mesh displacement        |
|                |   equations.                                     |
|                | * **MOM** - Applies to fluid momentum equations. |
+----------------+--------------------------------------------------+
|**VERTEX**      |Type of rotation specification.                   |
+----------------+--------------------------------------------------+
|<bc_id1>        |Side set ID number of the primary side set.       |
+----------------+--------------------------------------------------+
|<bc_id2>        |Side set ID number of the secondary side set.     |
+----------------+--------------------------------------------------+
|<bc_id3>        |Side set ID number of the tertiary side set.      |
+----------------+--------------------------------------------------+

The vertex is defined as the point at the intersection of the primary, secondary, tertiary side set. Note that it is possible for these three side sets to intersect at more than one discrete point. The VERTEX condition is applied to all such points.

The next six parameters dictate how the x, y, and z components of the vector equation are replaced by boundary conditions or rotated components using pairs of specifiers, e.g., <string_x> and <int_x> for the x-component of the equation.

+----------------+-------------------------------------------------------+
|<string_x>      |A character string that specifies what will replace the|
|                |x-component of the vector equation (MESH or            |
|                |MOMENTUM). This string may be the name of a            |
|                |boundary condition already specified in the boundary   |
|                |condition specification section or one of the rotation |
|                |strings listed in Table 4 (Valid VERTEX Tangent        |
|                |Equation Rotation Strings).                            |
+----------------+-------------------------------------------------------+
|<int_x>         |This is an integer parameter specified as follows:     |
|                |                                                       |
|                | * If <string_x> is a boundary condition name, then    |
|                |   <int_x> is the side set or node set designation     |
|                |   to which the appropriate boundary condition applies.|
|                |   This provides a means of distinguishing between     |
|                |   boundary conditions possessing the same string name |
|                |   but applied to different side sets or node sets.    |
|                | * If <string_x> is a rotation string from Table 4,    |
|                |   <int_x> should be specified as 0.                   |
+----------------+-------------------------------------------------------+
|<string_y>      |A character string that specifies what will replace the|
|                |y-component of the vector equation (MESH or MOMENTUM). |
|                |This string may be the name of a boundary condition    |
|                |already specified in the boundary condition            |
|                |specification section or one of the rotation strings   |
|                |listed in Table 4.                                     |
+----------------+-------------------------------------------------------+
|<int_y>         |This is an integer parameter specified as follows:     |
|                |                                                       |
|                | * If <string_y> is a boundary condition name, then    |
|                |   <int_y> is the side set or node set designation to  |
|                |   which the appropriate boundary condition applies.   |
|                |   This provides a means of distinguishing between     |
|                |   boundary conditions possessing the same string name |
|                |   but applied to different side sets or node sets.    |
|                | * If <string_y> is a rotation string from Table 4,    |
|                |   <int_y> should be specified as 0.                   |
+----------------+-------------------------------------------------------+
|<string_z>      |A character string that specifies what will replace the|
|                |z-component of the vector equation (MESH or MOMENTUM). |
|                |This string may be the name of a boundary condition    |
|                |already specified in the boundary condition            |
|                |specification section or one of the rotation strings   |
|                |listed in Table 4.                                     |
+----------------+-------------------------------------------------------+
|<int_z>         |This is an integer parameter specified as follows:     |
|                |                                                       |
|                | * If <string_z> is a boundary condition name, then    |
|                |   <int_z> is the side set or node set designation to  |
|                |   which the appropriate boundary condition applies.   |
|                |   This provides a means of distinguishing between     |
|                |   boundary conditions possessing the same string name |
|                |   but applied to different side sets or node sets.    |
|                | * If <string_z> is a rotation string from Table 4,    | 
|                |   <int_z> should be specified as 0.                   |
+----------------+-------------------------------------------------------+

**Table 4. Valid VERTEX Tangent Equation Rotation Strings**

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

In most cases, seeding of the tangent vectors is not needed along vertices, although it is possible to specify a seed method as defined in the *ROT SURFACE* card via the parameters {seed_method}, <float1>, <float2>, and <float3>. Note also that a seed vector must be specified to use the S rotation option.

------------
**Examples**
------------

The following is an example of a vertex rotation specification:
::

   ROT = MESH   VERTEX   3 4 6   PLANE   4 PLANE 3   PLANE   6   NONE

::

   ROT = MESH   VERTEX   5 4 6   PLANE   4 KINEMATIC 5   PLANE   6   NONE

In the first example, the vertex is at the intersection of side sets 3, 4 and 6, and the three mesh equations at this vertex are replaced by *PLANE* conditions from side sets 4, 3, and 6, respectively. In the second example, the vertex is at the intersection of side sets 4, 5, and 6, respectively. Since it is conceivable that side set 5 might represent a free surface that curves in three dimensions, the last *VERTEX* card might apply to more than one point.

-------------------------
**Technical Discussion**
-------------------------

* Despite the fact that *VERTEX* cards apply only at single points, 
  definitions  
  of the normal, tangent and binormal vectors are still operative. The normal vector, **N**, is the outward-pointing normal to the primary side set, the tangent vector, **T**, is defined to lie along the curve defined by the intersection of the primary and secondary side set, and the binormal vector, **B**, is defined simply as the cross product of the normal vector with the tangent vector. Note that the sense of the tangent vector is chosen so that the binormal vector will always point outwards from the domain.

* At a vertex, it is normally the case that all three rotated components will 
  be replaced by boundary conditions as suggested by the examples. However, it is not a rarity that a rotated component, usually **N** or **T**, will also appear.

* The same hierarchy of precedence is used to determine which rotation
  specification will be applied when more than one could apply to a node. The rule is as follows:

  The Rotation condition that will be applied is:

     *A>The first VERTEX condition in the input deck that could
     apply. If there is no contravening VERTEX condition then,*

     *B>The first EDGE condition in the input deck that could
     apply. If there is no contravening EDGE condition then,*

     *C>The first SURFACE condition that could apply.*

* Very often *VERTEX* cards are used to resolve ambiguities that arise at 
  points where multiple *SURFACE* or *EDGE* cards could apply.



--------------
**References**
--------------

GT-007.2: Tutorial on droplet on incline problem, July 30, 1999, T. A. Baer

GT-012.0: 3D Roll coating template and tutorial for GOMA, February 21, 2000, P.R. Schunk

GT-018.1: ROT card tutorial, January 22, 2001, T. A. Baer
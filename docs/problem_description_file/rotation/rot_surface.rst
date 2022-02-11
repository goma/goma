***************
**ROT SURFACE**
***************

::

	ROT = {MESH | MOM} SURFACE <bc_id> <string_x> <int_x> <string_y>
    <int_y> <string_z> <int_z> {seed_method} <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

This rotation specification card identifies a specific surface that requires rotation of equations for proper application of boundary conditions in three dimensional problems. It identifies the boundary conditions that are to be applied on that surface. It also identifies which equation components are to be replaced by boundary conditions, which are to be replaced by rotated equation components, and which are to be left alone. Equation components refer, currently, to only rotation of mesh and momentum equations. This card also identifies the manner in which two independent tangent vectors are to be determined on the surface.

Definitions of the first three input parameters are as follows:

+----------------+--------------------------------------------------+
|**{MESH | MOM}**|Equation type (<eq_type>) to which this rotation  |
|                |condition applies:                                |
|                |                                                  |
|                | * **MESH** - Applies to mesh equations.          |
|                | * **MOM** - Applies to fluid momentum equations. |
+----------------+--------------------------------------------------+
|**SURFACE**     |Type of rotation specification.                   |
+----------------+--------------------------------------------------+
|<bc_id>         |An integer identifying the side set designation   |
|                |of the surface to which this rotation condition   |
|                |applies.                                          |
+----------------+--------------------------------------------------+

The next six parameters dictate how the x, y, and z components of the vector equation are replaced by boundary conditions or rotated equations using pairs of specifiers, e.g., <string_x> and <int_x> for the x-component of the equation.

+----------------+-------------------------------------------------------+
|<string_x>      |A character string that specifies what will replace the|
|                |x-component of the vector equation (MESH or MOMENTUM). |
|                |This string may be the name of a boundary condition    |
|                |already specified in the boundary condition            |
|                |specification section or one of the rotation strings   |
|                |listed in Table 1 (Valid Equation Rotation Strings).   |
+----------------+-------------------------------------------------------+
|<int_x>         |This is an integer parameter specified as follows:     |
|                |                                                       |
|                | * If <string_x> is a boundary condition name, then    |
|                |   <int_x> is the side set or node set designation     |
|                |   to which the appropriate boundary condition applies.|
|                |   This provides a means of distinguishing between     |
|                |   boundary conditions possessing the same string name |
|                |   but applied to different side sets or node sets.    |
|                | * If <string_x> is a rotation string from Table 1,    |
|                |   <int_x> should be specified as 0.                   |
+----------------+-------------------------------------------------------+
|<string_y>      |A character string that specifies what will replace the|
|                |y-component of the vector equation (MESH or MOMENTUM). |
|                |This string may be the name of a boundary condition    |
|                |already specified in the boundary condition            |
|                |specification section or one of the rotation strings   |
|                |listed in Table 1.                                     |
+----------------+-------------------------------------------------------+
|<int_y>         |This is an integer parameter specified as follows:     |
|                |                                                       |
|                | * If <string_y> is a boundary condition name, then    |
|                |   <int_y> is the side set or node set designation to  |
|                |   which the appropriate boundary condition applies.   |
|                |   This provides a means of distinguishing between     |
|                |   boundary conditions possessing the same string name |
|                |   but applied to different side sets or node sets.    |
|                | * If <string_y> is a rotation string from Table 1,    |
|                |   <int_y> should be specified as 0.                   |
+----------------+-------------------------------------------------------+
|<string_z>      |A character string that specifies what will replace the|
|                |z-component of the vector equation (MESH or MOMENTUM). |
|                |This string may be the name of a boundary condition    |
|                |already specified in the boundary condition            |
|                |specification section or one of the rotation strings   |
|                |listed in Table 1.                                     |
+----------------+-------------------------------------------------------+
|<int_z>         |This is an integer parameter specified as follows:     |
|                |                                                       |
|                | * If <string_z> is a boundary condition name, then    |
|                |   <int_z> is the side set or node set designation to  |
|                |   which the appropriate boundary condition applies.   |
|                |   This provides a means of distinguishing between     |
|                |   boundary conditions possessing the same string name |
|                |   but applied to different side sets or node sets.    |
|                | * If <string_z> is a rotation string from Table 1,    | 
|                |   <int_z> should be specified as 0.                   |
+----------------+-------------------------------------------------------+

**Table 1. Valid Equation Rotation Strings**

======================= =====================================================
**{string_x|y|z}**      **Description of Equation Rotation Selections**
======================= =====================================================
**NONE, NA,** or **NO** No rotation is performed for this equation component.
**N**                   This equation component is replaced by normal
                        component of the residual tn • R
**T**                   This equation component is replaced by the tangential
                        component of the residual: t • R (EDGE and VERTEX only)
**T1**                  This equation component is replaced by the first
                        tangential component of the residual: T1 • R
**T2**                  This equation component is replaced by the second
                        tangential component of the residual: T2 • R
**X**                   This equation is replaced by the x-component of the
                        residual.
**Y**                   This equation is replaced by the y-component of the
                        residual.
**Z**                   This equation is replaced by the z-component of the
                        residual.
**S**                   The equation component is replaced by the projection 
                        of the equations in the direction of the seed vector:
                        S • R
**B**                   The equation component is replaced by the projection  
                        of the equations in the direction of the binormal vector.
======================= =====================================================

In most cases, only one of the three equations on a surface will be replaced by boundary conditions, and the remaining two equations will be rotated in the two tangent directions. Such a form constrains the normal motion of the solid or fluid while allowing tangential motions to occur stress-free.

The last four parameters in the card specify how to calculate the tangent vectors on the surface. In 3D, an infinite number of equally valid tangent pairs exist, so this card enables specifying how to choose those pairs. More specifically it identifies how to identify the first tangent vector (T1) since the second tangent vector is always be obtained via the cross product of the normal vector with the first tangent vector (T1).

============= ==========================================================
{seed_method} A character string that defines the method of tangent
              calculation. Valid options are listed in the Surface
              Tangent Calculation Method (Table 2).
<float1>      x-component of the seed vector, s. This parameter is
              only needed if {seed_method} is **SEED**.
<float2>      y-component of the seed vector, s. This parameter is
              only needed if {seed_method} is **SEED**.
<float3>      z-component of the seed vector, s. This parameter is
              only needed if {seed_method} is **SEED**.
============= ==========================================================

Note that the seed vector specified does not have to be a unit vector.

**Table 2. Surface Tangent Calculation Method**

======================= =====================================================
**{seed_method}**       **Description of Tangent Calculation Methods**
======================= =====================================================
**NONE**                Tangent vectors should not be calculated. This is the
                        usual choice for EDGE and VERTEX rotation types.
**SEED**                The first tangent vector (T1) is calculated from a 
                        surface projection of a seed vector, s: T1 = (I – nn) • s
**BASIS**               The first tangent is the direction of the first basis 
                        vector in the surface using a weighted average for adjacent elements.
**BASIS_FIRST**         The first tangent is the direction of the first  
                        finite element basis vector in the first element containing a given node.
**BASIS_RESEED**        The tangent resulting from BASIS_FIRST is used to
                        reseed tangent calculation in the adjacent elements. 
                        (This method is the most reliable.)
======================= =====================================================

------------
**Examples**
------------

The following are several examples of useful rotation specifications for surfaces:
::

   ROT = MESH   SURFACE   99   KINEMATIC   99   T2   0   T1   0   BASIS_RESEED

::

   ROT = MESH   SURFACE   16   T1   0   T2   0   PLANE   16   SEED   1.  0.  0.

::

   ROT = MOM   SURFACE   5   VELO_NORMAL 5   T1  0  T2   0   BASIS

The first example applies to the mesh equations at side set 99, the second to mesh equations at side set 16, and the third to the fluid momentum equations at side set 5. As described previously, the <string_x>, <string_y> and <string_z> parameters can be any boundary condition name or rotation string. Thus for the first example above, the x-component of the mesh equation is replaced by a KINEMATIC boundary condition on side set 99, the y-component of the mesh equation is replaced by the second tangential component (T2) of the mesh equation, and the z-component of the mesh equation is replaced by the first tangential component (T1) of the mesh equation. Since the rotation selections in the first example (T2 and T1) are rotated components instead of boundary conditions, a value of zero for the <int_y> and <int_z> parameters is appropriate. Finally, for the first example, BASIS_RESEED was chosen as the {seed_method}, and thus no subsequent parameters were required. The second example, however, uses SEED as the {seed_method} and thus is followed by the x, y, and z components of the tangent vector, respectively, as <float1> of 1., <float2> of 0., and <float3> of 0.

-------------------------
**Technical Discussion**
-------------------------

The necessary background discussing the nature and need for rotation procedures and rotation specifications is supplied in several of the references listed below. Briefly, however, in order to apply certain boundary conditions accurately it is necessary that the vector components of the solid mesh or fluid momentum equations be replaced by components that are tangent and normal to the surface in question. This procedure is referred to in this context as “rotation of equations.” It should be noted that explicitly specifying rotation conditions is really only necessary for three dimensional problems. Rotation also occurs in two-dimensional problems, but is sufficiently simpler that it can be automated and is therefore transparent to the user.

Not every boundary condition needs an accompanying rotation specification card and those that do are identified in the description of each boundary condition. Each rotated boundary condition will require at least one SURFACE rotation card be included for the boundary condition’s side set. Failure to do so is an error. The boundary conditions most often encountered that will require rotation cards are the *VELO_NORMAL* card applied to the fluid momentum equations and the *KINEMATIC, PLANE*, and *SPLINE* cards applied to the solid mesh equations.

In almost every case the boundary condition constraint will replace the normal rotated component so only the two tangential components of the rotated equation remain. All three examples shown above are just this situation. This has the effect of constraining the normal motion of the solid or fluid and imposing zero tangential forces due to the natural boundary conditions present in both fluid and solid momentum equations.

Specification of a seed vector method is needed so that a unique pair of tangent vectors may be determined at each point on the surface. The BASIS, BASIS_FIRST and BASIS_RESEED use the finite element grid in the surface as a means of defining the first tangent vector. They can employ averaging over elements that share a node. They should be employed on surfaces for which it is difficult to find a single consistent seed vector for every point on the surface. The SEED method finds the projection of the vector supplied in the surface at the point of interest. This projection vector is normalized to obtain the first tangent vector. It should be clear that only vectors that are never normal to any point on the surface will be suitable. In practice, this condition can sometimes be hard to meet for some surfaces. In these cases, the other seeding methods should be used.



--------------
**References**
--------------

GT-007.2: Tutorial on droplet on incline problem, July 30, 1999, T. A. Baer

GT-012.0: 3D Roll coating template and tutorial for GOMA, February 21, 2000, P.R. Schunk

GT-018.1: ROT card tutorial, January 22, 2001, T. A. Baer
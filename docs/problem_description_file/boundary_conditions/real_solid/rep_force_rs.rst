****************
**REP_FORCE_RS**
****************

::

	BC = REP_FORCE_RS SS <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR REALSOLID)**

This boundary condition card applies a force per unit area (traction) that varies as the
inverse of the fourth power of the distance from a planar surface (see Technical
Discussion below) on a *TALE* or Dynamic Lagrangian mesh region. This boundary
condition can be used to impose a normal contact condition (repulsion) or attraction
condition (negative force) between a planar surface and the surface of a *TALE* region. It
differs from REP_FORCE card only in the mesh-motion type to which it applies. The
force per unit area is applied uniformly over the boundary delineated by the side set ID.
The applied force is a vector in the normal direction to the Lagrangian interface.

Definitions of the input parameters are as follows, where <floatlist> has five
parameters:

================= =======================================================
**REP_FORCE_RS**  Name of the boundary condition (<bc_name>)
**SS**            Type of boundary condition (<bc_type>), where **SS**
                  denotes side set in the EXODUS II database.
<bc_id>           The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set
                  in EXODUS II) in the problem domain.
<float1>          Coefficient of repulsion, :math:`\lambda`.
<float2>          Coefficient *a* of plane equation.
<float3>          Coefficient *b* of plane equation.
<float4>          Coefficient *c* of plane equation.
<float5>          Coefficient *d* of plane equation.
================= =======================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = REP_FORCE_RS SS 10 1.e+03. 1.0 0.0 0.0 -3.0

This card results in a vector traction in the normal direction on surface side set 10
defined by :math:`\vec F` = –1.0e3 ⁄ :math:`h^4` where *F* is a force per unit area that varies with the distance
*h* from the plane specified by 1.0x – 3.0 = 0.0 .

-------------------------
**Technical Discussion**
-------------------------

The repulsive force is defined by :math:`\vec F` = F(:math:`\vec n`) where *F* is a force per unit area that varies
with the distance *h* from a plane defined by the equation ax + by + cz + d = 0. The
magnitude of the function :math:`\vec F` is defined as:

.. figure:: /figures/071_goma_physics.png
	:align: center
	:width: 90%

The normal vector is defined as the outward pointing normal to the surface. For internal
surfaces defined by side sets which include both sides of the interface, this condition
will result in exactly a zero traction, i.e., internal surface side sets must be attached to
one material only to get a net effect.

It is important to note that this boundary condition can only be applied to *TALE* mesh
motion types (cf. *Mesh Motion* card). As an example of how this boundary condition
card is used, consider the need to apply some load pressure uniformly on a surface that
is large enough such that this surface never penetrates a predefined planar boundary.
This condition hence can be use to impose an impenetrable contact condition.



--------------
**References**
--------------

GT-005.3: THE NEW TOTAL-ARBITRARY-LAGRANGIAN-EULERIAN (TALE)
CAPABILITY and its applicability to coating with/on deformable media, August 6,
1999, P. R. Schunk

.. 
	TODO - Image on line 62 needs to be taken out so the right equation can be written.
*************
**REP_FORCE**
*************

::

	BC = REP_FORCE SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MESH)**

This boundary condition card applies a force per unit area (traction) that varies as the
inverse of the fourth power of the distance from a planar surface to a Lagrangian or
dynamic Lagrangian mesh region. This boundary condition can be used to impose a
normal contact condition (repulsion) or attraction condition (negative force) between a
planar surface and the surface of a Lagrangian region. The force per unit area is applied
uniformly over the boundary delineated by the side set ID. The applied force is a vector
in the normal direction to the Lagrangian interface.

Definitions of the input parameters are as follows, with <float_list> having five
parameters:

=============== ==================================================================
**REP_FORCE**   Name of the boundary condition (<bc_name>)
**SS**          Type of boundary condition (<bc_type>), where **SS** denotes
                side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set in
                EXODUS II) in the problem domain.
<float1>        Coefficient of repulsion, λ.
<float2>        Coefficient *a* of plane equation.
<float3>        Coefficient *b* of plane equation.
<float4>        Coefficient *c* of plane equation.
<float5>        Coefficient *d* of plane equation.
=============== ==================================================================

Refer to the Technical Discussion for an explanation of the various coefficients.

------------
**Examples**
------------

The following sample card:
::

     BC = FORCE_REP SS   10  1.e+03.   1.0   0.0 0.0   -3.0

results in a vector traction of magnitude –1.0e3 ⁄ h\ :sup:`4` in the normal direction to 
surface
side set 10 and the distance h is measured from side set 10 to the plane defined by
1.0x – 3. = 0.

-------------------------
**Technical Discussion**
-------------------------

The *REP_FORCE* boundary condition produces a vector traction in the normal
direction to a surface side set, defined by:

.. math::

	\vec{F} = F (\vec{n}) = - \frac{\lambda}{h^4}

	

where *F* is a force per unit area that varies with the distance *h* from a plane defined 
by

.. math::

	ax + by + cz + d = 0



The normal vector is defined as the outward pointing normal to the surface. For internal
surfaces defined by side sets which include both sides of the interface, this condition
will result in exactly a zero traction, i.e., internal surface side sets must be attached 
to one element block only to get a net effect.

**Important note**: this boundary condition can only be applied to *LAGRANGIAN,
DYNAMIC_LAGRANGIAN* or *ARBITRARY* mesh motion types (cf. *Mesh Motion* card).
For real-solid mesh motion types, refer to *REP_FORCE_RS*. Furthermore, it is rare and
unlikely that this boundary condition be applied to *ARBITRARY* mesh motion regions.
An example application of this boundary condition card is to apply some load pressure
uniformly on a surface that is large enough such that this surface never penetrates a
predefined planar boundary. Hence, this condition can be use to impose an
impenetrable contact condition.


--------
**FAQs**
--------

On internal two-sided side sets, this boundary condition results in double the force in
the same direction.


.. 
	TODO - Equations need to replace the images in lines 63 and 70.
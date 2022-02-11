*********
**PLANE**
*********

::

	BC = PLANE SS <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MESH)**

This card is used to specify a surface (solid) boundary position of a planar surface. It is
applied as a rotated condition on the mesh equations (see *EQ* cards *mesh1, mesh2
mesh3*). The form of this equation is given by

.. math::

   f(x, y, z) = ax + by + cz + d = 0

Definitions of the input parameters are given below; note that <floatlist> has four
parameters corresponding to the four constants in the equation:

================== ========================================================================
**PLANE**          Name of the boundary condition name (<bc_name>).
**SS**             Type of boundary condition (<bc_type>), where **SS** denotes
                   side set in the EXODUS II database.
<bc_id>            The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (side set in
                   EXODUS II) in the problem domain.
<float1>           :math:`a` in function :math:`f(x, y, z)`
<float2>           :math:`b` in function :math:`f(x, y, z)`
<float3>           :math:`c` in function :math:`f(x, y, z)`
<float4>           :math:`d` in function :math:`f(x, y, z)`
================== ========================================================================

------------
**Examples**
------------

Following is a sample input card:
::

     BC = PLANE SS 3 0.0 1.0 0.0 -0.3

results in setting the side set elements along the side set 3 to a plane described by the
equation :math:`f(x, y, z, t) = y – 0.3 = 0` .

-------------------------
**Technical Discussion**
-------------------------

This, like most boundary conditions on geometry with arbitrary grid motion, is applied
to the weighted residuals of the mesh equation rotated into the normal-tangential basis
on the boundary. Specifically, this boundary condition displaces the normal component
after rotation of the vector residual equation, leaving the tangential component to
satisfy the natural mesh-stress free state. That is to say, this boundary condition allows
for mesh to slide freely in the tangential direction of the plane surface.

This boundary condition can be applied regardless of the *Mesh Motion* type, and is
convenient to use when one desires to move the plane with time normal to itself.



--------------
**References**
--------------

GT-001.4: GOMA and SEAMS tutorial for new users, February 18, 2002, P. R. Schunk
and D. A. Labreche

GT-013.2: Computations for slot coater edge section, October 10, 2002, T.A. Baer

.. 
	TODO - The image in line 19 needs to be repalced with the equation.
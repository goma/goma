****************
**MOVING_PLANE**
****************

::

	BC = MOVING_PLANE <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MESH)**

The *MOVING_PLANE* card is used to specify a surface (solid) boundary position
versus time for a planar surface (cf. *PLANE* boundary condition card). It is applied as a
rotated condition on the mesh equations (see *EQ* cards *mesh1, mesh2, mesh3*). The
form of the equation is given by

.. math::

   f(x, y, z, t) = ax + by + cz + d + g(t) = 0

and the function :math:`g(t)` is defined as

.. math::

   g(t) = \lambda_1t + \lambda_2t^2 + \lambda_3t^3

Definitions of the input parameters are given below; note that <floatlist> has seven
parameters corresponding to the seven constants in the above equations:

====================== ============================================================
**MOVING_ PLANE**      Name of the boundary condition name (<bc_name>).
**SS**                 Type of boundary condition (<bc_type>), where **SS** denotes
                       side set in the EXODUS II database
<bc_id>                The boundary flag identifier, an integer associated with
                       <bc_type> that identifies the boundary location (side set in
                       EXODUS II) in the problem domain.
<float1>               :math:`a` in function :math:`f(x, y, z, t)`
<float2>               :math:`b` in function :math:`f(x, y, z, t)`
<float3>               :math:`c` in function :math:`f(x, y, z, t)`
<float4>               :math:`d` in function :math:`f(x, y, z, t)`
<float5>               :math:`\lambda_1` coefficient in :math:`g(t)`
<float6>               :math:`\lambda_2` coefficient in :math:`g(t)`
<float7>               :math:`\lambda_3` coefficient in :math:`g(t)`
====================== ============================================================

------------
**Examples**
------------

The boundary condition card
::

     BC = MOVING_PLANE SS 3 0. 1. 0. -0.3 0.1 0.0 0.0

results in a plane originally positioned at y = 0.3 to move at a velocity of -*0.1*, viz.
the position of all nodes on the plane will follow:

.. math::

   f(x, y, z, t) = y - 0.3 + 0.1t = 0

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
convenient to use in place of *PLANE* when one desires to move the plane with time
normal to itself.




.. 
	TODO - The images in lines 20,26, and 61 need to be repalced with the equations.



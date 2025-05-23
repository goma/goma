*****************
**VELO_SLIP_ROT**
*****************

::

	BC = VELO_SLIP_ROT SS <bc_id> <float_list> [integer] [float5]

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is a variant of the *VELO_SLIP* boundary condition and serves
much the same function: to allow the fluid to slip relative to a solid substrate boundary.
The difference is that the assumed substrate is a rotating cylindrical surface with axis
parallel to the z-direction. Also as in the *VELO_SLIP* case, an optional variable slip
coefficient model is available that allows for slip to occur only in a region near to a
mesh node. This boundary condition is applicable generally only to two-dimensional
problems or very specialized three dimensional problems.

The <float_list> has four values and there are two optional values; definitions of the
input parameters are as follows:

================= ================================================================
**VELO_SLIP_ROT** Name of the boundary condition.
**SS**            Type of boundary condition (<bc_type>), where **SS**
                  denotes side set in the EXODUS II database.
<bc_id>           The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set
                  in EXODUS II) in the problem domain.
<float1>          :math:`\beta` the slip coefficient. The inverse of 
                  :math:`\beta` defines the
                  scaling between stress and slip. Hence, for small values
                  of :math:`\beta`, large shear stresses are needed 
                  for a given amount
                  of slip, and conversely, for large values of
                  :math:`\beta`, the amount
                  of stress needed for the same degree of slip decreases
                  (see below for a more rigorous description).
<float2>          :math:`\omega`, rotation rate of the cylindrical substrate 
                  surface in
                  radians/T. Positive values for this parameter correspond
                  to rotation in the clockwise direction.
<float3>          :math:`x_c`, the x-position of rotation axis.
<float4>          :math:`y_c`, the y-position of rotation axis.
[integer]         :math:`N_{cl}`, a single-node node set identification number.
                  When variable coefficient slip relation is used, distance is
                  measured relative to this node (see discussion below).
                  For problems involving dynamic contact lines, this
                  nodeset coincides with the location of the contact line.
[float5]          :math:`\alpha`, the distance scale in the variable slip model 
                  (see the discussion below). Both :math:`N_{cl}` and
                  :math:`\alpha` should be present to
                  activate the variable slip model.
================= ================================================================

------------
**Examples**
------------

The following is a sample card without the optional parameters:
::

     BC = VELO_SLIP_ROT SS 10 0.1 3.14 0.0 1.0

This condition specifies a moderate amount of slip (0.1) on a cylindrical surface
rotating at 3.14 rad/sec around the point (0.0,1.0).

-------------------------
**Technical Discussion**
-------------------------

The comments that appear in the Technical Discussion section of the *VELO_SLIP* card
apply equally well here. In particular, the discussion of the variable slip coefficient
model applies here as well. The only significant difference is that the velocity of the
substrate is not a fixed vector; instead, it is tangent to the cylindrical substrate with a
magnitude consistent with the radius of the cylinder and the rotation rate.



--------------
**References**
--------------

No References.
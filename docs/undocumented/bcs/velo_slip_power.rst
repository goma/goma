***************
VELO_SLIP_POWER
***************

::

	BC = VELO_SLIP_POWER SS <bc_id> <float_list>
	
-----------------------
Description / Usage
-----------------------

(WIC/VECTOR MOMENTUM)**

This boundary condition allows for slip between the fluid and a boundary using an
implementation of the Navier slip relation. This relation fixes the amount of slip as a
function of the applied shear stress. The scaling between stress and slip is a user
parameter.

The slip velocity is in the tangential direction :math:`t \cdot v_{slip}` and is
raised to a power (user param m)

There is an optional tangent vector to be input to the BC, the BC expects this vector
to be a unit tangent vector.

There are five required values in <float_list> and three optional values; definitions of
the input parameters are as follows:

**VELO_SLIP_POWER**
  Name of the boundary condition (<bc_name>).
  
**SS**
  Type of boundary condition (<bc_type>), where **SS** denotes
  side set in the EXODUS II database.
  
<bc_id>
  The boundary flag identifier, an integer associated with
  <bc_type> that identifies the boundary location (side set in
  EXODUS II) in the problem domain.

<float1>
  :math:`\beta`, the slip coefficient. The inverse of :math:`\beta` 
  defines the scaling between stress and slip. 

<float2>
  :math:`v_{s,x}`, the x-component of surface velocity vector. This would
  be the x-component of the fluid velocity if a no slip
  condition were applied.

<float3>
  :math:`v_{s,y}`, the y-component of surface velocity vector. This would
  be the y-component of the fluid velocity if a no slip
  condition were applied.

<float4>
  :math:`v_{s,z}`, the z-component of surface velocity vector. This would
  be the z-component of the fluid velocity if a no slip
  condition were applied.

<float5>
  :math:`m`, the power to raise the slip velocity to

[float6]
  :math:`t_x`, x-component of tangent vector

[float7]
  :math:`t_y`, y-component of tangent vector

[float8]
  :math:`t_z`, z-component of tangent vector

------------
Examples
------------

Following is a sample card without the optional parameters:
::

     BC = VELO_SLIP_POWER SS 10 0.1 0.0 0.0 0.0 2.0 

-------------------------
Technical Discussion
-------------------------

Boundary condition of the form:

.. math::

   \mathbf{n} \cdot \mathbf{\tau} = - \frac{1}{\beta}\left(\mathbf{t} \cdot (\mathbf{v} - \mathbf{v_s})\right)^m
   

----------
Theory
----------

No Theory.

--------
FAQs
--------

No FAQs.

--------------
References
--------------

*************************************
**Lower Velocity Function Constants**
*************************************

::

   Lower Velocity Function Constants = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card takes the specification of the Lower-wall velocity function for the confined
channel lubrication capability, or the lub_p equation. This function specifies the
velocity of the Lower channel wall as a function of time. Currently two models for
{model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model invokes a squeeze/separation velocity uniformly across the entire material|
|                          |region, viz. the two walls are brought together/apart at a constant rate. This option|
|                          |requires two floating point values                                                   |
|                          |                                                                                     |
|                          | * <float1> is the velocity component in the x-direction. L/t                        |
|                          | * <float2> is the velocity component in the y-direction. L/t                        |
|                          | * <float3> is the velocity component in the z-direction. L/t (NOTE: this is usually |
|                          |   taken as zero as it is set in the Lower Wall Height Function model)               |
+--------------------------+-------------------------------------------------------------------------------------+
|**SLIDER_POLY_TIME**      |This model implements a spatially-uniform velocity in the x-direction that is        |
|                          |specified as a polynomial in time. The value of time may be scaled by a given scaling|
|                          |factor and the polynomial may have an unlimited number of terms.                     |
|                          |                                                                                     |
|                          | * <float1> is the time scaling factor                                               |
|                          | * <float2-N> are the coefficients in front of the t^(i-2) term                      |
+--------------------------+-------------------------------------------------------------------------------------+

.. figure:: /figures/490_goma_physics.png
	:align: center
	:width: 90%

+--------------------------+-------------------------------------------------------------------------------------+
|**ROLL**                  |This model invokes a wall velocity which corresponds to a rolling-motion. This model |
|                          |takes nine constants ???? :                                                          |
|                          |                                                                                     |
|                          | * <float1> Roll radius, L.                                                          |
|                          | * <float2> x-coordinate of axis origin, L.                                          |
|                          | * <float3> y-coordinate of axis orgin, L.                                           |
|                          | * <float4> z-coordinate of axis origin, L.                                          |
|                          | * <float5> Direction angle 1 of rotation axis                                       |
|                          | * <float6> Direction angle 2of rotation axis                                        |
|                          | * <float7> Direction angle 3of rotation axis                                        |
|                          | * <float8> Squeeze rate.                                                            |
|                          | * <float9> rotation rate                                                            |
+--------------------------+-------------------------------------------------------------------------------------+
|**TANGENTIAL_ROTATE**     |This model allows a unique specification of tangential motion in a lubrication shell |
|                          |element. Previous implementations allowed specification only in terms of coordinate  |
|                          |direction, but this option can be used to rotate a cylinder. Five floats are required|
|                          |                                                                                     |
|                          | * <float1> x-comnponent of a vector tangential to the shell. This vector must never |
|                          |   be normal to the shell. It is then projected onto the shell.                      |
|                          | * <float2> y-comnponent of a vector tangential to the shell.                        |
|                          | * <float3> z-comnponent of a vector tangential to the shell.                        |
|                          | * <float4> U1, or scalar speed of wall velocity in a direction determined by the    |
|                          |   cross product ot the tangent vector and the normal vector to the shell. (L/t)     |
|                          | * <float5> U2 scalar speed component in direction normal to U1. (L/t)               |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Lower Velocity Function Constants = CONSTANT {v_x= -0.001} {vy=0.00} {vz=0}

This card results in an Lower wall speed of -0.001 in the x-direction which is tangential
to the substrate, thus generating a Couette component to the flow field.

-------------------------
**Technical Discussion**
-------------------------

For non-curved shell meshes, most of the time they are oriented with the x-, y-, or zplane.
This card is aimed at applying a tangential motion to that plane, and so one of
the three components is usually zero.




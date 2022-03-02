*************************************
**Upper Velocity Function Constants**
*************************************

::

   Upper Velocity Function Constants = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card takes the specification of the upper-wall velocity function for the confined
channel lubrication capability, or the lub_p equation. This function specifies the
velocity of the upper channel wall as a function of time. Currently two models for
{model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model invokes a squeeze/separation velocity uniformly across the entire material|
|                          |region, viz. the two walls are brought together/apart at a constant rate. This option|
|                          |requires two floating point values                                                   |
|                          |                                                                                     |
|                          | * <float1> is the velocity component in the x-direction. L/t.                       |
|                          | * <float2> is the velocity component in the y-direction. L/t                        |
|                          | * <float3> is the velocity component in the z-direction. L/t (NOTE: this is usually |
|                          |   taken as zero as it is set in the Upper Wall Height Function model)               |
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
|                          | * <float8> Squeeze rate                                                             |
|                          | * <float9> rotation rate                                                            |
+--------------------------+-------------------------------------------------------------------------------------+
|**TANGENTIAL_ROTATE**     |his model allows a velocity that is always tangential to a shell surface, not        |
|                          |necessarily aligned along the coordinate directions. It requires three               |
|                          |specifications. First, a vector (v) that is always non-colinear to the normal vector |
|                          |of the shell must be specified. This is used to make unique tangent vectors. The last|
|                          |two specifications are the two tangential components to the velocity. The first      |
|                          |velocity is applied in the direction of t1 = v×n. The second velocity is then applied|
|                          |in the t = t ×n  direction.                                                          |
|                          |                                                                                     |
|                          | * <float1> vx                                                                       |
|                          | * <float2> vy                                                                       |
|                          | * <float3> vz                                                                       |
|                          | * <float4> velocity in the t1 direction                                             |
|                          | * <float5> velocity in the t2 direction                                             |
+--------------------------+-------------------------------------------------------------------------------------+
|**CIRCLE_MELT**           |Model which allows a converging or diverging height that is like a circle. Also works|
|                          |for melting.                                                                         |
|                          |                                                                                     |
|                          | * <float1> - x-location of the circle center (circle is in x-y plane)               |
|                          | * <float2> - radius of circle                                                       |
|                          | * <float3>- minimum height of circle                                                |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Upper Velocity Function Constants = CONSTANT {v_x= -0.001} {vy=0.00} {vz=0}

This card results in an upper wall speed of -0.001 in the x-direction which is tangential
to the substrate, thus generating a Couette component to the flow field.

-------------------------
**Technical Discussion**
-------------------------

For non-curved shell meshes, most of the time they are oriented with the x-, y-, or zplane.
This card is aimed at applying a tangential motion to that plane, and so one of
the three components is usually zero.




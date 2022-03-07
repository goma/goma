***********************************
**Upper Height Function Constants**
***********************************

::

   Upper Height Function Constants = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card takes the specification of the upper-height function for the confined channel
lubrication capability, or the lub_p equation. This function specifies the height of the
channel versus distance and time. Currently three models for {model_name} are
permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT_SPEED**        |This model invokes a squeeze/separation velocity uniformly across the entire material|
|                          |region, viz. the two walls are brought together/apart at a constant rate. This option|
|                          |requires two floating point values                                                   |
|                          |                                                                                     |
|                          | * <float1> the separation velocity (rate) in units of length/time                   |
|                          | * <float2> the initial wall separation in units of length                           |
|                          | * <float3> An OPTIONAL parameter which scales the addition of an external field     |
|                          |   called “HEIGHT” which is read in using the External Field or External Pixel Field |
|                          |   capabilities. If this field is present, the value of it is added to the height    |
|                          |   calculated with this model.                                                       |
+--------------------------+-------------------------------------------------------------------------------------+
|**ROLL_ON**               |This model invokes a squeeze/separation velocity in a hinging-motion along one       |
|                          |boundary. The model is best explained with the figure in the technical discussion    |
|                          |section. The equation for the gap h as a function of time and the input parameters   |
|                          |(floats) is as follows:                                                              |
|                          |                                                                                     |
|                          | * <float1> is x0 in units of length                                                 |
|                          | * <float2> is hlow in units of length                                               |
|                          | * <float3> is h Δ, in units of length                                               |
|                          | * <float4> is the verticle separation velocity (if negative then squeeze velocity)  |
|                          |   in units of length/time                                                           |
|                          | * <float5> is the length of the plate, L.                                           |
+--------------------------+-------------------------------------------------------------------------------------+

.. figure:: /figures/486_goma_physics.png
	:align: center
	:width: 90%

+--------------------------+-------------------------------------------------------------------------------------+
|**ROLL**                  |This model is used for a roll coating geometry. This option requires 8 floats:       |
|                          |                                                                                     |
|                          | * <float1> x-coordinate of origin, L.                                               |
|                          | * <float2> y-coordinate of orgin, L.                                                |
|                          | * <float3> z-coordinate of origin, L.                                               |
|                          | * <float4> Direction angle 1 of rotation axis                                       |
|                          | * <float5> Direction angle 2of rotation axis                                        |
|                          | * <float6> Direction angle 3of rotation axis                                        |
|                          | * <float7> rotation speed L/t.                                                      |
+--------------------------+-------------------------------------------------------------------------------------+
|**FLAT_GRAD_FLAT**        |This model used two arctan functions to mimic a flat region, then a region of        |
|                          |constant slope, then another flat region. The transitions between the two regions are|
|                          |curved by the arctan function. This currently on works for changes in the x          |
|                          |direction. This option requires five floating point values                           |
|                          |                                                                                     |
|                          | * <float1> x location of the first transition (flat to grad)                        |
|                          | * <float2> height of the first flat region                                          |
|                          | * <float3> x location of the second transition (grad to flat)                       |
|                          | * <float4> height of the second flat region                                         |
|                          | * <float5> parameter controlling the curvature of the transitions                   |
+--------------------------+-------------------------------------------------------------------------------------+
|**POLY_TIME**             |This time applies a time-dependent lubrication height in the form of a polynomial. It|
|                          |can take as many arguments as GOMA can handle, and the resulting height function is  |
|                          |                                                                                     |
|                          | * <floati> value of Ci                                                              |
+--------------------------+-------------------------------------------------------------------------------------+

.. figure:: /figures/487_goma_physics.png
	:align: center
	:width: 90%

+--------------------------+-------------------------------------------------------------------------------------+
|**JOURNAL**               |This model simulates a journal bearing. It is intended to be run on a cylindrical    |
|                          |shell mesh aligned along the z axis and centered at (0,0). It could be extended to be|
|                          |more flexible, but this is all it is currently capable of. The height is defined by  |
|                          |                                                                                     |
|                          |h(θ ) = C(1+ε cos(0))                                                                |
|                          |                                                                                     |
|                          |Where C is the mean lubrication height and is the eccentricity of the two cylinders, |
|                          |with the smallest gap in the –y direction.                                           |
|                          |                                                                                     |
|                          | * <float1> C                                                                        |
|                          | * <float2> ε                                                                        |
+--------------------------+-------------------------------------------------------------------------------------+
|**EXTERNAL_FIELD**        |Not recognized. Oddly, this model is invoked with the extra optional float on the    |
|                          |CONSTANT_SPEED option.                                                               |
+--------------------------+-------------------------------------------------------------------------------------+

::

   External Field = HEIGHT Q1 name.exoII (see this card)

------------
**Examples**
------------

Following is a sample card:

::

   Upper Height Function Constants = CONSTANT_SPEED {v_sq = -0.001} {h_i=0.001}

This results in an upper wall speed of 0.001 in a direction which reduces the gap, which
is initial 0.001.

-------------------------
**Technical Discussion**
-------------------------

The material function model ROLL_ON prescribes the squeezing/separation motion of
two non-parallel flate plates about a hinge point, as shown in the figure below.




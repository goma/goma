***********************************
**Lower Height Function Constants**
***********************************

::

   Lower Height Function Constants = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card takes the specification of the lower-height function for the confined channel
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

.. figure:: /figures/488_goma_physics.png
	:align: center
	:width: 90%

+--------------------------+-------------------------------------------------------------------------------------+
|**ROLL**                  |This model is used for a roll coating geometry. This option requires 8 floats:       |
|                          |                                                                                     |
|                          | * <float1> x-coordinate of origin, L.                                               |
|                          | * <float2> y-coordinate of orgin, L.                                                |
|                          | * <float3> z-coordinate of origin, L.                                               |
|                          | * <float4> Direction angle 1 of rotation axis                                       |
|                          | * <float5> Direction angle 2 of rotation axis                                       |
|                          | * <float6> Direction angle 3 of rotation axis                                       |
|                          | * <float7> rotation speed L/t.                                                      |
+--------------------------+-------------------------------------------------------------------------------------+
|**TABLE** <integer1>      |Please see discussion at the beginning of the material properties Chapter 5 for input|
|<character_string1>       |description and options. Most likely *character_string1* will be **LOWER_DISTANCE**  |
|{LINEAR | BILINEAR}       |This option is good for inputing table geometry versus distance. Specifically, an    |
|[integer2] [FILE = filenm]|arbitrary lower height function model is input as a function of the x-direction      |
|                          |coordinate of the Lower Velocity Function model. This option in turn requires the use|
|                          |of SLIDER_POLY_TIME lower velocity function model. See example below.                |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Lower Height Function Constants = CONSTANT_SPEED {v_sq = -0.001} {h_i=0.001}

This results in an lower wall speed of 0.001 in a direction which reduces the gap, which
is initial 0.001.

In another example:

Lower Height Function Constants = TABLE 2 LOWER_DISTANCE 0
LINEAR FILE=shell.dat

where shell.dat is a table with 2 columns, the first the position, the second the height.

-------------------------
**Technical Discussion**
-------------------------

The material function model ROLL_ON prescribes the squeezing/separation motion of
two non-parallel flate plates about a hinge point, as shown in the figure below.





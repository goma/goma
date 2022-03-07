**************************
**Film Evaporation Model**
**************************

::

   Film Evaporation Model = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card takes the specification of the evaporation rate for the film-flow equation
capability, specifically the shell_filmp equation. This function specifies the rate
of evaporation in the unit of length per time. Currently two models for
{model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model specifies a constant evaporation rate. This option requires one floating  |
|                          |point values                                                                         |
|                          |                                                                                     |
|                          | * <float1> is the evaporation rate in the unit of length/time                       |
+--------------------------+-------------------------------------------------------------------------------------+
|**CONC_POWER**            |This model specifies evaporation rate function of particles volume fraction and the  |
|                          |input parameters (floats). This model is proposed by Schwartz et al (2001). The      |
|                          |functional form is:                                                                  |
|                          |                                                                                     |
|                          | * <float1> is the pure liquid evaporation rate in units of length per time          |
|                          | * <float2> is exponent v and it should satisfy 0 <ν <1                              |
|                          | * <float3> is the maximum packing volume fraction 0max                              |
+--------------------------+-------------------------------------------------------------------------------------+

.. figure:: /figures/494_goma_physics.png
	:align: center
	:width: 90%

------------
**Examples**
------------

Following is a sample card:

::

   Film Evaporation Model = CONC_POWER 1.0e-3 0.5 0.64

This results in film evaporation with the pure liquid evaporation rate of 1.0e-3, exponent
of 0.5, and maximum packing volume fraction of 0.64.




--------------
**References**
--------------

Leonard W. Schwartz, R. Valery Roy, Richard R. Eley, and Stanislaw Petrash,
“Dewetting Patterns in a Drying Liquid Film”, Journal of Colloid and Interface
Science 234, 363–374 (2001)

****************
Thermal Exponent
****************

::

   Thermal Exponent = CONSTANT <float> [T]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify a thermal exponential factor for **CARREAU_WLF, BINGHAM, THERMAL, EPOXY, FILLED_EPOXY, POWERLAW_SUSPENSION** and **CARREAU_SUSPENSION**
viscosity models, as selected in the *Liquid Constitutive Equation* card. The value represented by the thermal exponent varies between these liquid constitutive 
models; the appropriate values for each model is indicated below.

Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the thermal exponent.                                                                 |
|                 |                                                                                                            |
|                 | * <float> - the value of the thermal exponent for the viscosity model specified in the *Liquid Constitutive|
|                 |   Equation* card.                                                                                          |
|                 |                                                                                                            |
|                 | * for the **BINGHAM, THERMAL, EPOXY,** or **FILLED_EPOXY** model,                                          |
|                 |                                                                                                            |
|                 |   * <float> - the Eμ/R parameter. This has the dimensions of temperature in whatever units are consistent  |
|                 |     with the problem and describes the thinning of viscosity with temperature.                             |
|                 |                                                                                                            |
|                 | * for the **CARREAU_WLF** model model,                                                                     |
|                 |                                                                                                            |
|                 |   * <float> - the c1 constant of the equation for the temperature-dependent shift factor.                  |
|                 |                                                                                                            |
|                 | * for the **POWERLAW_SUSPENSION** or **CARREAU_SUSPENSION** model,                                         |
|                 |                                                                                                            |
|                 |   * <float> - the exponent for the Krieger viscosity model, *m*.                                           |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**LEVEL_SET**    |Name of the model for level-set dependent thermal exponent factor. Allows for this exponent level to be a   |
|                 |function of the level-set field. Specifically used for changing the thermal exponent from one constant value|
|                 |on the negative side of the interface to another constant value on the positive side. The model requires    |
|                 |three floats:                                                                                               |
|                 |                                                                                                            |
|                 | * <float1> - the value of thermal exponent in the negative regions of the level set function.              |
|                 | * <float2> - the value of thermal exponent in the positive regioons of the level-set function.             |
|                 | * <float3> Length scale over which the transition occurs. If this parameter is set to zero, it will default|
|                 |   to one-half the Level-Set Length Scale value specified.                                                  |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the thermal exponent to 0.5.

::

   Thermal Exponent = CONSTANT 0.5

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.





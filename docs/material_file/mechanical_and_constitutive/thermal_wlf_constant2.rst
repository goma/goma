*********************
Thermal WLF Constant2
*********************

::

   Thermal WLF Constant2 = CONSTANT <float> [T]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the thermal constant 2 of the **CARREAU_WLF** viscosity
model in the *Liquid Constitutive Equation* card. Definitions of the input parameters are
as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for Thermal Constant2.                                                                    |
|                 |                                                                                                            |
|                 | * <float> - the value of :math:`c_2`, in the equation representing the temperature-dependent shift factor  |
|                 |   for the **CARREAU_WL** constitutive model.                                                               |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**LEVEL_SET**    |Name of the model for level-set dependent WLF thermal constant 2. Allows for this thermal constant 2 level  |
|                 |to be a function of the level-set field. Specifically used for changing the thermal constant 2 from one     |
|                 |constant value on the negative side of the interface to another constant value on the positive side. The    |
|                 |model requires three floats:                                                                                |
|                 |                                                                                                            |
|                 | * <float1> - the value of thermal constant 2 in the negative regions of the level set function.            |
|                 |                                                                                                            |
|                 | * <float2> - the value of thermal constant 2 in the positive regioons of the level-set function.           |
|                 |                                                                                                            |
|                 | * <float3> Length scale over which the transition occurs. If this parameter is set to zero, it will default|
|                 |   to one-half the Level-Set Length Scale value specified.                                                  |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the *Thermal WLF Constant2* to 0.1.

::

   Thermal WLF Constant2 = CONSTANT 0.1

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.




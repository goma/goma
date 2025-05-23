******************
Low Rate Viscosity
******************

::

   Low Rate Viscosity = CONSTANT <float> [M/Lt]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the low-rate viscosity parameter for the **POWER_LAW, CARREAU, CARREAU_WLF, BINGHAM, SUSPENSION,**
**THERMAL, CURE, EPOXY, FILLED_EPOXY, POWERLAW_SUSPENSION** and **CARREAU_SUSPENSION** model options of the *Liquid Constitutive Equation*
card. This is also the reference viscosity value in the **HERSCHEL_BULKLEY** constitutive equation.

Definitions of the input parameters are as follows:

+-----------------+----------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the low-rate viscosity.                                                       |
|                 |                                                                                                    |
|                 | * <float> - the value of the low-rate viscosity. This value is also called the zero strain-rate    |
|                 |   limit of the viscosity and in models is normally called μ0.                                      |
+-----------------+----------------------------------------------------------------------------------------------------+
|**LEVEL_SET**    |Name of the model for level-set dependent low-rate viscosity. Allows for this viscosity level to be |
|                 |a function of the level-set field. Specifically used for changing the low-rate viscosity from one   |
|                 |constant value on the negative side of the interface to another constant value on the positive side.|
|                 |The model requires three floats:                                                                    |
|                 |                                                                                                    |
|                 | * <float1> - the value of viscosity in the negative regions of the level set function.             |
|                 | * <float2> - the value of viscosity in the positive regioons of the level-set function.            |
|                 | * <float3> Length scale over which the transition occurs. If this parameter is set to zero, it will|
|                 |   default to one-half the Level-Set Length Scale value specified.                                  |
+-----------------+----------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the low rate viscosity to 10:

::

   Low Rate Viscosity = CONSTANT 10.

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.



--------------
**References**
--------------

No References.
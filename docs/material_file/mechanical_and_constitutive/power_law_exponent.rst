******************
Power Law Exponent
******************

::

   Power Law Exponent = CONSTANT <float> []

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the power-law exponent parameter of the **POWER_LAW, CARREAU, BINGHAM, CARREAU_WLF, CURE,**
**SUSPENSION, FILLED_EPOXY, POWERLAW_SUSPENSION, CARREAU_SUSPENSION,** and **HERSCHEL_BULKLEY** fluid options of the *Liquid Constitutive Equation* card.

Definitions of the input parameters are as follows:

+-----------------+---------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the power-law exponent.                                                            |
|                 |                                                                                                         |
|                 | * <float> - the value of the power-law exponent. This variable is normally *n* in the constitutive laws.|
+-----------------+---------------------------------------------------------------------------------------------------------+
|**LEVEL_SET**    |Name of the model for level-set dependent power law exponent. Specifically used for changing the exponent|
|                 |from one constant value on the negative side of the interface to another constant value on the positive  |
|                 |side. The model requires three floats:                                                                   |
|                 |                                                                                                         |
|                 | * <float1> - the value of power-law exponent in the negative regions of the level set function.         |
|                 | * <float2> - the value of power-law exponent in the positive regioons of the level-set function.        |
|                 | * <float3> Length scale over which the transition occurs. If this parameter is set to zero, it will     |
|                 |   default to one-half the Level-Set Length Scale value specified.                                       |
+-----------------+---------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the power law exponent to 0.2:

::

   Power Law Exponent = CONSTANT 0.2

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.



--------------
**References**
--------------

No References.
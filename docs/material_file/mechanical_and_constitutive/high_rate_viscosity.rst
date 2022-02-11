*******************
High Rate Viscosity
*******************

::

   High Rate Viscosity = CONSTANT <float> [M/Lt]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the high-rate viscosity parameter of the **CARREAU, BINGHAM, CARREAU_WLF** and **CARREAU_SUSPENSION** fluid
options of the *Liquid Constitutive Equation* card. Definitions of the input parameters are as follows:

+-----------------+---------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the high-rate viscosity.                                                           |
|                 |                                                                                                         |
|                 | * <float> - the value of the high-rate viscosity. This value is normally called μinf in models.         |
+-----------------+---------------------------------------------------------------------------------------------------------+
|**LEVEL_SET**    |Name of the model for level-set dependent high-rate viscosity. Allows for this viscosity level to be a   |
|                 |function of the level-set field. Specifically used for changing the high-rate viscosity from one constant|
|                 |value on the negative side of the interface to another constant value on the positive side. The model    |
|                 |requires three floats:                                                                                   |
|                 |                                                                                                         |
|                 | * <float1> - the value of viscosity in the negative regions of the level set function.                  |
|                 | * <float2> - the value of viscosity in the positive regioons of the level-set function.                 |
|                 | * <float3> Length scale over which the transition occurs. If this parameter is set to zero, it will     |
|                 |   default to one-half the Level-Set Length Scale value specified.                                       |
+-----------------+---------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the high rate viscosity to 10.:

::

   High Rate Viscosity = CONSTANT 10.

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.



--------------
**References**
--------------

No References.
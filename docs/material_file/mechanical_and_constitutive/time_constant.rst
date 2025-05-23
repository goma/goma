*************
Time Constant
*************

::

   Time Constant = CONSTANT <float> [t]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the time constant parameter of the **CARREAU, BINGHAM, CARREAU_WLF** and **CARREAU_SUSPENSION** fluid
options of the *Liquid Constitutive Equation* card. Definitions of the input parameters are as follows:

+-----------------+---------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the time constant.                                                                 |
|                 |                                                                                                         |
|                 | * <float> - the value of the time constant, λ.                                                          |
+-----------------+---------------------------------------------------------------------------------------------------------+
|**LEVEL_SET**    |Name of the model for level-set dependent time constant. Allows for this time constant level to be a     |
|                 |function of the level-set field. Specifically used for changing the time constant from one constant value|
|                 |on the negative side of the interface to another constant value on the positive side. The model requires |
|                 |three floats:                                                                                            |
|                 |                                                                                                         |
|                 | * <float1> - the value of time constant in the negative regions of the level set function.              |
|                 | * <float2> - the value of time constant in the positive regioons of the level-set function.             |
|                 | * <float3> Length scale over which the transition occurs. If this parameter is set to zero, it will     |
|                 |   default to one-half the Level-Set Length Scale value specified.                                       |
+-----------------+---------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the time constant to 0.2.

::

   Time Constant = CONSTANT 0.2

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.




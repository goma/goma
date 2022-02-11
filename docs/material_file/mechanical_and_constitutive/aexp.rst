****
Aexp
****

::

   Aexp = CONSTANT <float> []

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the Aexp parameter of the **CARREAU,BINGHAM, CARREAU_WLF** and **CARREAU_SUSPENSION** model options of
the *Liquid Constitutive Equation* card. Definitions of the input parameters are as follows:

+-----------------+---------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for Aexp.                                                                              |
|                 |                                                                                                         |
|                 | * <float> - the value of the a exponent in the liquid constitutive models; also, a dimensionless        |
|                 |   parameter that describes the transition between the low-rate and the power-law region for the Carreau |
|                 |   model (see Bird, et. al., 1987).                                                                      |
+-----------------+---------------------------------------------------------------------------------------------------------+
|**LEVEL_SET**    |Name of the model for level-set dependent Aexp parameter. Allows for this parameter level to be a        |
|                 |function of the level-set field. Specifically used for changing the Aexp parameter from one constant     |
|                 |value on the negative side of the interface to another constant value on the positive side. The model    |
|                 |requires three floats:                                                                                   |
|                 |                                                                                                         |
|                 | * <float1> - the value of Aexp parameter in the negative regions of the level set function.             |
|                 | * <float2> - the value of Aexp parramete in the positive regioons of the level-set function.            |
|                 | * <float3> Length scale over which the transition occurs. If this parameter is set to zero, it will     |
|                 |   default to one-half the Level-Set Length Scale value specified.                                       |
+-----------------+---------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets Aexp to 3.0:

::

   Aexp = CONSTANT 3.0

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.



--------------
**References**
--------------

Bird, R. B., Armstrong, R. C., and Hassager, O. 1987. Dynamics of Polymeric Liquids,
2nd ed., Wiley, New York, Vol. 1.

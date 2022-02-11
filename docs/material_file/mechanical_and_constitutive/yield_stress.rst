************
Yield Stress
************

::

   Yield Stress = CONSTANT <float> [M/Lt2]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the yield stress parameter, :math:`τ_y`, of the
**BINGHAM** and **HERSCHEL_BULKLEY** model options of the *Liquid Constitutive
Equation* card. Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the yield stress.                                                                     |
|                 |                                                                                                            |
|                 | * <float> - the value of the yield stress, τy, which has the dimensions of stress, in whatever units are   |
|                 |   consistent with the problem and marks the transition from solid-like to fluid-like behavior for the      |
|                 |   Bingham-Carreau-Yasuda model.                                                                            |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the yield stress to 100:

::

   Yield Stress = CONSTANT 100.

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.




***************************
**Thermodynamic Potential**
***************************

::

   Thermodynamic Potential = {model_name} {float_list}

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the anodic or cathodic thermodynamic potential in a
thermal battery cell.

+----------------------+-------------------------------------------------------------------------------------+
|{model_name}          |Name of thermodynamic potential model. Currently, two thermodynamic potential models |
|                      |are available, namely **LiSi** and **FeS2**. Each of these and their accompanying    |
|                      |input parameters (the <float_list>) is given below:                                  |
+----------------------+-------------------------------------------------------------------------------------+
|**LiSi**              |This model requires seven floating-point parameters:                                 |
|                      |                                                                                     |
|                      | * <float1> - Limit of electrode utilization for the first anode reaction.           |
|                      | * <float2> - Limit of electrode utilization for the second anode reaction.          |
|                      | * <float3> - Anode thickness.                                                       |
|                      | * <float4> - Anode porosity.                                                        |
|                      | * <float5> - Molar volume of active anode material.                                 |
|                      | * <float6> - Current density output by the thermal battery cell.                    |
|                      | * <float7> - Number of electrons involved in anode reactions.                       |
+----------------------+-------------------------------------------------------------------------------------+
|**FeS2**              |This model requires eight floating-point parameters:                                 |
|                      |                                                                                     |
|                      | * <float1> - Limit of electrode utilization for the first cathode reaction.         |
|                      | * <float2> - Limit of electrode utilization for the second cathode reaction.        |
|                      | * <float3> - Limit of electrode utilization for the third cathode reaction.         |
|                      | * <float4> - Cathode thickness.                                                     |
|                      | * <float5> - Cathode porosity.                                                      |
|                      | * <float6> - Molar volume of active cathode material.                               |
|                      | * <float7> - Current density output by the thermal battery cell.                    |
|                      | * <float8> - Number of electrons involved in cathode reactions.                     |
+----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following are two sample input cards:

::

   Thermodynamic Potential = LiSi 0.283 0.474 0.088 0.275 54.61 0.0246 3.25

::

   Thermodynamic Potential = FeS2 0.375 0.434 0.5 0.046 0.244 23.93 0.0246 4.0

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




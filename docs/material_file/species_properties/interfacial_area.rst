********************
**Interfacial Area**
********************

::

   Interfacial Area = {model_name} {float_list}

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the product of interfacial area per unit volume by exchange
current density (i.e., ai0) in the Butler-Volmer kinetic model of current density.

+----------------------+-------------------------------------------------------------------------------------+
|{model_name}          |Name of the model for interfacial area, of which there are currently two available,  |
|                      |namely **CONSTANT** and **THERMAL_BATTERY**.                                         |
+----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**          |constant value of Interfacial Area                                                   |
|                      |                                                                                     |
|                      | * <float1> - the value of the product of interfacial area per unit volume and       |
|                      |   exchange current density.                                                         |
+----------------------+-------------------------------------------------------------------------------------+
|**THERMAL_BATTERY**   |this option requires the following nine parameters:                                  |
|                      |                                                                                     |
|                      | * <float1> - Initial value of the product of interfacial area per unit volume by    |
|                      |   exchange current density.                                                         |
|                      | * <float2> - Limit of electrode utilization beyond which ai0 = 0.                   |
|                      | * <float3> - Activation energy for the Arrhenius dependency of ai0 on temperature.  |
|                      | * <float4> - Initial electrode/electrolyte temperature.                             |
|                      | * <float5> - Cathode thickness.                                                     |
|                      | * <float6> - Cathode porosity.                                                      |
|                      | * <float7> - Molar volume of active cathode material.                               |
|                      | * <float8> - Current density output by the thermal battery cell.                    |
|                      | * <float9> - Number of electrons involved in cathode reactions.                     |
+----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following are two sample input cards:

::

   Interfacial Area = CONSTANT 1.0

::

   Interfacial Area = THERMAL_BATTERY 20.0 0.375 20000.0 846.0 0.046 0.244 23.93 0.0246 4.0

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
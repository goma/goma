************************
**Solution Temperature**
************************

::

   Solution Temperature = {model_name} <float_list>

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the temperature of an electrolyte solution (i.e., when
diffusion and migration transport of charged species is involved).

+----------------------+-------------------------------------------------------------------------------------+
|{model_name}          |Name of the electrolyte-solution model, for which there are currently two options:   |
|                      |**CONSTANT** and **THERMAL_BATTERY**; the former model has a single parameter in the |
|                      |<float_list> while the latter has six.                                               |
+----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**          |A constant model of the solution temperature.                                        |
|                      |                                                                                     |
|                      | * <float1> - the value of electrolyte-solution temperature.                         |
+----------------------+-------------------------------------------------------------------------------------+
|**THERMAL_BATTERY**   |A specialized model of electrolyte solutions for Thermal Batteries                   |
|                      |(Chen, et. al., 2000).                                                               |
|                      |                                                                                     |
|                      | * <float1> - Initial electrolyte solution temperature ( K )                         |
|                      | * <float2> - Ambient temperature ( K )                                              |
|                      | * <float3> - Cross-sectional area from which heat is lost to ambient ( m2           |
|                      | * <float4> - Heat transfer coefficient ( W/m2/K )                                   |
|                      | * <float5> - Mass of battery cell ( kg )                                            |
|                      | * <float6> - Heat capacity of electrolyte solution ( J/kg/K)                        |
+----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following are two sample input cards:

::

   Solution Temperature = CONSTANT 313.0

::

   Solution Temperature = THERMAL_BATTERY 846. 298. 0.0316 7.7 0.6 1030.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

SAND2000-0207: Final Report on LDRD Project: A Phenomenological Model for
Multicomponent Transport with Simultaneous Electrochemical Reactions in
Concentrated Solutions, Chen, K. S., Evans, G. H., Larson, R. S., Noble, D. R., and
Houf, W. G., January 2000.
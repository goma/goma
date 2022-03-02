******************
**Vapor Pressure**
******************

::

   Vapor Pressure = {model_name} <species> {float_list} [varies]

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for the vapor pressure for each species;
it has two main classes of use. The first class regards multiphase flow in porous media,
which is activated when the media type is set to **POROUS_UNSATURATED** or
**TWO_PHASE** (cf. the *Media Type* card). The second class of use of this data card is
for specification of vapor pressure at the external boundary of a liquid domain, for
which the bounding gas phase is modeled with a lumped parameter approach, or at an
internal interface between a liquid and a gas. No curvature effects are included here.
Eventually the models in this class will be supported in the porous-media cases.
Definitions of the input parameters are as follows:

+-----------------------+-------------------------------------------------------------------------------------+
|{model_name}           |Name of the model for the vapor pressure, based on the class of use.                 |
|                       |                                                                                     |
|                       |For the first class of multiphase flows in porous media, {model_name} can be one of  |
|                       |the following:                                                                       |
|                       |                                                                                     |
|                       | * **KELVIN** - for a volatile liquid                                                |
|                       | * **IDEAL_GAS** - for a non-condensable gas                                         |
|                       | * **FLAT** - for a volatile liquid                                                  |
|                       |                                                                                     |
|                       |For the second class regarding specification of vapor pressure at the external       |
|                       |boundary of a liquid domain or the interface between a gas and a liquid, {model_name}|
|                       |can be one of the following:                                                         |
|                       |                                                                                     |
|                       | * **CONSTANT** - for a constant vapor pressure model                                |
|                       | * **ANTOINE** - for temperature-dependent, nonideal gases                           |
|                       | * **RIEDEL** -for temperature-dependent, nonideal gases                             |
+-----------------------+-------------------------------------------------------------------------------------+
|<species>              |An integer designating the species equation. Typically this value is zero if the     |
|                       |problem is one of a single solvent in a partially saturated medium.                  |
+-----------------------+-------------------------------------------------------------------------------------+
|{float_list}           |One or more floating point numbers (<float1> through <floatn>) whose values are      |
|                       |determined by the selection for {model_name}.                                        |
+-----------------------+-------------------------------------------------------------------------------------+

.. figure:: /figures/454_goma_physics.png
	:align: center
	:width: 90%

Vapor pressure model choices and their parameters are discussed below.

*Models in the first class...*

+-----------------------+------------------------------------------------------------------------------------------+
|**KELVIN** <species>   |The <float_list> for the **KELVIN** option specifies input values for seven parameters:   |
|<float_list>           |                                                                                          |
|                       | * <float1> - Equilibrium vapor pressure across a flat interface                          |
|                       | * <float2> - Liquid density                                                              |
|                       | * <float3> - Molecular weight of the liquid                                              |
|                       | * <float4> - Gas law constant                                                            |
|                       | * <float5> - Operating temperature                                                       |
|                       | * <float6> - Molecular weight of air or gas phase                                        |
|                       | * <float7> - Ambient pressure of that gas phase                                          |
|                       |                                                                                          |
|                       |The **KELVIN** option is used to include the effect of vapor-pressure lowering that       |
|                       |results in equilibrium over high curvature menisci, i.e., small pores. The equation form  |
|                       |of this is                                                                                |
+-----------------------+------------------------------------------------------------------------------------------+
|**FLAT** <species>     |The **FLAT** option requires the same seven parameters as the **KELVIN** model but leaves |
|<float_list>           |off the exponential function, i.e., the vapor pressure is independent of the level of     |
|                       |capillary pressure. The constants are still needed so that the gas-phase concentration can|
|                       |be calculated with the ideal gas law. See the **KELVIN** option above for definition of   |
|                       |the <float_list> values.                                                                  |
+-----------------------+------------------------------------------------------------------------------------------+
|**IDEAL_GAS** <species>|The <float_list> for this model has three values, where:                                  |
|<float_list>           |                                                                                          |
|                       | * <float1> - Molecular weight of the gas                                                 |
|                       | * <float2> - Gas law constant                                                            |
|                       | * <float3> - Operating temperature                                                       |
+-----------------------+------------------------------------------------------------------------------------------+

*Models in the second class..*

+-----------------------+------------------------------------------------------------------------------------------+
|**CONSTANT** <species> |This model is used for a constant species source such as a homogeneous reaction term. The |
|<float1>               |<float_list> has a single value:                                                          |
|                       |                                                                                          |
|                       | * <float1> - Vapor pressure                                                              |
+-----------------------+------------------------------------------------------------------------------------------+
|**ANTOINE** <species>  |The **ANTOINE** model for vapor pressure is used in conjunction with the *VL_EQUIL*       |
|<float_list>           |boundary condition. If specified, a temperature-dependent vapor pressure for species i is |
|                       |calculated.                                                                               |
|                       |                                                                                          |
|                       |The model requires six values in the <float_list>, where:                                 |
|                       |                                                                                          |
|                       | * <float1> - A, the unit conversion factor for pressure based on the units in the        |
|                       |   material file                                                                          |
|                       | * <float2> - Bi, Antoine coefficient for species i                                       |
|                       | * <float3> - Ci, Antoine coefficient for species i                                       |
|                       | * <float4> - Di, Antoine coefficient for species i                                       |
|                       | * <float5> - Tmin, Minimum temperature of the range over which the Antoine relation will |
|                       |   hold                                                                                   |
|                       | * <float6> - Tmax, Maximum temperature of the range over which the Antoine relation will |
|                       |   hold                                                                                   |
+-----------------------+------------------------------------------------------------------------------------------+
|**RIEDEL** <species>   |The **RIEDEL** model for vapor pressure is used in conjunction with the *VL_EQUIL*        |
|<float_list>           |boundary condition card. If specified, a temperature-dependent vapor pressure for species |
|                       |i is calculated.                                                                          |
|                       |                                                                                          |
|                       |The model requires eight values in the <float_list>, where:                               |
|                       |                                                                                          |
|                       | * <float1> - A, the unit conversion factor for pressure based on the units in the        |
|                       |   material file                                                                          |
|                       | * <float2> - Bi, Riedel constant for species i                                           |
|                       | * <float3> - Ci, Riedel constant for species i                                           |
|                       | * <float4> - Di, Riedel constant for species i                                           |
|                       | * <float5> - Ei, Riedel constant for species i                                           |
|                       | * <float6> - Fi, Riedel constant for species i                                           |
|                       | * <float7> - Tmin, Minimum temperature of the range over which the relation will hold    |
|                       | * <float8> - Tmax, Maximum temperature of the range over which the relation will hold    |
+-----------------------+------------------------------------------------------------------------------------------+

.. figure:: /figures/455_goma_physics.png
	:align: center
	:width: 90%

.. figure:: /figures/456_goma_physics.png
	:align: center
	:width: 90%

------------
**Examples**
------------

An example use of the Antoine model for vapor pressure follows:

::

   Vapor Pressure = ANTOINE 0 1 9.380340229 3096.516433 -53.668 0.1 1000

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
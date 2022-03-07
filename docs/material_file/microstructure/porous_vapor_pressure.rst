*************************
**Porous Vapor Pressure**
*************************

::

   Porous Vapor Pressure = {model_name} {integer} {float_list} [M/L-t2]

-----------------------
**Description / Usage**
-----------------------

Used to specify the model for the vapor pressure for each multiphase flow component
in the porous medium that is activated for *Media Type* **POROUS_UNSATURATED** or
**POROUS_TWO_PHASE**.

Definitions of the input parameters are as follows:

+-------------------+-------------------------------------------------------------------------------------+
|{model_name}       |The permissible values for the model in this class are **KELVIN** and **FLAT** for a |
|                   |volatile liquid, and **NON_VOLATILE** for a non-volatile liquid.                     |
+-------------------+-------------------------------------------------------------------------------------+
|{integer}          |All models require an integer field after the model name which is the species_number;|
|                   |always set to zero until a multicomponent capability exists.                         |
+-------------------+-------------------------------------------------------------------------------------+
|{float_list}       |One or more floating point numbers (<float1> through <floatn>) whose values are      |
|                   |determined by the selection for {model_name}.                                        |
+-------------------+-------------------------------------------------------------------------------------+

Porous vapor pressure model choices and their parameters are presented below; consult
the Technical Discussion for relevant details.

+-----------------------------+-------------------------------------------------------------------------------------+
|**KELVIN** <integer> <float1>|For the **KELVIN** porous vapor pressure model, the {float_list} has a five values:  |
|<float2>... <float5>         |                                                                                     |
|                             | * <float1> - p*v, vapor pressure on a flat interface                                |
|                             | * <float2> - ρl, the liquid density                                                 |
|                             | * <float3> - Mw , molecular weight of liquid                                        |
|                             | * <float4> - R, the gas law constant                                                |
|                             | * <float5> - T, the operating temperature                                           |
+-----------------------------+-------------------------------------------------------------------------------------+
|**FLAT** <integer> <float1>  |For the **FLAT** porous vapor pressure model, the {float_list} has a five values     |
|<float2>... <float5>         |(same as **KELVIN** above):                                                          |
|                             |                                                                                     |
|                             | * <float1> - p*v, vapor pressure on a flat interface                                |
|                             | * <float2> - ρl, the liquid density                                                 |
|                             | * <float3> - Mw , molecular weight of liquid                                        |
|                             | * <float4> - R, the gas law constant                                                |
|                             | * <float5> - T, the operating temperature                                           |
|                             |                                                                                     |
|                             |The **FLAT** option requires the same parameters as the **KELVIN** model but leaves  |
|                             |out the exponential function.                                                        |
+-----------------------------+-------------------------------------------------------------------------------------+
|**NON_VOLATILE** <integer>   |The **NON_VOLATILE** model requires no additional input.                             |
+-----------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The sample input card:

::

   Porous Vapor Pressure = FLAT 0 {Vap_Pres} {density} {30.} {Rgas} {T}

applies the **FLAT** model as described above to vapor-liquid equilibrium (assumed to be
single component for now) using all APREPRO-defined parameters.

-------------------------
**Technical Discussion**
-------------------------

The **KELVIN** option is used to include the effect of vapor-pressure lowering that
results in equilibrium over high curvature menisci, i.e., small pores. The equation form
of this is:

.. figure:: /figures/421_goma_physics.png
	:align: center
	:width: 90%

The **FLAT** option requires the same parameters but leaves out the exponential function.
The constants are still needed so that the gas-phase concentration can be calculated
with the ideal gas law. The functional form is

.. figure:: /figures/422_goma_physics.png
	:align: center
	:width: 90%

where S is the local saturation, and ρgv is the gas phase density of vapor. This model is
ad-hoc but nonetheless leads to some interesting results. It basically says that as
saturation increases, the gas-liquid menisci, and correspondingly the interfacial area
available for evaporation, become more concentrated and hence the gas-phase vapor
concentration increases.

The **NON_VOLATILE** option should be set if no gas-phase transport of vapor of the
liquid phase component is desired, as if the liquid phase were non-volatile. Goma, with
this choice, sets the gas phase concentration of liquid vapor to zero.

For nonvolatile pore liquids, the vapor pressure on a flat interface, viz. the first required
floating point on this card, should be set to zero. As of 6/13/02 this card has only been
implemented for pure liquid solvents, so that no equilibrium solvent partitioning across
the interface is present.


--------
**FAQs**
--------

Sometimes system aborts can happen with the Kelvin model because of real large,
negative capillary pressures. In this case, the exponential term can exceed the machine
limit. This can happen well into a transient run. The user should be aware of this;
consult GT-009.3 for tips related to dealing with this problem.

--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

SAND96-2149: Drying in Deformable Partially-Saturated Porous Media: Sol-Gel
Coatings, Cairncross, R. A., P. R. Schunk, K. S. Chen, S. S. Prakash, J. Samuel, A. J.
Hurd and C. Brinker (September 1996)

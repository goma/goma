************************
**Porous Gas Constants**
************************

::

   Porous Gas Constants = IDEAL_GAS <float_list> [varies]

-----------------------
**Description / Usage**
-----------------------

This required card is used for *Media Types* of **POROUS_UNSATURATED** and
**POROUS_TWO_PHASE**, and is used to input some standard thermodynamic gas
constants needed for vapor-liquid equilibrium calculations (see *Media Type* card).
Eventually more than one model may be allowed for nonideal gas situations.

The **IDEAL_GAS** model is the only model currently requiring standard constants; they
are defined as follows:

+-------------------+-------------------------------------------------------------------------------------+
|**IDEAL_GAS**      |the model name requiring constants for the thermodynamic ideal gas law.              |
|                   |                                                                                     |
|                   | * <float1> - MWair, the molecular weight of the insoluble                           |
|                   |   gas in the gas phase [g/mole].                                                    |
|                   | * <float2> - R, the universal gas law constant [M-L2/t2/K]                          |
|                   | * <float3> - T, the temperature [deg K]                                             |
|                   | * <float4> - pamb, the ambient gas pressure.                                        |
+-------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The sample input card follows:

::

   Porous Gas Constants = IDEAL_GAS 28.0 8. 315 275 1.06e+5

-------------------------
**Technical Discussion**
-------------------------

For **POROUS_UNSATURATED** media types the ambient pressure dictates the
equilibrium pressure for the calculation of the gas-phase density of solvent (viz. the
total ambient pressure minus the vapor pressure will be the gas partial pressure, from
which the concentration of gas can be computed based on the other gas constants). In
**POROUS_TWO_PHASE** media types, the gas partial pressure is a dependent
variable and computed as a part of the Darcy law mass balance. In this case the
dynamic pressure is used instead of <float4> here for the calculation of the gas-phase
concentrations.

It is important to realize that setting the ambient pressure on this card for Media Types
of **POROUS_UNSATURATED** will potentially affect your saturation curve and the
appropriate values of your liquid phase pressure boundary conditions. If possible, you
should set this value to zero, and base your Saturation versus vapor pressure curve
accordingly. Also, in that case your liquid pressure boundary conditions can all be
referenced to zero. However, if you choose a gauge pressure, or thermodynamic
pressure, you Saturation/capillary pressure curve must be shifted accordingly, as do
your boundary conditions. Also, remember these pressures will affect your solid
pressure state in poroelastic problems.



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
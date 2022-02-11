*********************
**YFLUX_SULFIDATION**
*********************

::

	BC = YFLUX_SULFIDATION SS <bc_id> {char_string} <integer> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

The *YFLUX_SULFIDATION* card enables computation of the molar flux of the
diffusing species (e.g. copper vacancy) using copper-sulfidation kinetics at the
specified boundary (gas/ :math:`Cu^2S` or Cu/ :math:`Cu^2S` interface). When used in conjunction with
the *KIN_LEAK* card, it also enables the determination of velocity normal to the moving
gas/ :math:`Cu^2S` interface.

The <float_list> contains ten values to be defined; these and all input parameter
definitions are as follows:

+---------------------+-----------------------------------------------------+
|**YFLUX_SULFIDATION**| Name of the boundary condition (<bc_name>).         |
+---------------------+-----------------------------------------------------+
|**SS**               | Type of boundary condition (<bc_type>), where **SS**|
|                     | denotes side set in the EXODUS II database.         |
+---------------------+-----------------------------------------------------+
|<bc_id>              | The boundary flag identifier, an integer associated |
|                     | with <bc_type> that identifies the boundary location| 
|                     | (side set in EXODUS II) in the problem domain.      |
+---------------------+-----------------------------------------------------+
|{char_string}        | Name of sulfidation kinetic models. Allowable names |
|                     | are:                                                |
|                     |                                                     |
|                     |   * **SOLID_DIFFUSION_SIMPLIFIED**                  |
|                     |   * **SOLID_DIFFUSION**                             |
|                     |   * **SOLID_DIFFUSION_ELECTRONEUTRALITY**           |
|                     |   * **SOLID_DIFFUSION_ELECTRONEUTRALITY_LINEAR**    |
|                     |   * **GAS_DIFFUSION**                               |
|                     |   * **FULL**                                        |
|                     |   * **ANNIHILATION**                                |
|                     |   * **ANNIHILATION_ELECTRONEUTRALITY**              |
|                     |                                                     |
|                     | Detailed description of kinetic models with these   |
|                     | name key words are presented in the Technical       |
|                     | Discussion section below.                           |
+---------------------+-----------------------------------------------------+
|<integer>            | Species number of concentration.                    |
+---------------------+-----------------------------------------------------+
|<float1>             | Stoichiometric coefficient.                         |
+---------------------+-----------------------------------------------------+
|<float2>             | Rate constant for forward copper sulfidation        |
|                     | reaction.                                           |
+---------------------+-----------------------------------------------------+
|<float3>             | Activation energy for forward copper sulfidation    |
|                     | reaction.                                           |
+---------------------+-----------------------------------------------------+
|<float4>             | Rate constant for backward copper sulfidation       |
|                     | reaction.                                           |
+---------------------+-----------------------------------------------------+
|<float5>             | Activation energy for backward copper sulfidation   |
|                     | reaction.                                           |
+---------------------+-----------------------------------------------------+
|<float6>             | Temperature.                                        |
+---------------------+-----------------------------------------------------+
|<float7>             | Bulk concentration of :math:`H^2S`                  |
+---------------------+-----------------------------------------------------+
|<float8>             | Bulk concentration of :math:`O^2`                   |
+---------------------+-----------------------------------------------------+
|<float9>             | Molecular weight of copper sulfide ( :math:`Cu_2S` )|
+---------------------+-----------------------------------------------------+

------------
**Examples**
------------

Examples of this input card follow:
::

   BC = YFLUX_SULFIDATION SS 3 SOLID_DIFFUSION_ELECTRONEUTRALITY 0 -2.0
   1.46e+7 6300.0 1.2e+14 6300.0 303.0 1.61e-11 8.4e-6 159.14 5.6

::

   BC = YFLUX_SULFIDATION SS 1 ANNIHILATION_ELECTRONEUTRALITY 0 1.0   10.0
   0.0   0.0   0.0   303.0 1.61e-11 8.4e-6 159.14 5.6

::

   BC = YFLUX_SULFIDATION SS 3 SOLID_DIFFUSION 1 -2.0 1.46e7 6300.0 1.2e+14
   6300.0 303.0 1.61e-11 8.4e-6 159.14 5.6

::

   BC = YFLUX_SULFIDATION SS 1 ANNIHILATION 1 1.0 10.0   0.0   0.0   0.0
   303.0 1.61e-11 8.4e-6 159.14 5.6

-------------------------
**Technical Discussion**
-------------------------

Key word **SOLID_DIFFUSION_SIMPLIFIED** refers to the following simplified
kinetic model of copper sulfidation in which gas-phase diffusion is neglected and Cu 
is taken to be the diffusing species:

.. figure:: /figures/143_goma_physics.png
	:align: center
	:width: 90%

where *r* is molar rate of formation of sulfidation-corrosion product, :math:`Cu^2S`, per unit area, :math:`^cH_2S`
is the molar concentration of :math:`H_2S` taken to be fixed at its bulk value,
:math:`^cCu` is the
molar concentration of Cu at the sulfidation surface ( :math:`Cu_2S` /gas interface), *k* is the rate
constant, *E* is the activation energy, *R* is the universal gas constant, and *T* is the temperature.

Key word **SOLID_DIFFUSION** refers to the following kinetic model of copper
sulfidation in which gas-phase diffusion is neglected and Cu vacancies and electron
holes are taken as the diffusing species:

.. figure:: /figures/144_goma_physics.png
	:align: center
	:width: 90%

where *r* is molar rate of formation of :math:`Cu_2S` per unit area, :math:`^cH_2S` and :math:`^cO_2` are the molar
concentrations of :math:`H_2S` and :math:`O_2`, respectively, taken to be fixed at their bulk values, :math:`^c_V`
and :math:`^c_h` are the molar concentrations of Cu vacancies and electron holes, respectively,
at the sulfidation surface, :math:`k_1` and :math:`k_{-1}` are rate constants, respectively, for the forward and
backward sulfidation reactions, :math:`E_1` and :math:`E_{-1}` are activation energies, respectively, for the
forward and the backward sulfidation reactions.

Key word **SOLID_DIFFUSION_ELECTRONEUTRALITY** refers to the following
kinetic model of copper sulfidation in which Cu vacancies and electron holes are taken
as the diffusing species and the electroneutrality approximation is applied such that
concentrations of Cu vacancies and electron holes are equal to each other:

.. figure:: /figures/145_goma_physics.png
	:align: center
	:width: 90%

Key word **SOLID_DIFFUSION_ELECTRONEUTRALITY_LINEAR** refers to the
following kinetic model of copper sulfidation:

.. figure:: /figures/146_goma_physics.png
	:align: center
	:width: 90%

Key word **GAS_DIFFUSION** refers to the following simplified kinetic model of
copper sulfidation in which solid-phase diffusion is neglected, and :math:`H_2S` and 
:math:`O_2` are taken to be the diffusing species:

.. figure:: /figures/147_goma_physics.png
	:align: center
	:width: 90%

Key word **FULL** refers to the following kinetic model in which diffusion in both the
gas phase and the solid phase are important, and :math:`H_2S`, :math:`O_2`,
Cu vacancies, and electron holes are taken as the diffusing species:

.. figure:: /figures/148_goma_physics.png
	:align: center
	:width: 90%

where :math:`^cH_2S` and :math:`^cO_2` are the time-dependent molar concentrations of 
:math:`H_2S` and :math:`O_2`, respectively, at the sulfidation surface.

Key word **ANNIHILATION** refers to the following kinetic model in which diffusion
in both the gas phase and the solid phase are important, and 
:math:`H_2S`, :math:`O_2`, Cu vacancies, and
electron holes are taken as the diffusing species:

.. figure:: /figures/149_goma_physics.png
	:align: center
	:width: 90%

where :math:`k_2` are :math:`E_2` are the rate constant and activation energy, 
respectively, for the annihilation reaction.

Key word **ANNIHILATION_ELECTRONEUTRALITY** is similar to
**ANNIHILATION** except that, here, the electroneutrality approximation is applied and
concentrations of Cu vacancies and electron holes are taken to be equal to each other:

.. figure:: /figures/150_goma_physics.png
	:align: center
	:width: 90%




.. TODO - Lines 108, 122, 138, 145, 153, 161, 173, and 184 is a photo that needs to be replaced with the correct equation.
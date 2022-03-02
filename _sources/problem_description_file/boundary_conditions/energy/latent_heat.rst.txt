***************
**LATENT_HEAT**
***************

::

	BC = LATENT_HEAT SS <bc_id> <integer> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition card is used for latent heat release/adsorption at an external
interface. The flux quantity is specified on a per mass basis so the heat and mass
transfer coefficients are in units of L/t.

The <float_list> has three values to be specified; definitions of the input parameters 
are as follows:

=============== =================================================================
**LATENT_HEAT** Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS**
                denotes side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set
                in EXODUS II) in the problem domain.
<integer>       Species number.
<float1>        Latent Heat for the pure w+1 species component.
<float2>        Mass transfer coefficient for the w+1 species component.
<float3>        Sink concentration for the w+1 species component.
=============== =================================================================

The float values on this card apply to the bulk species (i.e., the w+1 component) in a
multi-species problem and in the case of a pure fluid. Important usage comments are
contained in the Technical Discussion below.

------------
**Examples**
------------

The following is a sample input card:
::

   BC = LATENT_HEAT SS 3 0   540. 0.1   0.

Two more detailed examples are contained in the Technical Discussion section.

-------------------------
**Technical Discussion**
-------------------------

The *LATENT_HEAT* boundary condition has the form

.. figure:: /figures/135_goma_physics.png
	:align: center
	:width: 90%

where :math:`\underline{n}` is the outward normal to the surface, :math:`\underline{q}` is the heat flux vector, :math:`\Delta` :math:`H_\nu` is the 
heat of vaporization, :math:`\rho` is density, :math:`h_i` is the heat-transfer coefficient for species *i*, and :math:`y^0_i`
is the reference concentration of species *i* at locations remote from the boundary;
the summation is over the number of species in the material block. The manner of
usage of this boundary condition depends on the set of conditions characterizing the
problem; example conditions are described below.

This card is used for external surfaces for which heat transfer and mass transfer beyond
it’s surfaces are governed by heat and mass transfer coefficients. The *LATENT_HEAT*
BC is applied to the energy equation so a heat flux can be specified for thermal
problems alone. The mass transfer portion of the vaporization phenomenon is handled
by the *KIN_LEAK* and *YFLUX BC* cards; these boundary conditions are applied to the
mesh equations. The *LATENT_HEAT_INTERNAL* card should be used for internal
surfaces, or interfaces, at which transfer is governed by actual physics being modeled
as a part of the problem.

When vaporization of a pure liquid is being modeled, there is only a ’single species’,
the bulk volatile liquid. In the single species case, the *Species Properties* of the
corresponding material file (which includes the *Heat of Vaporization* card) is not even read so the actual value of the latent heat of vaporization must be entered on the*LATENT_HEAT* card (<float1>). If multiple species are present, the latent heat value
for each species is entered in the material file and the *LATENT_HEAT* card does for the energy equation the same thing the *KIN_LEAK* card does for the mesh equation (i.e., collects the flux conditions that apply for each species).

For mass transfer in the single species/pure liquid case, the mass transfer coefficient is specified on the *KIN_LEAK* card. When multiple species are present, the mass transfer coefficient and driving concentration on the *KIN_LEAK* card are set to zero and the appropriate coefficient and driving concentration are set for each species on the *YFLUX* card, one for each species. The *KIN_LEAK* card (or the LATENT HEAT for energy flux) must be present to signal *Goma* to look for multiple *YFLUX* cards.

The latent heat quantity is specified on a per mass basis and the transfer coefficients are in units of L/t. Some examples of *LATENT_HEAT* application follow:

Pure Liquid Case
::

   BC = LATENT HEAT   SS 3 0   540.   0.1   0.

::

   BC = KIN_LEAK   SS 3   0.1   0.

Two-Species Case
::

   BC = LATENT HEAT   SS 3 0   0.   0.1   0.

::

   BC = KIN_LEAK   SS 3   0.   0.

::

   BC = YFLUX SS 3 0   0.12   0.

::

   BC = YFLUX   SS 3 1   0.05   0.

plus, in the corresponding material file:

::

   ---Species Properties
   Diffusion Constitutive Equation= FICKIAN
   Diffusivity = CONSTANT   0   1.e-8
   Latent Heat Vaporization = CONSTANT   0   540.
   Latent Heat Fusion = CONSTANT   0   0.
   Vapor Pressure = CONSTANT   0   0.
   Species Volume Expansion = CONSTANT   0   1.
   Reference Concentration = CONSTANT   0   0.
   Diffusivity = CONSTANT   1   1.e-6
   Latent Heat Vaporization = CONSTANT   1   125.
   Latent Heat Fusion = CONSTANT   1   0.
   Vapor Pressure = CONSTANT   1   0.
   Species Volume Expansion = CONSTANT   1   1.
   Reference Concentration = CONSTANT   1   0.



--------------
**References**
--------------

No References.
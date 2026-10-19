********************************************
**Gas Phase Darcy Velocity in Porous Media**
********************************************

::

   Gas phase Darcy velocity in porous media = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This post-processing option will lead to the explicit calculation and storage of the
Darcy velocity components in the gas phase, viz. the velocity of the gas phase due to
gas-phase pressure gradients. This option is only available for
*POROUS_TWO_PHASE* media types (cf. *Media Type* card). The velocity components
appear in the output EXODUS II file as the nodal variables **Darcy_Vel_g_0, Darcy_Vel_g_1** and 
**Darcy_Vel_g_2**.

The permissible values for this postprocessing option are:

============= ================================================================
**yes**       Calculate the gas-phase Darcy velocity components and
              write to the output EXODUSII file.
**no**        Do not calculate the gas phase velocity components.
============= ================================================================

------------
**Examples**
------------

This input example turns on calculation of the gas phase velocity components:
::

   Gas phase Darcy velocity in porous media =yes

-------------------------
**Technical Discussion**
-------------------------

The gas-phase Darcy velocity is given by the extended Darcy law, which accounts for
the relative reduced flow due to the presence of another phase, viz.

.. math::

   \boldsymbol{F}_g=\rho_g\boldsymbol{v}_g=-\frac{\rho_gkk_g}{\mu_g}(\nabla p_g-\rho_g\boldsymbol{g})

Here :math:`\nu_g` represents the Darcy flux, or Darcy velocity, in the gas phase, k is the
permeability of the porous medium, :math:`k_g` is the relative permeabilities for the gas and
liquid phases respectively, :math:`\mu_g` are the gas viscosity, :math:`p_g` is the pressure in the gas phase, and g is the gravitational force vector. :math:`\rho_g` is the density of the gas phase and is equal to the sum of the partial densities of air and solvent vapor,

.. math::

   \rho_g=\rho_{gv}+\rho_{ga}



--------------
**References**
--------------

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

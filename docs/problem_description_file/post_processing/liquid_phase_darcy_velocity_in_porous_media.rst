***********************************************
**Liquid Phase Darcy Velocity in Porous Media**
***********************************************

::

   Liquid phase Darcy velocity in porous media = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This post-processing option will lead to the explicit calculation and storage of the
Darcy velocity components in the liquid phase, viz. the velocity of the liquid phase due
to liquid phase pressure gradients. This option is available for all porous media types
(cf. *Media Type* card). The velocity components appear in the output EXODUS II file
as the nodal variables **Darcy_Vel_l_0, Darcy_Vel_l_1** and **Darcy_Vel_l_2**.

The permissible values for this postprocessing option are:

============= ================================================================
**yes**       Calculate the liquid-phase Darcy velocity components and
              write to the output EXODUSII file.
**no**        Do not calculate the liquid phase velocity components.
============= ================================================================

------------
**Examples**
------------

This input example turns on calculation of the liquid phase velocity components:
::

   Liquid phase Darcy velocity in porous media = yes

-------------------------
**Technical Discussion**
-------------------------

The liquid-phase Darcy velocity is given by the extended Darcy law, which accounts
for the relative reduced flow due to the presence of another phase, viz.

.. figure:: /figures/332_goma_physics.png
	:align: center
	:width: 90%

Here :math:`v_l` represents the Darcy flux, or Darcy velocity, in the gas phase, k is the
permeability of the porous medium, :math:`k_l` is the relative permeabilities for the liquid and
liquid phases respectively, :math:`\mu_l` are the liquid viscosity, :math:`p_l` is the pressure in the liquid phase, and g is the gravitational force vector. :math:`\rho_l` is the density of the liquid phase.



--------------
**References**
--------------

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

.. 
	TODO - Line 43 is a photo that needs to be swapped with the equation.
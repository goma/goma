*********************************************
**Total Density of Solvents in Porous Media**
*********************************************

::

   Total density of solvents in porous media = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This post processing option can be used to trigger the computation and output of the
total density of solvents as a nodal field variable to the output EXODUS II file. Three
nodal variables are written, **Rho_Total_Liq, Rho_Total_air** and **Rho_Total_solid**.
The mathematical details are given below in the technical discussion. This option
applies to media types of *POROUS_SATURATED, POROUS_UNSATURATED*, and
*POROUS_TWO_PHASE* (see *Media Type* card). The options are:

============= ================================================================
**yes**       Calculate and write the total solvent densities as a
              postprocessing variable to the output EXODUSII file.
**no**        Do not calculate the total solvent densities.
============= ================================================================

------------
**Examples**
------------

::

   Total density of solvents in porous media = yes

This card will result in the calculation and output of the mixture density of solvent
(viz., phase mixture of liquid solvent in vapor form, liquid form, and the form adsorbed
in the solid skeleton for partially saturated porous flows). The form of that mixture
density is given in the technical discussion.

-------------------------
**Technical Discussion**
-------------------------

In saturated flow cases, viz. for *Media Type* selection *POROUS_SATURATED*, the total
solvent density is

.. figure:: /figures/325_goma_physics.png
	:align: center
	:width: 90%

where :math:`\rho_1` is the pure liquid density and :math:`\phi` is the porosity. Here we have assumed that no liquid solvent is adsorbed into the solid struts (currently the assumption used
throughout *Goma*).

For partially saturated flows, viz. for *Media Type* selection *POROUS_UNSATURATED*
or *POROUS_TWO_PHASE*, the total density is given by

.. figure:: /figures/326_goma_physics.png
	:align: center
	:width: 90%

where :math:`\rho_{gv}` is the density of solvent vapor in the total gas-solvent vapor mixture (see
*Density of solvents in gas phase in porous media card*), *S* is the saturation (see *Porous
Saturation* card), and :math:`\chi_{ls}` is the volume fraction of solvent in liquid phase (including any condensed species component). When calculating the total density of the liquid
(**Rho_Total_liq**), the liquid vapor density comes from a Kelvin vapor-liquid
equilibrium relation. The total density of the gas phase (**Rho_Total_gas**) will use a
vapor density fro air and a volume fraction of zero (0) since air is insoluble. The total
density of the solid in the gas (**Rho_Total_solid**) is zero (0).



--------------
**References**
--------------

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

.. 
	TODO - Lines 46 and 56 are photos that need to be swapped with the equations.
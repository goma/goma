****************************************************
**Density of Solvents in Gas Phase in Porous Media**
****************************************************

::

   Density of solvents in gas phase in porous media = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This post processing option can be used to trigger the computation and output of the
density of solvents in the gas phase only, including the volume occupied by the
assumed insoluble gas (e.g. air), as a nodal field variable to the output EXODUS II file.
The nodal variables are called **RhoSolv_g_liq, RhoSolv_g_air** and **RhoSolv_g_solid**.
The mathematical details are given below in the technical discussion. This option
applies to media types of *POROUS_UNSATURATED*, and *POROUS_TWO_PHASE*
(see *Media Type* card). The options are:

============= ================================================================
**yes**       Calculate and write the gas phase solvent density as a
              postprocessing variable to the output EXODUSII file.
**no**        Do not calculate the total solvent density.
============= ================================================================

------------
**Examples**
------------

The following input card turns off writing solvent densities to the EXODUS II file:
::

   Density of solvents in gas phase in porous media = no

-------------------------
**Technical Discussion**
-------------------------

The air and solid components are insoluble in the gas phase so the **RhoSolv_g_air** and
**RhoSolv_g_solid** variables will be zero. The gas-density of liquid solvents
**(RhoSolv_g_liq)** is determined from the vapor-liquid equilibrium relationship at a
liquid-vapor meniscus. Specifically,

.. figure:: /figures/327_goma_physics.png
	:align: center
	:width: 90%

where :math:`M_w` is the average molecular weight of solvents in the mixture, R is the ideal gas
constant, T is the temperature, and :math:`\rho_v` is the equilibrium vapor pressure. Note that this
vapor pressure can be affected by local meniscus curvature through the Kelvin equation
(cf. Schunk, 2002 and *Porous Vapor Pressure* card).



--------------
**References**
--------------

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

..
	TODo - Line 45 is a photo that needs to be exchanged with the correct equation.
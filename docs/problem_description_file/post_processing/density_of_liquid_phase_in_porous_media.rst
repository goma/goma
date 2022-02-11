*******************************************
**Density of Liquid Phase in Porous Media**
*******************************************

::

   Density of liquid phase in porous media = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This post processing option can be used to trigger the computation and output of the
density of solvents in the liquid phase only, averaged over the mixture, as a nodal field
variable to the output EXODUS II file. The nodal variable is called **Rho_Liq_Phase**.
The mathematical details are given below in the technical discussion. This option
applies to media types of *POROUS_SATURATED, POROUS_UNSATURATED*, and
*POROUS_TWO_PHASE* (see *Media Type* card). The options are:

============= ================================================================
**yes**       Calculate and write the density of solvent in the liquid phase
              as a postprocessing variable to the output EXODUSII file.
**no**        Do not calculate the liquid solvent density.
============= ================================================================

------------
**Examples**
------------

An example of an input card which activates writing of the density to the EXODUS II
file is:
::

   Density of liquid phase in porous media = yes

-------------------------
**Technical Discussion**
-------------------------

In liquid-saturated flow cases, viz. for *Media Type* selection *POROUS_SATURATED*,
the total solvent density in the liquid phase is

.. figure:: /figures/328_goma_physics.png
	:align: center
	:width: 90%

where :math:`\rho_l` is the pure liquid density and :math:`\phi` is the porosity. Here we have assumed that no liquid solvent is adsorbed into the solid struts (currently the assumption used
throughout *Goma*).

For partially saturated flows, viz. for *Media Type* selection *POROUS_UNSATURATED*
or *POROUS_TWO_PHASE*, the density of solvent in the liquid phase only is given by

.. figure:: /figures/329_goma_physics.png
	:align: center
	:width: 90%

where S is the saturation (see *Porous Saturation* card). Compare this with the quantity
computed with the *Total density of solvents in porous media* card.



--------------
**References**
--------------

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

.. 
	TODO - Lines 43 and 53 are photos that need to be swapped with the equations.
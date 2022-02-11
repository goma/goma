************
**VL_EQUIL**
************

::

	BC = VL_EQUIL SS <bc_id> <integer_list> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC/MASS)**

This boundary condition card enforces vapor-liquid equilibrium between a gas phase
and a liquid phase using Raoult’s law. The condition only applies to interphase mass,
heat, and momentum transfer problems with discontinuous (or multivalued) variables
at an interface, and it must be invoked on fields that employ the **Q1_D** or **Q2_D**
interpolation functions to “tie” together or constrain the extra degrees of freedom at 
the interface in question.

The <integer_list> has three values and the <float_list> has five values; definitions 
of the input parameters are as follows:

============ =============================================================
**VL_EQUIL** Name of the boundary condition (<bc_name>).
**SS**       Type of boundary condition (<bc_type>), where **SS** denotes
             side set in the EXODUS II database.
<bc_id>      The boundary flag identifier, an integer associated with
             <bc_type> that identifies the boundary location (side set in
             EXODUS II) in the problem domain.
<integer1>   Species number.
<integer2>   Element block ID of liquid phase.
<integer3>   Element block ID of gas phase.
<float1>     Base ambient pressure in gas phase.
<float2>     Molecular weight of first volatile species.
<float3>     Molecular weight of second volatile species.
<float4>     Molecular weight of condensed phase.
<float5>     Molecular weight of insoluble gas phase.
============ =============================================================

This boundary condition is applied to ternary, two-phase flows that have rapid mass
exchange between phases, rapid enough to induce a diffusion velocity at the interface,
and to thermal contact resistance type problems. The best example of this is rapid
evaporation of a liquid component into a gas. In the current discontinuous mass 
transfer model, we must require the same number of components on either side of 
interface. In this particular boundary, two of three components are considered 
*volatile*, so they participate in both vapor and liquid phases. The third component 
is  considered either non-volatile or non-condensable, so it remains in a single phase.

------------
**Examples**
------------

A sample input card follows for this boundary condition:
::

   BC = VL_EQUIL SS 4 0 1 2 1.e+06 28. 18. 1800. 18.

The above card demonstrates these characteristics: species number is “0”; liquid phase
block id is 1; gas phase block id is 2; ambient pressure is 1.e6 Pa; the molecular
weights of the volatile species are 28 and 18; of the condensed phase and insoluble
portion of the gas phase, 1800 and 18, respectively.

-------------------------
**Technical Discussion**
-------------------------

One of the simplest forms of the equilibrium relation is the Raoult’s law, where the
mole fraction of a species is equal to its mole fraction in the liquid multiplied by the ratio of its pure component vapor pressure to the total pressure in the system.

.. figure:: /figures/154_goma_physics.png
	:align: center
	:width: 90%

where :math:`y_i` are the mole fraction of species *i* in the gas phase and 
:math:`x_i` is the mole fraction in
the liquid phase. The molecular weights required in this boundary card are used for
converting mass fractions to mole fractions. The temperature dependency in the
equilibrium expression comes from a temperature-dependent vapor pressure model.
Either Riedel or Antoine temperature-dependent vapor pressure model can be specified
in the *VAPOR PRESSURE* material card in order to link temperature to Raoult’s law.



--------------
**References**
--------------

GTM-007.1: New Multicomponent Vapor-Liquid Equilibrium Capabilities in GOMA,
December 10, 1998, A. C. Sun

Schunk, P.R. and Rao, R.R. 1994. “Finite element analysis of multicomponent twophase
flows with interphase mass and momentum transport,” IJNMF, 18, 821-842.

.. TODO - Line 72 has a photo that needs to be replaces with the proper equation.
**********************
**IS_EQUIL_PSEUDORXN**
**********************

::

	BC = IS_EQUIL_PSEUDORXN SS <bc_id> <integer_list> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary condition card enforces equilibrium between a species component in
two ideal solution phases via a finite-rate kinetics formalism. The condition only
applies to problems containing internal interfaces with discontinuous (or multilevel)
species unknown variables. The species unknown variable must employ the **Q1_D** or
**Q2_D** interpolation functions in both adjacent element blocks. This boundary
condition constrains the species equations on both sides of the interface (i.e., 
supplies a boundary condition) by specifying the interfacial mass flux on both sides.

*IS_EQUIL_PSEUDORXN* is equivalent to the *VL_EQUIL_PSEUDORXN* except for
the fact that we do not assume that one side of the interface is a gas and the other is 
a liquid. Instead, we assume that both materials on either side of the interface are 
ideal solutions, then proceed to formulate an equilibrium expression consistent with 
that.

The <integer_list> requires three values; definitions of the input parameters are as
follows:

====================== ==============================================================
**IS_EQUIL_PSEUDORXN** Name of the boundary condition (<bc_name>).
**SS**                 Type of boundary condition (<bc_type>), where **SS**
                       denotes side set in the EXODUS II database.
<bc_id>                The boundary flag identifier, an integer associated with
                       <bc_type> that identifies the boundary location (side set 
                       in EXODUS II) in the problem domain.
<integer1>             Species number.
<integer2>             Element Block ID of the first phase, the “+” phase.
<integer3>             Element Block ID of the second phase, the “-” phase.
<float>                Rate constant for the forward reaction in units of length
                       divided by time.
====================== ==============================================================

------------
**Examples**
------------

The sample card:
::

   BC = IS_EQUIL_PSEUDORXN SS 4 0 1 2 100.

demonstrates the following characteristics: species number is “0”; the “+” phase
element block id is “1”; the “-” phase element block id is “2”; a forward rate constant of 100. cm :math:`s^{-1}`.

-------------------------
**Technical Discussion**
-------------------------

The *IS_EQUIL_PSEUDORXN* boundary condition uses the following equations
representing a kinetic approach to equilibrium expressed by an ideal solution model for
thermodynamics on either side of the interface. Initially, we relate species *k* on the + side to species *k* on the - side of the interface via a kinetic formulation, whose rate constant is fast enough to ensure equilibrium in practice. However, later we may extend the capability to more complicated stoichiometric formulations for equilibrium, since the formulation for the equilibrium expression is readily extensible, unlike *Goma’s* previous treatment.

.. figure:: /figures/163_goma_physics.png
	:align: center
	:width: 90%

.. figure:: /figures/164_goma_physics.png
	:align: center
	:width: 90%

where

.. figure:: /figures/165_goma_physics.png
	:align: center
	:width: 90%

The “-” phase is defined as the reactants, while the “+” phase is defined to be the
products. The expression for the concentration equilibrium constant, :math:`K_k^c` , is based on the ideal solution expression for the chemical potentials for species *k* in the two phases [Denbigh, p. 249],

.. figure:: /figures/166_goma_physics.png
	:align: center
	:width: 90%

where :math:`\mu_k^{+*}` (T,P) is defined as the chemical potential of species *k* in its pure state (or a
hypothetical pure state if a real pure state doesn’t exist) at temperature *T* and pressure
*P*. :math:`\mu_k^{+*}` (T,P) is related to the standard state of species *k* in phase +, :math:`\mu_k^+`, :math:`\underline{o}` + (T) , which is
independent of pressure, through specification of the pressure dependence of the pure
species k. Two pressure dependencies are initially supported:

.. figure:: /figures/167_goma_physics.png
	:align: center
	:width: 90%

.. figure:: /figures/168_goma_physics.png
	:align: center
	:width: 90%

With these definitions, :math:`K_k^c` can be seen to be equal to

.. figure:: /figures/169_goma_physics.png
	:align: center
	:width: 90%

where

.. figure:: /figures/170_goma_physics.png
	:align: center
	:width: 90%

The chemical potential for a species in a phase will be calculated either from
CHEMKIN or from the *Chemical Potential, Pure Species Chemical Potential, and
Standard State Chemical Potential cards* in the materials database file.

The choice for the independent variable for the species unknown is relatively arbitrary,although it does change the actual equation formulation for the residual and Jacobian terms arising from the boundary condition. The internal variable Species_Var_Type in the Uniform_Problem_Description structure is queried to determine what the actual species independent variable is. A choice of SPECIES_UNDEFINED_FORM is
unacceptable. If either mole fractions or molar concentration is chosen as the
independent variable in the problem, the convention has been to formulate terms of the
residuals in units of moles, cm, and seconds. Therefore, division of the equilibrium
equations :math:`W_k` by occurs before their inclusion into the residual.
:math:`J^l_k` and :math:`J^g_k` are the
diffusive flux of species *k* 
( gm :math:`cm^{-2}` 
:math:`s^{-1}` ) relative to the mass-averaged velocity. 
:math:`u^s` is
the velocity of the interface. A typical value of 
:math:`k^f` that would lead to good numerical
behavior would be 100 cm 
:math:`s^{-1}`, equivalent to a reaction with a reactive sticking
coefficient of 0.01 at 1 atm and 300 K for a molecule whose molecular weight is near to
:math:`N_2` or
:math:`H_2 S`.



--------------
**References**
--------------

Denbigh, K., The Principles of Chemical Equilibrium, Cambridge University Press,
Cambridge, 1981

.. TODO - Lines 66, 70, 76, 83, 93, 97, 103, and 109 have photos that need to be replaced with the proper equations.

**********************
**VL_EQUIL_PSEUDORXN**
**********************

::

	BC = VL_EQUIL_PSEUDORXN SS <bc_id> <integer_list> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary condition card enforces vapor-liquid equilibrium between a gas phase
and a liquid phase species component using Raoult’s law expressed via a finite-rate
kinetics formalism. The condition only applies to problems containing internal
interfaces with discontinuous (or multilevel) species unknown variables. The species
unknown variable must employ the **Q1_D** or **Q2_D** interpolation functions in both
adjacent element blocks. This boundary condition constrains the species equations on
both sides of the interface (i.e., supplies a boundary condition) by specifying the
interfacial mass flux on both sides.

Definitions of the input parameters are as follows:

====================== =============================================================
**VL_EQUIL_PSEUDORXN** Name of the boundary condition (<bc_name>).
**SS**                 Type of boundary condition (<bc_type>), where **SS**
                       denotes side set in the EXODUS II database.
<bc_id>                The boundary flag identifier, an integer associated with
                       <bc_type> that identifies the boundary location (side set
                       in EXODUS II) in the problem domain.
<integer1>             Species number.
<integer2>             Element Block ID of the liquid phase.
<integer3>             Element Block ID of the gas phase.
<float>                Rate constant for the forward reaction in units of length
                       divided by time.
====================== =============================================================

This boundary condition is typically applied to multicomponent two-phase flows that
have rapid mass exchange between phases. The best example of this is rapid
evaporation of a liquid component into a gas.

------------
**Examples**
------------

The following sample input card
::

   BC = VL_EQUIL_PSEUDORXN SS 4 0 1 2 100.

demonstrates the following characteristics: species number is “0”; liquid phase element
block id is “1”; gas phase element block id is “2”; a forward rate constant of 100.0 cm
:math:`s^{-1}`.

-------------------------
**Technical Discussion**
-------------------------

The *VL_EQUIL_PSEUDORXN* boundary condition uses the following equations
representing a kinetic approach to equilibrium expressed by Raoult’s law, relating
species *k* on the liquid side to species *k* on the gas side.

.. figure:: /figures/159_goma_physics.png
	:align: center
	:width: 90%

.. figure:: /figures/160_goma_physics.png
	:align: center
	:width: 90%

where

.. figure:: /figures/161_goma_physics.png
	:align: center
	:width: 90%

and where

.. figure:: /figures/162_goma_physics.png
	:align: center
	:width: 90%

The usage of the same index, *k*, on either side of the interface is deliberate and
represents a stoichiometric limitation to this type of boundary condition.
:math:`Y_k^l` and :math:`Y_k^g`
are the mass fraction of species *k* on the liquid and gas sides of the interface,
respectively. :math:`W_k^l` is the molecular weight of species *k*. :math:`S_k^l` is the source term for
creation of species *k* in the liquid phase at the interface (mol :math:`cm^{-2}` 
:math:`s^{-1}` ). is the pseudo
reaction rate ( :math:`cm s^{-1}` ) input from the boundary condition card. 
:math:`K_k^c` is the concentration
equilibrium constant, which for the restricted stoichiometry cases covered by this
boundary condition, is unitless. :math:`p_k^v` is the vapor pressure of gas species *k* above a liquid
entirely consisting of liquid species *k*. It is a function of temperature. 
:math:`\tilde{C}_l` is the average
concentration in the liquid (mol :math:`cm^{-3}`). :math:`C_k^l` and :math:`C_k^g` are the liquid and gas concentrations
of species *k* (mol :math:`cm^{-3}`).

The choice for the independent variable is arbitrary, although it does change the actual equation formulation for the residual and Jacobian terms arising from the boundary condition. The internal variable Species_Var_Type in the
Uniform_Problem_Description structure may be queried to determine what the
actual species independent variable is. Also note, if mole fractions or molar
concentration are chosen as the independent variable in the problem, the convention
has been to formulate terms of the residuals in units of moles, cm, and seconds.
Therefore, division of the equilibrium equations by :math:`W_k` would occur before their inclusion into the residual. :math:`j_k^l` and :math:`j_k^g` are the diffusive flux of species k (gm :math:`cm^{-2}` :math:`s^{-1}` )
relative to the mass averaged velocity. :math:`u_s` is the velocity of the interface. A typical
value of :math:`k^f` that would lead to good numerical behavior would be 100 cm 
:math:`s^{-1}`,
equivalent to a reaction with a reactive sticking coefficient of 0.01 at 1 atm and 300 K for a molecule whose molecular weight is near to :math:`N_2` or :math:`H_2S`.




.. TODO - Lines 65, 69, 75, and 81 have photos that need to be replaced with the proper equations.

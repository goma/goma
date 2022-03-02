**************
**Media Type**
**************

::

   Media Type = {model_name}

-----------------------
**Description / Usage**
-----------------------

This card is used to designate the characteristic medium type for solid materials so that
the proper microstructural features/models may be imposed. Basically, the choices are
dictated by whether the medium is to be modeled as *porous* (viz. a medium in which
flow will be determined relative to the motion of a porous solid skeleton) or as
*continuous* (viz., in which the mechanics equations apply to all parts of the medium and
not weighted by a solid fraction). If porous flow through Darcy or Brinkman
formulations are desired in the material, then the phase is designated as *continuous*.

The input parameter is a {model_name} and has the following possible values:

+-----------------------------+----------------------------------------------------------------------------+
|{model_name}                 |{model_name} Name of the media model; the choices are                       |
|                             |                                                                            |
|                             | * **CONTINUOUS**                                                           |
|                             | * **POROUS_SATURATED**                                                     |
|                             | * **POROUS_UNSATURATED**                                                   |
|                             | * **POROUS_TWO_PHASE**                                                     |
|                             | * **POROUS_BRINKMAN**                                                      |
|                             | * **POROUS_SHELL_UNSATURATED**                                             |
+-----------------------------+----------------------------------------------------------------------------+

Specific characteristics of these types are identified below, including other cards that
must be present.

* If the type chosen is **CONTINUOUS**, then the material is assumed to be
  amorphous and no further microstructure properties need to be specified (next
  required card is the *Diffusion Constitutive Equation*).

* In a porous medium with one phase in the pores (i.e. a saturated medium), use
  **POROUS_SATURATED** then only the Porosity and Permeability cards are
  required. A **POROUS_SATURATED** medium model enables the user to solve the
  simplest porous flow equation for the liquid phase pressure only for rigid porous
  media (see the porous_sat or porous_liq equation cards). For deformable porous
  saturated media, one can employ a stress balance and porosity equation for
  deformable porous media (see mesh* equation cards and porous_deform equation
  card).

* In a porous medium with two phases in the pores (such as air-water, i.e., an
  unsaturated medium), two options exist - **POROUS_UNSATURATED**, a
  formulation of the porous flow problem using the capillary pressure as the field
  variable (gas pressure assumed to be uniform), and **POROUS_TWO_PHASE**, a
  formulation of the porous flow problem using the liquid pressure and gas pressure
  as field variables. All the cards in this Microstructure porous flow section, except
  the Brinkman cards (*FlowingLiquid Viscosity and Inertia Coefficient*), are needed
  for the unsaturated or two-phase models. As in the saturated case above, these
  options can also be chosen for deformable porous media, for which the Lagrangian
  mesh stress equations and the porosity equation are used to complete the effective
  stress principle formulation.

* The **POROUS_BRINKMAN** model is an extension of the Navier-Stokes
  equation for porous media. In addition, it has an inertia term intended to account
  for boundary and interface deficiencies at Reynold’s numbers greater than one
  (Re > 1), a deficiency in all Darcy flow models (see, e.g., Gartling, et. al., 1996). It
  is a vector formulation (the momentum equations) of saturated flow in a porous
  medium which reduces to the Navier-Stokes equations as the porosity increases to
  one (φ → 1). For Brinkman flow, the input parameters (i.e., cards) that must be
  specified from this section are *Porosity, Permeability, FlowingLiquid Viscosity*,
  and *Inertia Coefficient*. Please note the use of two viscosities; for the Brinkman
  media type, the viscosity entered via the (*Mechanical Properties and Constitutive
  Equations*) Viscosity card is interpreted to be the Brinkman viscosity (μB) and is
  used to calculate the viscous stresses (see Gartling, et. al., 1996) while the
  *FlowingLiquid Viscosity* (μ) is used in the correction term for nonlinear drag
  forces in porous media. Brinkman viscosity is an effective value and can be taken
  as the porosity weighted average of the matrix and fluid. It is generally not correct
  to set it equal to the liquid viscosity (Martys, et. al., 1994; Givler and Altobelli,
  1994).

* The **POROUS_SHELL_UNSATURATED** model is used for thin shell, open
  pore, porous media, viz. the shell_sat_open equation. This media type
  instructs GOMA to obtain most of the media properties from the bulk continuum
  specifications just like **POROUS_UNSATURATED**. Exceptions are the
  Porous Shell Cross Permeability model and the Porous Shell
  Height material models. Please see the porous shell tutorial

------------
**Examples**
------------

Following is a sample card:

::

   Media Type = POROUS_TWO_PHASE

This card will require a plethora of material models for Darcy flow of liquid and gas in
a porous medium. It also will require the use of two Darcy flow mass balances in the
*Problem Description EQ* specification section, specifically *porous_liq* and *porous_gas*
equations. See references below for details.

-------------------------
**Technical Discussion**
-------------------------

In solving porous medium problems, it is important to understand that each
conservation equation represents a component, or species balance. The *porous_liq*
equation is actually a species balance for the liquid phase primary component (e.g.
water) for all phases in the medium, viz. liquid, gas, and solid. This is the case even
though the dependent variable is the liquid phase pressure. This is the only required
equation for rigid **POROUS_SATURATED** media. The same holds true for rigid
**POROUS_UNSATURATED** media, as the liquid solvent is present in liquid and gas
vapor form (it is actually taken as insoluble in the solid). For deformable media, one
must add a stress balance through the *mesh** equations (in *LAGRANGIAN* form, as
described on the *Mesh Motion* card) and a solid phase “solvent” balance which is used
to solve for the porosity, viz. the *porous_deform* equation. In these cases, the gas is
taken to be at constant pressure. If pressure driven Darcy flow is important in the gas,
an additional species balance for the primary gas component is required through the
porous_gas equation. This last case is the so-called **POROUS_TWO_PHASE** media
type.

Options for representing the solid medium as rigid or deformable are discussed under
the *Saturation, Permeability* and *Porosity* cards. When rigid porous media are
modeled, both porosity and permeability are constant. In *Goma* 4.0, these concepts
were being researched and improved, with much of the usage documentation residing
in technical memos.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s capabilities for partially saturated flow in porous media,
September 1, 2002, P. R. Schunk

Gartling D. K., C. E. Hickox and R. C. Givler 1996. "Simulations of Coupled Viscous
and Porous Flow Problems", Comp. Fluid Dynamics, 7, 23-48.

Givler, R. C. and S. A. Altobelli 1994. “A Determination of the Effective Viscosity for
the Brinkman-Forchheimer Flow Model.” J. Fluid Mechanics, 258, 355-370.

Martys, N., D. P. Bantz and E. J. Barboczi 1994. “Computer Simulation Study of the
Effective Viscosity in Brinkman’s Equation.” Phys. Fluids, 6, 1434-1439


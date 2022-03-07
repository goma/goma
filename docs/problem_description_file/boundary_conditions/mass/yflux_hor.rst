*************
**YFLUX_HOR**
*************

::

	BC = YFLUX_HOR SS <bc_id> <integer> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

The **YFLUX_HOR** card enables computation of the molar flux of the specified species
at the specified boundary (i.e., at the electrode surface) using the linearized Butler-
Volmer kinetics such as that for the hydrogen-oxidation reaction in polymerelectrolyte-
membrane fuel cells.

The <floatlist> consists of 10 values; definitions of the input parameters are as follows:

============== =================================================================
**YFLUX_HOR**  Name of the boundary condition (<bc_name>).
**SS**         Type of boundary condition (<bc_type>), where **SS** denotes
               side set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set in
               EXODUS II) in the problem domain.
<integer>      Species number of concentration.
<float1>       Product of interfacial area per unit volume by exchange
               current density, :math:`ai_0`, in units of :math:`A/cm^3`.
<float2>       Catalyst layer or catalyzed electrode thickness, 
               *H*, in unit of *cm*.
<float3>       Reference concentration, :math:`c_{ref}`, in units of 
               *moles*/ :math:`cm^3`.
<float4>       Anodic direction transfer coefficient, :math:`\alpha_a`.
<float5>       Cathodic direction transfer coefficient, :math:`\alpha_c`.
<float6>       Temperature, *T*, in unit of *K*.
<float7>       Theoretical open-circuit potential, :math:`U_0`, in unit of V.
<float8>       Reaction order, :math:`\beta`.
<float9>       Number of electrons involved in the reaction, *n*.
<float10>      Electrode potential, *V*, in unit of *V*.
============== =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = YFLUX_HOR SS 14 0 1000. 0.001 4.e-5 1. 1. 353. 0. 0.5 2. 0.

-------------------------
**Technical Discussion**
-------------------------

For electrochemical reactions such as the hydrogen-oxidation reaction (HOR), surface
overpotential is relatively small such that the Butler-Volmer kinetic model can be
linearized to yield:

.. figure:: /figures/151_goma_physics.png
	:align: center
	:width: 90%

where *r* is the surface reaction rate in units of *moles*/ :math:`cm^2-s`; 
:math:`ai_0` denotes the product of
interfacial area per unit volume by exchange current density, 
which has units of :math:`A/cm^3`;
*H* is the catalyst layer or catalyzed electrode thickness in unit of *cm*; *n* is the number of
electrons involved in the electrochemical reaction; *R* is the universal gas constant
( :math:`\equiv` 8.314 J/mole-K); T is temperature in unit of *K*; *c* and 
:math:`c_{ref}` are, respectively, species
and reference molar concentrations in units of *moles*/ :math:`cm^3`; :math:`\beta` is reaction order; :math:`\alpha_a` and :math:`\alpha_c`
are, respetively, the anodic and cathodic transfer coefficients; *V* and :math:`\phi` are,
respectively, the electrode and electrolyte potentials in unit of *V*; :math:`U_0` and is the opencircuit
potential in unit of *V*.



--------------
**References**
--------------

J. Newman, Electrochemical Systems, 2nd Edition, Prentice-Hall, NJ (1991).

K. S. Chen and M. A. Hickner, “Modeling PEM fuel cell performance using the finiteelement
method and a fully-coupled implicit solution scheme via Newton’s technique”,
in ASME Proceedings of FUELCELL2006-97032 (2006).

.. TODO - Line 62 has a photo that needs to be replaces with the proper equation.
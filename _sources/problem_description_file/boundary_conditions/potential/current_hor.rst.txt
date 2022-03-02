***************
**CURRENT_HOR**
***************

::

	BC = CURRENT_HOR SS <bc_id> <integer> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(WIC/POTENTIAL)**

The **CURRENT_HOR** card enables the specification of the variable current density as
given by linearized Butler-Volmer kinetics (such as that for the hydrogen-oxidation
reaction in polymer-electrolyte-membrane fuel cells) at the specified boundary (i.e., at
the electrode surface).

The <floatlist> consists of 9 values; definitions of the input parameters are as follows:

=============== ====================================================================
**CURRENT_HOR** Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS**
                denotes side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set
                in EXODUS II) in the problem domain.
<integer>       Species number of concentration.
<float1>        Product of interfacial area per unit volume by exchange
                current density, 
                :math:`ai_0`, in units of A/ :math:`cm^3`.
<float2>        Catalyst layer or catalyzed electrode thickness, *H*, 
                in unit of *cm*.
<float3>        Reference concentration, 
                :math:`c_{ref}`, in units of moles/ :math:`cm^3`.
<float4>        Anodic direction transfer coefficient, :math:`\alpha_a`.
<float5>        Cathodic direction transfer coefficient, :math:`\alpha_c`.
<float6>        Temperature, *T*, in unit of *K*.
<float7>        Theoretical open-circuit potential, :math:`U_0`, in unit of V.
<float8>        Reaction order, :math:`\beta`.
<float9>        Electrode potential, *V*, in unit of *V*.
=============== ====================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = CURRENT_HOR SS 14 0 1000. 0.001 4.e-5 1. 1. 353. 0. 0.5 0.

-------------------------
**Technical Discussion**
-------------------------

For electrochemical reactions such as the hydrogen-oxidation reaction (HOR), surface
overpotential is relatively small such that the Butler-Volmer kinetic model can be
linearized to yield a simplified equation for computing current density:

.. figure:: /figures/186_goma_physics.png
	:align: center
	:width: 90%

where *j* is current density in units of A/ :math:`cm^2`; :math:`ai_0` denotes the product of interfacial area
per unit volume by exchange current density, which has units of A/ :math:`cm^3`; *H* is the
catalyst layer or catalyzed electrode thickness in unit of *cm*; *c* and :math:`c_{ref}` are, respectively,
species and reference molar concentrations in units of moles/ :math:`cm^3`; :math:`\beta` is reaction order;
:math:`\alpha_a` and :math:`\alpha_c` are, respetively, the anodic and cathodic transfer coefficients; *F* is the
Faraday’s constant ( 96487 *C/mole*); *R* is the universal gasl constant ( 8.314
*J/mole-K*); *T* is temperature in unit of *K*; *V* and :math:`\phi` are, respectively, the electrode and
electrolyte potentials in unit of *V*; :math:`U_0` and is the open-circuit potential in unit of *V*.



--------------
**References**
--------------

J. Newman, Electrochemical Systems, 2nd Edition, Prentice-Hall, NJ (1991).

K. S. Chen and M. A. Hickner, “Modeling PEM fuel cell performance using the finiteelement
method and a fully-coupled implicit solution scheme via Newton’s technique”,
in ASME Proceedings of FUELCELL2006-97032 (2006).






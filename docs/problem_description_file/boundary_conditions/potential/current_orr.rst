***************
**CURRENT_ORR**
***************

::

	BC = CURRENT_ORR SS <bc_id> <integer> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(WIC/POTENTIAL)**

The **CURRENT_ORR** card enables the specification of the variable current density as
given by the Tafel kinetics (such as that for the oxygen-reduction reaction in polymerelectrolyte-
membrane fuel cells) at the specified boundary (i.e., at the electrode
surface).

The <floatlist> consists of 8 values; definitions of the input parameters are as follows:

=============== =================================================================
**CURRENT_ORR** Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS**
                denotes side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set
                in EXODUS II) in the problem domain.
<integer>       Species number of concentration.
<float1>        Product of interfacial area per unit volume by exchange
                current density, 
                :math:`ai_0`, in units of A/ :math:`cm^3`.
<float2>        Catalyst layer or catalyzed electrode thickness, *H*, in unit
                of *cm*.
<float3>        Reference concentration, :math:`c_{ref}`, in units of 
                \moles/ :math:`cm^3`.
<float4>        Cathodic direction transfer coefficient, :math:`\alpha_c`.
<float5>        Temperature, *T*, in unit of *K*.
<float6>        Electrode potential, *V*, in unit of *V*.
<float7>        Theoretical open-circuit potential, :math:`U_0`, in unit of *V*.
<float8>        Reaction order, :math:`\beta`.
=============== =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = CURRENT_ORR SS 15 1 0.01 0.001 4.e-5 1. 353. 0.7 1.18 1.

-------------------------
**Technical Discussion**
-------------------------

For electrochemical reactions such as the oxygen-reduction reaction (ORR), surface
overpotential is large and negative such that the first exponential term in the Butler-
Volmer kinetic model is much smaller than the second term and thus can be dropped to
yield the Tafel kinetic model for computing current density:

.. figure:: /figures/187_goma_physics.png
	:align: center
	:width: 90%

here *j* is current density in units of A/ :math:`cm^2`; :math:`ai_0` denotes the product of interfacial area
per unit volume by exchange current density, which has units of A/ :math:`cm^3`; *H* is the
catalyst layer or catalyzed electrode thickness in unit of *cm*; *c* and :math:`c_{ref}` are, respectively,
species and reference molar concentrations in units of moles/ :math:`cm^3`; :math:`\beta` is reaction order; :math:`\alpha_c` is the anodic and cathodic transfer coefficient; *F* is the Faraday’s constant ( 96487 *C/mole*); *R* is the universal gasl constant ( 8.314 *J/mole-K*); *T* is temperature in unit of *K*; *V* and :math:`\phi` are, respectively, the electrode and electrolyte potentials in unit of *V*; and :math:`U_0` is the open-circuit potential in unit of *V*.



--------------
**References**
--------------

J. Newman, Electrochemical Systems, 2nd Edition, Prentice-Hall, NJ (1991).

K. S. Chen and M. A. Hickner, “Modeling PEM fuel cell performance using the finiteelement
method and a fully-coupled implicit solution scheme via Newton’s technique”,
in ASME Proceedings of FUELCELL2006-97032 (2006).
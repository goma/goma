*****************
**ELEC_TRACTION**
*****************

::

	BC = ELEC_TRACTION SS <bc_id> <integer> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition card is used to add to the momentum equation the electric, or
Maxwell, stress at a free-surface. Definitions of the input parameters are as follows:

================= =============================================================
**ELEC_TRACTION** Name of the boundary condition (<bc_name>).
**SS**            Type of boundary condition (<bc_type>), where **SS**
                  denotes side set in the EXODUS II database.
<bc_id>           The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set
                  in EXODUS II) in the problem domain.
<integer>         Integer value indicating the element block ID from
                  which to apply the boundary condition.
<float>           A term-multiplier.
================= =============================================================

Since this boundary condition only adds the electric stress, it is commonly used with
one of the *CAPILLARY, CAP_RECOIL_PRESS* or *CAP_REPULSE* boundary
conditions, viz. the capillary stress must be added separately.

------------
**Examples**
------------

For a system consisting of an insulating liquid (element block ID 1) and an insulating,
passive gas (element block ID 2) with a free-surface designated by side set 12, the
following is a sample usage:
::

     BC = ELEC_TRACTION SS 12 1 1.0

::

     BC = ELEC_TRACTION SS 12 2 1.0

::

     BC = CAPILLARY SS 12 1.0 0.0 0.0 1

The first and second lines adds the electric stress due to the electric field in the liquid
and gas phases, respectively. The third line adds the capillary stress due to surface
tension. **IMPORTANT NOTE: the optional element block ID argument to the
CAPILLARY card is used to make sure that the capillary stress is added from within a
phase where the momentum equations are defined. The same holds for the KINEMATIC
boundary condition**.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition adds the electric, or Maxwell, stress contribution to the
traction condition. To use this boundary condition there must be a *VOLTAGE* equation
present in one or both of the materials neighboring the interface, i.e., one or both of the
neighboring materials must be a dielectric. The electrical permittivity of each dielectric
material must be supplied via the *Electrical Conductivity* card (yes, this is a kludge) in
the material property file.

In its most general form, the traction condition is written

.. figure:: /figures/104_goma_physics.png
	:align: center
	:width: 90%

where **T** is the stress tensor, the superscripts ( *o* ) and ( *i* ) denote the outer and inner
phases, *n* is a unit normal pointing into the outer phase, -*H* is the local mean curvature,
and :math:`\sigma` is the surface tension.

The stress tensor can be written as the sum of the mechanical stress :math:`T_m` (e.g., the
Newtonian stress tensor) and an electrical stress 
:math:`T_e`, viz. **T** **=** :math:`T_m` **+** :math:`T_e`. The electric
stress tensor provided through this boundary condition applies to incompressible,
polarizable materials:

.. figure:: /figures/105_goma_physics.png
	:align: center
	:width: 90%

where :math:`\varepsilon` is the electrical permittivity, 
**E** **=** – :math:`\Delta` *V* is the electric field and *V* is the voltage
or electric potential.

In expanded form, the traction condition becomes

.. figure:: /figures/106_goma_physics.png
	:align: center
	:width: 90%

The *ELEC_TRACTION* boundary condition is responsible for applying either the first
or second terms on the right hand side (specified through the element block ID
parameter) whereas the *CAPILLARY* (or related boundary condition) is responsible for
the third and fourth terms.

The term multiplier supplied by the <float> input is used in the elec_surf_stress()
function (*mm_ns_bc.c*) which applies the *ELEC_TRACTION* boundary condition. It is
the etm function argument. The normal term multipliers couldn’t be used because this
boundary condition can be applied from within a material that doesn’t have the
momentum equations defined (or properly set).




.. TODO - Lines 87 and 97 have photos that needs to be replaced with the real equation.
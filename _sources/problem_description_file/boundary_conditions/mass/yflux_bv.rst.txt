************
**YFLUX_BV**
************

::

	BC = YFLUX_BV SS <bc_id> <integer1> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

The *YFLUX_BV* card enables computation of the molar flux of the specified species
using Butler-Volmer kinetics at the specified boundary (namely, the electrode surface).
When used in conjunction with the *KIN_LEAK* card, it also enables the determination
of velocity normal to the moving solid-electrode surface.

The <floatlist> consists of nine values; definitions of the input parameters are as
follows:

============== =============================================================
**YFLUX_BV**   Name of the boundary condition (<bc_name>).
**SS**         Type of boundary condition (<bc_type>), where **SS** denotes
               side set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set in
               EXODUS II) in the problem domain.
<integer1>     Species number of concentration.
<float1>       Stoichiometric coefficient.
<float2>       Kinetic rate constant.
<float3>       Reaction order.
<float4>       Anodic direction transfer coefficient.
<float5>       Cathodic direction transfer coefficient.
<float6>       Electrode potential or applied voltage.
<float7>       Theoretical open-circuit potential.
<float8>       Molecular weight of solid deposit.
<float9>       Density of solid deposit.
============== =============================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = YFLUX_BV SS 1 0 -1. 0.00001 1. 0.21 0.21 -0.8 -0.22 58.71 8.9

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




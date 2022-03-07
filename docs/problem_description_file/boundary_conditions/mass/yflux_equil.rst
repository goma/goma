***************
**YFLUX_EQUIL**
***************

::

	BC = YFLUX_EQUIL SS <bc_id> {char_string} <integer> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary card is used when equilibrium-based mass transfer is occurring at an
vapor-liquid external boundary; i.e.,

.. figure:: /figures/142_goma_physics.png
	:align: center
	:width: 90%

This is different from an internal boundary since only one phase is represented in the
computational domain. This boundary condition then describes the rate of mass
entering or leaving the boundary via vapor-liquid equilibria. The :math:`w_i^v` is the mass
fraction of component *i* in :math:`\underline{vapor}` that is in equilibrium with the liquid phase. The :math:`w_i^{v \infty}`
is the bulk concentration of component *i* in vapor.

The <float_list> requires three input values; definitions of the input parameters are as follows:

=============== =================================================================
**YFLUX_EQUIL** Name of the boundary condition.
**SS**          Type of boundary condition (<bc_type>), where **SS** denotes
                side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set in
                EXODUS II) in the problem domain.
{char_string}   This refers to the equilibrium model for mass transfer; the
                options are either **FLORY** or **RAOULT**.
<integer1>      Species id.
<float1>        Total system pressure.
<float2>        Mass transfer coefficient.
<float3>        Bulk concentration in vapor ( :math:`w_i^{v \infty}` )
=============== =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = YFLUX_EQUIL SS 1   FLORY 0   1. 5.4e-3   0.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is very similar to *VL_EQUIL* and *VL_POLY* except that it is
only applied at an external boundary where vapor phase is not modeled in the problem.



--------------
**References**
--------------

GTM-007.1: New Multicomponent Vapor-Liquid Equilibrium Capabilities in GOMA,
December 10, 1998, A. C. Sun

.. TODO - Line 18 is a photo that needs to be replaced with the correct equation.
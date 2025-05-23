************************
**CAPILLARY_SHEAR_VISC**
************************

::

	BC = CAPILLARY_SHEAR_VISC SS <bc_id> <float_list> [integer]

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition card is used to apply additional capillary forces beyond
surface tension and surface tension gradients (as applied with use of the CAPILLARY
BC) to the momentum equation on a free-surface. These additional forces are caused
by surface deformation (surface expansion/contraction/shear) in the presence of
surface-active species. Microstructural layers of surfactants in a capillary free surface
can lead to significant dissipation of mechanical energy due to an effective surface
viscosity. These additional properties are specified as inputs to this boundary
condition.

Definitions of the input parameters are as follows:

======================== ==============================================================
**CAPILLARY_SHEAR_VISC** Name of the boundary condition.
**SS**                   Type of boundary condition (<bc_type>), where **SS**
                         denotes side set in the EXODUS II database.
<bc_id>                  The boundary flag identifier, an integer associated with
                         <bc_type> that identifies the boundary location (side set in
                         EXODUS II) in the problem domain.
<float1>                 :math:`\mu_s`, surface shear viscosity.
<float2>                 :math:`\kappa_s` , surface extensional/dilatational viscosity.
[integer]                Optional integer value indicating the element block id from
                         which to apply the boundary condition. This is used to force
                         the capillary stresses to be applied from within a phase
                         where the momentum equations are defined.
======================== ==============================================================

------------
**Examples**
------------

Following is a sample card:
::

    BC = CAPILLARY SS 12   1.0 10.0 0.0

::

    BC = CAPILLARY_SHEAR_VISC SS 12   0.001 0.01

These cards specifies that capillary forces be applied to the free surface on side set 12.
If a surface tension material parameter value or model is supplied, this is the surface
tension value used. If not, the surface tension value used is 1.0. An external isotropic
pressure of 10.0 is applied from the surrounding environment. The second card adds a
surface viscosity effect. Note that you must solve the shell equation 
*EQ = n_dot_curl_v* to pick up this term.

-------------------------
**Technical Discussion**
-------------------------

The functional form of this boundary condition is

.. figure:: /figures/117_goma_physics.png
	:align: center
	:width: 90%

where *n* is the outward normal to the surface, *T* is the fluid stress tensor, 
:math:`\P_{ex}` is the
external applied pressure described above, *H* is the surface curvature defined as,
H = – :math:`\Delta_s` ⋅ n ⁄ 2, :math:`\sigma` is the surface tension, and
:math:`\Delta_s` is the surface divergence
operator defined as :math:`\Delta_s` f = (I – nn) ⋅ :math:`\Delta` f .

The Boussinesq-Scriven surface rheological constitutive equation is as follows:

.. figure:: /figures/118_goma_physics.png
	:align: center
	:width: 90%

Here, :math:`\Delta_s` = ( :math:`\underline{I}` – :math:`\underline{nn}` ) ⋅ :math:`\Delta`  is the surface gradient operator, :math:`I_s` = (:math:`\underline{I}` - (:math:`\underline{nn}` )
unit tensor. :math:`\mu_s` and :math:`\kappa_s` is the surface
are the surface shear viscosity and surface dilatational viscosity,
respectively. The terms beyond the first three on the right are added by this boundary
condition card. Note that the first three terms on the right are balance of the stress in
the standard goma CAPILLARY condition, with surface tension gradients being
accommodated through variable surface tension. The boundary condition
CAPILLARY_SHEAR_VISC is used to set the additional terms of this constitutive
equation. *As of January 2006 only the 7th term on the right hand side is implemented,
as it is the only nonzero term in a flat surface shear viscometer*. The building blocks
for the other terms are available through additional shell equations (specifically you
must solve *EQ = n_dot_curl_v* equation on the same shell surface). . These remaining
terms actually represent additional dissipation caused by surface active species
microstructures flowing in the surface. The best source of discussion of this equation
is a book by Edwards et al. (1991. *Interfacial Transport Processes and Rheology*.
Butterworth-Heinemann, Boston).




.. TODO - Lines 67 and 80 have photos that needs to be replaced with the real equation.
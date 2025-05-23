**********************
**SHARP_COX_VELOCITY**
**********************

::

	BC = SHARP_COX_VELOCITY SS <bc_id> <floatlist

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is used to induce fluid velocity at a wetting line when using
*Level Set Interface Tracking*. Its formulation is identical to the
WETTING_SPEED_COX boundary condition, but it is applied as a single point source
on the boundary instead of a distributed stress.

This boundary condition can only be used for two-dimensional problems.

A description of the input parameters follows:

====================== ===================================================================
**SHARP_COX_VELOCITY** Name of the boundary condition.
**SS**                 Type of boundary condition (<bc_type>), where **SS**
                       denotes side set in the EXODUS II database.
<bc_id>                The boundary flag identifier, an integer associated with
                       <bc_type> that identifies the boundary location (side
                       set in EXODUS II) in the problem domain.
<float1>               :math:`\theta_s`, the static contact angle in degress.
<float2>               :math:`\sigma` is the surface tension.
<float3>               :math:`\varepsilon_s` is the dimensionless slip length,
                       i.e. the ratio of the slip length to the characteristic 
                       length scale of the macroscopic flow.
<float4>               :math:`\beta`, slip coefficient.
<float5>               :math:`t_{relax}` is a relaxation time which can be used
                       to smooth the imposed contact point velocity for transient problems. Set to zero for no smoothing.
<float6>               :math:`V_{old}` is an initial velocity used in velocity smoothing 
                       for transient problems. Set to zero when smoothing is not used.
====================== ===================================================================

------------
**Examples**
------------

An example:
::

   BC = SHARP_COX_VELOCITY SS 10 30.0 72.0 0.01 0.1   0   0

-------------------------
**Technical Discussion**
-------------------------

The implementation for this wetting condition is identical to that of
SHARP_WETLIN_VELOCITY, but the wetting velocity dependence is different.
Because the wetting stress is not applied at a point, it is most appropriate for use when
using subelement integration which similarly collapses the surface tension sources
associated with the interface onto the interfacial curve.

Note also that this boundary condition is strictly for use with two-dimensional
problems. Attempting to apply it to a three dimensional problem will result in an error
message.

----------
**Theory**
----------

Derivation of the force condition for this boundary condition starts with a relation for
wetting line velocity

.. figure:: /figures/211_goma_physics.png
	:align: center
	:width: 90%

where :math:`V_{Cox}` is computed from the Cox hydrodynamic wetting theory;

.. figure:: /figures/212_goma_physics.png
	:align: center
	:width: 90%

See VELO_THETA_COX for details of the Cox functions f and g. Note that the
parameters :math:`\lambda`, :math:`q_{inner}`, and :math:`q_{outer}` are currently not accessible from the input card and are
hard-set to zero. :math:`\lambda` is the ratio of gas viscosity to liquid viscosity whereas :math:`q_{inner}` and
:math:`q_{outer}` represent influences from the inner and outer flow regions.

Note that the convention for contact angles in this relation is that values of θ near to
zero indicate a high degree of wetting and values of :math:`\theta` near 180 ° indicate the opposite.
This is mapped to a stress value by analogy with Navier’s slip relation and has the
following form when the velocity smoothing is not used, 

.. figure:: /figures/213_goma_physics.png
	:align: center
	:width: 90%

The Cox wetting velocity requires evaluation of integrals for the function 
g(:math:`\theta`, :math:`\lambda`) which
is currently done numerically using 10-point Gaussian quadrature. As such the
evaluation of the integrals is expected to become inaccurate as either 
:math:`\theta_s` tends toward
zero or :math:`\theta` tends toward 180 degrees. Note that the integrand becomes singular as :math:`\theta` tends toward 0 or 180 degrees.


--------------
**References**
--------------

No References. 

.. TODO -Lines 73, 79 and 93 have pictures that need to be swapped with the correct equations.

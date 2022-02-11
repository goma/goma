************************
**SHARP_BLAKE_VELOCITY**
************************

::

	BC = SHARP_BLAKE_VELOCITY SS <bc_id> <floatlist

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is used to induce fluid velocity at a wetting line when using
*Level Set Interface Tracking*. Its formulation is identical to the
WETTING_SPEED_BLAKE boundary condition, but it is applied as a single point
source on the boundary instead of a distributed stress.

This boundary condition can only be used for two-dimensional problems.

A description of the input parameters follows:

======================== ======================================================
**SHARP_BLAKE_VELOCITY** Name of the boundary condition.
**SS**                   Type of boundary condition (<bc_type>), where **SS**
                         denotes side set in the EXODUS II database.
<bc_id>                  The boundary flag identifier, an integer associated with
                         <bc_type> that identifies the boundary location 
                         (side set in EXODUS II) in the problem domain.
<float1>                 :math:`\theta_s`, the static contact angle in degress.
<float2>                 V_0 is a pre-exponential velocity factor (see functional
                         form below).
<float3>                 g is a thermally scaled surface tension, i.e. :math:`\sigma`
                         /2nkT.
<float4>                 :math:`\beta`, slip coefficient.
<float5>                 :math:`t_{relax}` is a relaxation time which can be used to  
                         smooth the imposed contact point velocity for transient problems. Set to zero for no smoothing.
<float6>                 :math:`V_{old}` is an initial velocity used in velocity 
                         smoothing for transient problems. Set to zero when 
                         smoothing is not used.
======================== ======================================================

------------
**Examples**
------------

An example:
::

   BC = SHARP_BLAKE_VELOCITY SS 10 30.0 0.1 8. 0.001 0   0

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

Derivation of the force condition for this boundary condition starts with a simple
relation for wetting line velocity

.. figure:: /figures/209_goma_physics.png
	:align: center
	:width: 90%

Note that the convention for contact angles in this relation is that values of θ near to
zero indicate a high degree of wetting and values of θ near 180 ° indicate the opposite.
This is mapped to a stress value by analogy with Navier’s slip relation and has the
following form when the velocity smoothing is not used,

.. figure:: /figures/210_goma_physics.png
	:align: center
	:width: 90%


--------------
**References**
--------------

No References. 

.. TODO -Lines 73 and 82 have pictures that need to be swapped with the correct equations.
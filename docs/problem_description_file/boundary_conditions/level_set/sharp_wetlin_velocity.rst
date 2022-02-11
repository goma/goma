*************************
**SHARP_WETLIN_VELOCITY**
*************************

::

	BC = SHARP_WETLIN_VELOCITY SS <bc_id> <floatlist

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is used to induce fluid velocity at a wetting line when using
*Level Set Interface Tracking*. Its formulation is identical to the
WETTING_SPEED_LINEAR boundary condition, but it is applied as a single point
source on the boundary instead of a distributed stress.

This boundary condition can only be used for two-dimensional problems.

A description of the input parameters follows:

========================= ======================================================
**SHARP_WETLIN_VELOCITY** Name of the boundary condition.
**SS**                    Type of boundary condition (<bc_type>), where **SS**
                          denotes side set in the EXODUS II database.
<bc_id>                   The boundary flag identifier, an integer associated with
                          <bc_type> that identifies the boundary location (side 
                          set in EXODUS II) in the problem domain.
<float1>                  :math:`\theta_s`, the static contact angle in degress.
<float2>                  :math:`c_T`, proportionality constant as defined below
<float3>                  Currently not used.
<float4>                  :math:`\beta`, slip coefficient.
========================= ======================================================

------------
**Examples**
------------

An example:
::

   BC = SHARP_WETLIN_VELOCITY SS 10 30.0 0.1 0. 0.001

-------------------------
**Technical Discussion**
-------------------------

As noted above, this boundary condition imposes the same wetting stress dependence
as the WETTING_SPEED_LINEAR boundary condition. However, its application in
the FEM context is different. Instead of the wetting stress, 
:math:`\tau_w`, being applied according to the formula:

.. figure:: /figures/217_goma_physics.png
	:align: center
	:width: 90%

as is the case for the WETTING_SPEED_LINEAR condition, the Dirac function is
used to remove the integral and replace it with a point stress at the location where
:math:`\phi` = 0 on the boundary. Designating this point as
:math:`X_{cl}`, the vector applied to the momentum equation is given by

.. figure:: /figures/218_goma_physics.png
	:align: center
	:width: 90%

Because the wetting stress is not applied at a point, it is most appropriate for use when
using subelement integration which similarly collapses the surface tension sources
associated with the interface onto the interfacial curve. Note that this method of
application is identical to the SHARP_CA_2D boundary condition discussed
elsewhere.

Note also that this boundary condition is strictly for use with two-dimensional
problems. Attempting to apply it to a three dimensional problem will result in an error
message.

----------
**Theory**
----------

Derivation of the force condition for this boundary condition starts with a simple
relation for wetting line velocity

.. figure:: /figures/219_goma_physics.png
	:align: center
	:width: 90%

Note that the convention for contact angles in this relation is that values of 
:math:`\theta` near to
zero indicate a high degree of wetting and values of :math:`\theta` near 180 ° indicate the opposite.
This is mapped to a stress value by analogy with Navier’s slip relation,

.. figure:: /figures/220_goma_physics.png
	:align: center
	:width: 90%

It should be noted that there is no distinction for this model in the function of 
:math:`\beta` or :math:`c_T`.
The two parameters are interchangeable. In non-linear models, (see
WETTING_SPEED_BLAKE) this is no longer true.


--------------
**References**
--------------

No References. 

.. TODO -Lines 55, 64, 85 and 94 have pictures that need to be swapped with the correct equations.

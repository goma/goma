*******************
**BLAKE_DIRICHLET**
*******************

::

	BC = BLAKE_DIRICHLET SS <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(SIC/VECTOR MOMENTUM)**

This boundary condition is used to induce fluid velocity at a wetting line when using
*Level Set Interface Tracking*. It is an alternative to the WETTING_SPEED_BLAKE
BC which does not require the VELO_SLIP_LS or VELO_SLIP_FILL BC's on the
wetting boundary. It uses a Blake-DeConninck relationship between apparent contact
angle and wetting velocity. As the name implies, this boundary condition differs from
WETTING_SPEED_BLAKE in that wetting velocity is in a strong fashion on the
wetting boundary.

A description of the input parameters follows:

=================== ============================================================
**BLAKE_DIRICHLET** Name of the boundary condition.
**SS**              Type of boundary condition (<bc_type>), where **SS** denotes
                    side set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (side set in
                    EXODUS II) in the problem domain.
<float1>            :math:`\theta_s`, the static contact angle in degrees.
<float2>            :math:`V_0` is a pre-exponential velocity factor 
                    (see functional form below). (L/T)
<float3>            *g* is a thermally scaled surface tension, i.e. 
                    :math:`\sigma`/2nkT. Note that this parameter will be multiplied by the surface tension supplied in the material file when its used in the wetting velocity relation.
<float4>            *w*, is the width of the interface wetting region. It defaults to
                    the level set length scale if zero of less.
<float5>            :math:`\tau`, stability parameter (T).
<float6>            :math:`v_{sx}`, x-component of substrate velocity.
<float7>            :math:`v_{sy}`, y-component of substrate velocity (L/T).
<float8>            :math:`v_{sz}`, z-component of substrate velocity.
=================== ============================================================

------------
**Examples**
------------

Here is an example card:
::

   BC = BLAKE_DIRICHLET SS 10 30.0 20.1 7.0 0.0 0.001 0. 0. 0.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is an additional means to impose a wetting line velocity at the
contact line for level set interface tracking problems. It is related to the
WETTING_SPEED_BLAKE condition in that it uses the same Blake-DeConninck
relationship between contact angle and wetting speed, but it applies this relation to the
computational setting in a different way.

In this case, the following vector constraint is added to fluid momentum equation on
the sideset to which this boundary condition is applied:

.. figure:: /figures/237_goma_physics.png
	:align: center
	:width: 90%

The factor *P* is a large penalty parameter which swamps any contributions from the
volumetric momentum equations. Thus, the velocity, :math:`\sec{v}` , on this boundary will be set
solely by the preceding constraint. In this sense, it is a Dirichlet condition (strictly
speaking, Dirichlet conditions involve direct substitution of nodal degrees-of-freedom
with corresponding elimination of its equation from the matrix which this boundary
condition does NOT do).

In the preceding, the vector :math:`\sec{t}`, is a tangent vector to the surface and always points in the
same direction of the level set gradient on the boundary (that is, from negative to
positive). In three dimensions, :math:`\sec{t}` will also be normal to the contact line curve as it
intersects the surface itself.

The masking function, f(:math:`\phi`;w) , is used to limit the application of the wetting line
velocity to only that region of the boundary that is the in immediate vicinity of the
contact line. We use a simple “hat” function:

.. figure:: /figures/238_goma_physics.png
	:align: center
	:width: 90%

Needless to say, f(:math:`\phi`;w) is identically zero for level set values outside the 
interval (-*w*,*w*).

The stabilization term, -:math:`\tau` :math:`\frac{∂v}{∂t}` , is intended to introduce something like inertia to the
wetting line. That is to say, it’s primary effect is to limit the rate of change of the
wetting line velocity to “reasonable” values. The :math:`\tau` parameter should be chosen to be
on the order of the smallest anticipated time step size in the problem. Setting it at zero,
of course, will remove this term entirely.

In general, this boundary condition can be used to exclusively to set both the wetting
speed velocity and the no slip requirement on the indicated sideset set. This would also
include the no penetration requirement. The user may, however, find it advantageous
to apply this constraint directly with the VELO_NORMAL condition on the same side
set.

An additional note is that the “scaled viscosity” parameter g will be multiplied by the
surface tension value supplied with Surface Tension card in the material file.

----------
**Theory**
----------

The wetting speed model for this boundary condition is the same used by the
WETTING_SPEED_BLAKE card:


.. figure:: /figures/239_goma_physics.png
	:align: center
	:width: 90%


--------------
**References**
--------------

T. D. Blake and J. De Coninck 2002. “The Influence of Solid-Liquid Interactions on
Dynamic Wetting”, Advances in Colloid and Interface Science, 96, 21-36.

.. TODO -Lines 67, 87 and 117 have pictures that need to be swapped with the correct equations.

*****************
**COX_DIRICHLET**
*****************

::

	BC = COX_DIRICHLET SS <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(SIC/VECTOR MOMENTUM)**

This boundary condition is used to induce fluid velocity at a wetting line when using
*Level Set Interface Tracking*. It is an alternative to the WETTING_SPEED_COX
which does not require the VELO_SLIP_LS or VELO_SLIP_FILL BC's on the wetting
boundary. It implements a version of the Cox hydrodynamic model of wetting (see
below). As the name implies, this boundary condition differs from
WETTING_SPEED_COX in that wetting velocity is applied in a strong fashion on the
wetting boundary.

A description of the input parameters follows:

================= ===========================================================
**COX_DIRICHLET** Name of the boundary condition.
**SS**            Type of boundary condition (<bc_type>), where **SS**
                  denotes side set in the EXODUS II database.
<bc_id>           The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side 
                  set in EXODUS II) in the problem domain.
<float1>          :math:`\theta_s`, the static contact angle in degress.
<float2>          :math:`\varepsilon_s` is the dimensionless slip length, 
                  i.e. the ratio of the slip
                  length to the characteristic length scale of the 
                  macroscopic flow. (L)
<float3>          :math:`\sigma` is the surface tension. Note this value 
                  will be scaled by the surface tension value supplied in
                  the material file.(F/L)
<float4>          *w*, is the width of the interface wetting region. It 
                  defaults to the level set length scale if zero of less (L).
<float5>          :math:`\tau`, stability parameter (T).
<float6>          :math:`v_{sx}`, x-component of substrate velocity.
<float7>          :math:`v_{sy}`, y-component of substrate velocity (L/T).
<float8>          :math:`v_{sz}`, z-component of substrate velocity.
================= ===========================================================

------------
**Examples**
------------

Here is an example card:
::

   BC = COX_DIRICHLET SS 10 30.0 0.01 72.0   0. 0.001 0. 0. 0.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is applied in exactly the same manner as the
BLAKE_DIRICHLET boundary condition. The only substantial difference is the
model used to derive the wetting speeds relation to the local apparent contact angle.
The reader is referred to the BLAKE_DIRICHLET section of the manual for further
reference.

----------
**Theory**
----------

This boundary condition uses this relation that represents the Cox hydrodynamic
wetting model

.. figure:: /figures/240_goma_physics.png
	:align: center
	:width: 90%

See VELO_THETA_COX for details of the Cox functions f and g. Note that the
parameters :math:`\lambda`, :math:`q_{inner}`, and :math:`q_{outer}` are currently not accessible from the input card and are
hard-set to zero. :math:`\lambda` is the ratio of gas viscosity to liquid viscosity whereas :math:`q_{inner}` and
:math:`q_{outer}` represent influences from the inner and outer flow regions.


--------------
**References**
--------------

Stephan F. Kistler 1993. “Hydrodynamics of Wetting” in Wettability, edited by John
Berg, Surfactant Science Series, 49, Marcel Dekker, NewYork, NY, pp. 311-429.

.. TODO -Line 74 has am image that needs to be switched with the equation.
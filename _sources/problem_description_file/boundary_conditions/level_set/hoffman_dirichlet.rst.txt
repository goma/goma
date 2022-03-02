*********************
**HOFFMAN_DIRICHLET**
*********************

::

	BC = HOFFMAN_DIRICHLET SS <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(SIC/VECTOR MOMENTUM)**

This boundary condition is used to induce fluid velocity at a wetting line when using
*Level Set Interface Tracking*. It is an alternative to the
WETTING_SPEED_HOFFMAN boundary condition which does not require the
VELO_SLIP_LS or VELO_SLIP_FILL BC's on the wetting boundary. As the name
implies, this boundary condition differs from WETTING_SPEED_HOFFMAN in that
wetting velocity is in a strong fashion on the wetting boundary.

A description of the input parameters follows:

===================== =======================================================
**HOFFMAN_DIRICHLET** Name of the boundary condition.
**SS**                Type of boundary condition (<bc_type>), where **SS**
                      denotes side set in the EXODUS II database.
<bc_id>               The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location 
                      (side set in EXODUS II) in the problem domain.
<float1>              :math:`\theta_s`, the static contact angle in degress.
<float2>              not used.
<float3>              :math:`\sigma`, the surface tension. Note that this 
                      parameter will be scaled by the surface tension value 
                      supplied in the material file. (F/L)
<float4>              *w*, is the width of the interface wetting region. It 
                      defaults to the level set length scale if zero of less.
<float5>              :math:`\tau`, stability parameter (T).
<float6>              :math:`v_{sx}`, x-component of substrate velocity.
<float7>              :math:`v_{sy}`, y-component of substrate velocity (L/T).
<float8>              :math:`v_{sz}`, z-component of substrate velocity.
===================== =======================================================

------------
**Examples**
------------

Here is an example card:
::

   BC = HOFFMAN_DIRICHLET SS 10 30.0 0 72.0 0.0 0.001 0. 0. 0.

-------------------------
**Technical Discussion**
-------------------------

The technical details of the application of this boundary differ not at all from those
described for the BLAKE_DIRICHLET boundary condition. The user is referred to
that section for further details. This boundary condition differs only in the model used
to determine the wetting velocity. This is described below and in the
VELO_THETA_HOFFMAN card.

----------
**Theory**
----------

Derivation of this boundary condition starts with a relation that represents the Hoffman
wetting correlation

.. figure:: /figures/241_goma_physics.png
	:align: center
	:width: 90%

See VELO_THETA_HOFFMAN for details of the Hoffman function g. Note that the
convention for contact angles in this relation is that values of :math:`\theta` near to zero indicate a
high degree of wetting and values of :math:`\theta` near 180 ° indicate the opposite. This is mapped
to a stress value by analogy with Navier’s slip relation,

.. figure:: /figures/242_goma_physics.png
	:align: center
	:width: 90%


--------------
**References**
--------------

Stephan F. Kistler 1993. “Hydrodynamics of Wetting” in Wettability, edited by John
Berg, Surfactant Science Series, 49, Marcel Dekker, NewYork, NY, pp. 311-429.

.. TODO -Lines 70 and 79 have pictures that need to be swapped with the correct equations.

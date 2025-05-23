***********************
**WETTING_SPEED_BLAKE**
***********************

::

	BC = WETTING_SPEED_BLAKE SS <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is used to induce fluid velocity at a wetting line when using
*Level Set Interface Tracking*. It implements a version of the Blake-DeConinck
molecular-kinetic theory wetting model.

A description of the input parameters follows:

======================= ==========================================================
**WETTING_SPEED_BLAKE** Name of the boundary condition.
**SS**                  Type of boundary condition (<bc_type>), where **SS**
                        denotes side set in the EXODUS II database.
<bc_id>                 The boundary flag identifier, an integer associated with
                        <bc_type> that identifies the boundary location (side set in
                        EXODUS II) in the problem domain.
<float1>                :math:`\theta_s`, the static contact angle, degrees .
<float2>                V_0 is a pre-exponential velocity factor (see functional
                        form below).
<float3>                g is a thermally scaled surface tension, i.e. :math:`\sigma`/2nkT.
<float4>                *w*, width of interfacial region near contact line. Defaults to
                        level set length scale if zero or less.
<float5>                :math:`\beta`, slip coefficient.
<float6>                currently not used.
<float7>                currently not used.
<float8>                currently not used.
======================= ==========================================================

------------
**Examples**
------------

An example:
::

   BC = WETTING_SPEED_BLAKE SS 10 30.0 20.1 7.0 0.   0.001 0. 0. 0.

-------------------------
**Technical Discussion**
-------------------------

The implementation for this wetting condition is identical to that of
WETTING_SPEED_LINEAR, but the wetting velocity dependence is different.

Note that it is a requirement that when using this boundary condition that slip to some
extent be allowed on this boundary. This is most often done by applying a
VELO_SLIP_LS boundary condition in conjunction with this boundary condition. In
addition, a no penetration condition on the velocity is need in either the form of a
Dirichlet condition or a VELO_NORMAL condition. It is important to note that the
slipping condition need not relax the no slip requirement completely. In fact, its
parameters should be set so that no slip is for the most part satisfied on the boundary in
regions away from the contact line. Near the contact line however the parameters in
the slip condition and the WETTING_SPEED_BLAKE condition need to be fixed so
that appreciable fluid velocity is induced. This is a trial and error process at the current
time.

----------
**Theory**
----------

Derivation of this boundary condition starts with a relation propose by Blake and
DeConinck for wetting line motion

.. figure:: /figures/221_goma_physics.png
	:align: center
	:width: 90%

This is mapped to a stress value by analogy with Navier’s slip relation,

.. figure:: /figures/222_goma_physics.png
	:align: center
	:width: 90%

This relation contrasts with the “linear” relation applied by the
WETTING_SPEED_LINEAR relation in that the rate of change of the wetting velocity
with the contact angle decreases as the wetting angle deviates more and more from its
static value. This is more consisten with physical behaviors that the linear model.

In point of fact this condition is a vector condition so this scalar stress value multiplies
the unit vector tangent to the surface and normal to the contact line, 
:math:`\vec{t}` . This stress is
then weighted by smooth Dirac function to restrict its location to being near the
interface, weighted by a FEM shape function, integrated over the boundary sideset and
added to the fluid momentum equation for the corresponding node j, vis:

.. figure:: /figures/223_goma_physics.png
	:align: center
	:width: 90%


--------------
**References**
--------------

No References. 

.. TODO -Lines 75, 81 and 97 have pictures that need to be swapped with the correct equations.

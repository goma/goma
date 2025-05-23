************************
**WETTING_SPEED_LINEAR**
************************

::

	BC = WETTING_SPEED_LINEAR SS <bc_id> <floatlist

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is used to induce fluid velocity at a wetting line when using
*Level Set Interface Tracking*.

A description of the input parameters follows:

======================== =======================================================
**WETTING_SPEED_LINEAR** Name of the boundary condition.
**SS**                   Type of boundary condition (<bc_type>), where **SS**
                         denotes side set in the EXODUS II database.
<bc_id>                  The boundary flag identifier, an integer associated with
                         <bc_type> that identifies the boundary location 
                         (side set in EXODUS II) in the problem domain.
<float1>                 :math:`\theta_s`, the static contact angle in degrees.
<float2>                 :math:`c_T`, proportionality constant as defined below
<float3>                 *w* ,width of interfacial region near contact line. 
                         Defaults to level set length scale if zero or less (L).
<float4>                 :math:`\beta`, slip coefficient.
<float5>                 currently not used.
<float6>                 currently not used.
<float7>                 currently not used.
======================== =======================================================

------------
**Examples**
------------

An example:
::

   BC = WETTING_SPEED_LINEAR SS 10 30.0 0.1 0. 0.001 0. 0. 0.

-------------------------
**Technical Discussion**
-------------------------

The prescence of wetting or contact lines in problems using level set interface tracking
introduces the problem of modeling the motion of the wetting line. This boundary
condition presents one potential means for doing this. It adds a wall stress value at the
boundary in a region near to the wetting line (this is set with the Level Set Length Scale
discussed previously). This wall stress value depends upon the deviation of the
apparent contact angle determined from the level set function and a set static contact
angle. The bigger the deviation in principle the bigger the induced stress. The stress is
modeled by analogy with Navier’s slip relation (with slip coefficient :math:`\beta`).
The stress will induce a fluid velocity at the boundary which it is hoped will move the contact line at a velocity that is consistent with the rest of the flow.

An important note is that it is a requirement that when using this boundary condition
that slip to some extent be allowed on this boundary. This is most often done by
applying a VELO_SLIP_LS boundary condition in conjunction with this boundary
condition. In addition, a no penetration condition on the velocity is need in either the
form of a Dirichlet condition or a VELO_NORMAL condition. It is important to note
that the slipping condition need not relax the no slip requirement completely. In fact,
its parameters should be set so that no slip is for the most part satisfied on the boundary
in regions away from the contact line. Near the contact line however the parameters in
the slip condition and the WETTING_SPEED_LINEAR condition need to be fixed so
that appreciable fluid velocity is induced. This is a trial and error process at the current
time.

----------
**Theory**
----------

Derivation of this boundary condition starts with a simple relation for wetting line
velocity

.. figure:: /figures/227_goma_physics.png
	:align: center
	:width: 90%

Note that the convention for contact angles in this relation is that values of 
:math:`\theta` near to
zero indicate a high degree of wetting and values of :math:`\theta` near 180 ° indicate the opposite.
This is mapped to a stress value by analogy with Navier’s slip relation

.. figure:: /figures/228_goma_physics.png
	:align: center
	:width: 90%

It should be noted that there is no distinction for this model in the function of 
:math:`\beta` or :math:`c_T`.
The two parameters are interchangeable. In non-linear models, (see
WETTING_SPEED_BLAKE) this is no longer true.

In point of fact this condition is a vector condition so this scalar stress value multiplies
the unit vector tangent to the surface and normal to the contact line,
:math:`\vec{t}` . This stress is
then weighted by smooth Dirac function to restrict its location to being near the
interface, weighted by a FEM shape function, integrated over the boundary sideset and
added to the fluid momentum equation for the corresponding node j, vis:

.. figure:: /figures/229_goma_physics.png
	:align: center
	:width: 90%


--------------
**References**
--------------

No References. 

.. TODO -Lines 79, 88 and 104 have pictures that need to be swapped with the correct equations.

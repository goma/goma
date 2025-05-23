**********************
**LINEAR_WETTING_SIC**
**********************

::

	BC = LINEAR_WETTING_SIC SS <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(SIC/VECTOR MOMENTUM)**

This boundary condition is used to induce fluid velocity at a wetting line when using
*Level Set Interface Tracking*. It is an alternative to the WETTING_SPEED_LINEAR
BC which does not require the VELO_SLIP_LS or VELO_SLIP_FILL BC's on the
wetting boundary.

A description of the input parameters follows:

====================== =======================================================
**LINEAR_WETTING_SIC** Name of the boundary condition.
**SS**                 Type of boundary condition (<bc_type>), where **SS** 
                       denotes side set in the EXODUS II database.
<bc_id>                The boundary flag identifier, an integer associated with
                       <bc_type> that identifies the boundary location 
                       (side set in EXODUS II) in the problem domain.
<float1>               :math:`\theta_s`, the static contact angle in degress.
<float2>               :math:`c_T`, proportionality constant as defined below.
<float3>               *w*, width of interfacial region near contact line. Defaults to
                       level set length scale if zero or less (L).
<float4>               :math:`\beta`, slip coefficient.
<float5>               :math:`v_{sx}`, x-component of substrate velocity.
<float6>               :math:`v_{sy}`, y-component of substrate velocity.
<float7>               :math:`v_{sz}`, z-component of substrate velocity.
<float8>               :math:`\tau`, stability parameter.
====================== =======================================================

------------
**Examples**
------------

Here is an example card:
::

   BC = LINEAR_WETTING_SIC SS 10 30.0 0.1 0. 0.001 0. 0. 0. 0.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is an additional means to impose a wetting line velocity at the
contact line for level set interface tracking problems. The boundary condition uses a
form of the Navier-Stokes slip condition to impose a boundary shear stress term to the
momentum equation:

.. figure:: /figures/233_goma_physics.png
	:align: center
	:width: 90%

where :math:`\vec{n}` and :math:`\vec{t}` are the normal and tangent boundary vectors, respectively, :math:`\beta` is the
“slipping” parameter which in this context is used actually as a penalty parameter,
:math:`\vec{v}_s` is the substrate velocity, :math:`\tau` is a stabilization parameter, 
:math:`V_{wet}` is the wetting velocity given by the following relation

.. figure:: /figures/234_goma_physics.png
	:align: center
	:width: 90%

The masking function f(F) is given by the following relation as well:

.. figure:: /figures/235_goma_physics.png
	:align: center
	:width: 90%

where :math:`\alpha` is the width of the interfacial region near the contact line itself. It has the
effect of “turning off” the wetting velocity at points on the boundary away from the
interface.

This constraint is then introduced into the fluid momentum equation via the weak
natural boundary condition term:

.. figure:: /figures/236_goma_physics.png
	:align: center
	:width: 90%

When applying this boundary condition, the user should choose a value for 
:math:`\beta` which is
relatively small. Its size is dictated by the requirement that away from the interface this
boundary condition should be imposing a no-slip condition on the fluid velocity.
Conversely, in the vicinity of the wetting line this boundary condition will impose the
wetting velocity as computed from the preceding equation.

This boundary condition probably should be used in conjunction with a no penetration
boundary condition, for example, a VELO_NORMAL condition on the same sideset or
potentially a Dirichlet condition on velocity if the geometry permits this. In theory, this
boundary condition can be used to impose no penetration as well, but this will require a
very small value for :math:`\beta`. The user should experiment with this.

The stability parameter, :math:`\tau`, as requires commentary. It is helpful to imagine that this
parameter introduces a certain amount of inertia to motion of the contact line. With this
term active (non-zero value for :math:`\tau`), large changes of the contact line velocity with time
are restricted. This can be quite helpful during startup when the intial contact angle is
often very different from its equilibrium value and there can be very large velocities
generated as a result. These may in turn lead to low time step size and other numerical
problems.

Although every situation is different, one should choose values for :math:`\tau` which are on the
order of 1 to 10 times the starting time step size of the simulation. One should also
recognize that this term is not consistent from a physical standpoint and therefore one
should endeavor to keep :math:`\tau` as small as possible if not in fact equal to zero.



--------------
**References**
--------------

No References. 

.. TODO -Lines 58, 67, 73 and 84 have pictures that need to be swapped with the correct equations.

******************
**VELO_THETA_TPL**
******************

::

	BC = VELO_THETA_TPL NS <bc_id> <float_list> [integer]

-----------------------
**Description / Usage**
-----------------------

**(PCC/R_MOM_TANG1)**

This boundary condition card applies a dynamic contact angle dependent velocity
condition in place of one component of the fluid momentum equation (unlike CA_BC
which applies a fixed contact angle condition ot the mesh equation). The functional
form of this condition is given below and stems from the Blake-DeConinck psuedomolecular
kinetics theory of wetting. It is recommended that this condition or other
forms of it (cf. VELO_THETA_HOFFMAN or VELO_THETA_COX) be used for
steady and transient ALE problems. If you are deploying level-set technology to track
moving capillary surfaces and three-phase wetting lines then the counterpart to this
condition is WETTING_SPEED_LINEAR. It is noteworthy that this condition is
applied to the fluid momentum equation, so that the velocity of the wetting line and the
cosine of the current measured contact angle difference with the specified static value
are related in a linear way.

The <float_list> for this boundary condition has eight values; definitions of the input
parameters are as follows:

================== =======================================================================
**VELO_THETA_TPL** Name of the boundary condition.
**NS**             Type of boundary condition (<bc_type>), where **NS** denotes
                   node set in the EXODUS II database.
<bc_id>            The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (node set in
                   EXODUS II) in the problem domain.
<float1>           :math:`\theta`, equilibrium (static) contact angle subtended by wall
                   normal and free surface normal, in units of degrees.
<float2>           :math:`n_x` , x-component of normal vector to the geometry
                   boundary (see important note below regarding variable wall
                   normals, viz. non-planar solid walls).
<float3>           :math:`n_y` , y-component of normal vector to the geometry
                   boundary. (see important note below regarding variable wall
                   normals, viz. non-planar solid walls).
<float4>           :math:`n_z` , z-component of normal vector to the geometry
                   boundary. (see important note below regarding variable wall
                   normals, viz. non-planar solid walls).
<float5>           V_0 is a pre-exponential velocity factor (see functional
                   form below)
<float6>           *g* is a thermally scaled surface tension, i.e. 
                   :math:`\sigma` /2nkT. This value is multiplied by the surface tension value stipulated by the surface tension material model.
<float7>           :math:`t_{relax}` is a relaxation time which can be used to smooth the
                   imposed contact point velocity for transient problems. Set
                   to zero for no smoothing.
<float8>           :math:`V_{old}` is an initial velocity used in velocity smoothing for
                   transient problems. Set to zero when smoothing is not used.
[integer]          *blk_id*, an optional parameter that is the element block
                   number in conjugate problems that identifies the material
                   region to which the contact angle applies (usually the liquid
                   element block in solid/liquid conjugate problems). NOTE/
                   WARNING: As of 1/13/2013 this option seems not to work
                   with TALE problems.
================== =======================================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = VELO_THETA_TPL NS   100   {45}   0.   1.  0.  1000.0   5.e-4   0   0

This condition applies a contact angle of 45 degrees between the free surface normal at
the 100 nodeset and the vector (0,1,0). The velocity scale is 1000 and the sensitivity
scale is 5.e-4. Normally, this latter vector is the normal to the solid surface in contact
with the free surface at this point.

-------------------------
**Technical Discussion**
-------------------------

* The constraint that is imposed at the nodeset node is as follows:

.. figure:: /figures/128_goma_physics.png
	:align: center
	:width: 90%

.. figure:: /figures/129_goma_physics.png
	:align: center
	:width: 90%

When smoothing is not used, i.e. :math:`t_{relax}` = 0 , the imposed velocity is equal to that
stipulated by the Blake-DeConinck equation. Also see
WETTING_SPEED_LINEAR and WETTING_SPEED_BLAKE for level-set
versions of this and VELO_THETA_HOFFMAN and VELO_THETA_COX for
other functional forms.

* We recommend use of this condition over CA_BC for all transient problems. In
  this case this condition displaces a momentum equation component, with the other
  component being used to enforce no substrate penetration. The kinematic
  condition is applied to the mesh motion a this node so as to conserve mass.

* For steady problems, the substrate velocity will be extracted from adjoining
  VELO_TANGENT, VELO_SLIP, or VELO_SLIP_ROT boundary conditions.

* **Important: Variable Wall Normals**. Situations for which the wall shape is nonplanar,
  meaning that the normal vector is not invariant as the contact line moves,
  there is an option to leave all of the normal-vector components zero. In this case
  *Goma* then seeks to determine the local wall normal vector from the geometry it is
  currently on, using the element facets. It is recommended that this option not be
  used unless the geometry is truly nonplanar, as the logic is complex and not 100%
  reliable. See documentation for CA_BC for an example.

* This condition was motivated by T. D. Blake and is the so-called Blake-
  DeConinck condition (T. D. Blake, J. De Coninck 2002. “The influence of solidliquid
  interactions on dynamic wetting”, Advances in Colloid and Interface
  Science 96, 21-36. ). See this article for some options for the form of the preexponential velocity, V_0.

* **Important: Wall Normal convention**. The wall normal vector on an external
  solid boundary is defined in goma as the inward facing normal to the mesh, and the
  free surface normal to the liquid (or wetting phase for two-liquid systems) is
  defined as the outward facing normal to the free surface. Put another way and
  referring to the picture below, the wall normal is directed from the “solid phase” to
  the “liquid phase”, and the free surface normal is directed from the “liquid phase”
  or “wetting phase” to the “vapor phase” or “Non-wetting phase”. Note that for
  zero contact angle the liquid is “perfectly wetting”. The air-entrainment limit (viz.
  the hydrodynamic theory interpretation) would occure at a 180 degree contact
  angle. Recall that the angle is specified in radians on this card.

.. figure:: /figures/130_goma_physics.png
	:align: center
	:width: 90%




.. TODO - Lines 86, 90, and 132 have photos that needs to be replaced with the real equation.

**********************
**VELO_THETA_HOFFMAN**
**********************

::

	BC = VELO_THETA_HOFFMAN NS <bc_id> <float_list> [integer]

-----------------------
**Description / Usage**
-----------------------

**(PCC/R_MOM_TANG1)**

This boundary condition card applies a dynamic contact angle dependent velocity
condition in a similiar fashion to VELO_THETA_TPL, but the functional form of the
velocity is different. The functional form stems not from a theory of wetting, but
instead, from a correlation of many empirical measurements.

The <float_list> for this boundary condition has eight values; definitions of the input
parameters are as follows:

====================== ==============================================================
**VELO_THETA_HOFFMAN** Name of the boundary condition.
**NS**                 Type of boundary condition (<bc_type>), where **NS**
                       denotes node set in the EXODUS II database.
<bc_id>                The boundary flag identifier, an integer associated with
                       <bc_type> that identifies the boundary location 
                       (node set in EXODUS II) in the problem domain.
<float1>               :math:`\theta`, equilibrium (static) contact angle 
                       subtended by wall normal and free surface normal,
                       in units of degrees.
<float2>               :math:`n_x` , x-component of normal vector to the geometry
                       boundary (see important note below regarding variable wall
                       normals, viz. non-planar solid walls).
<float3>               :math:`n_y` , y-component of normal vector to the geometry 
                       boundary. (see important note below regarding variable wall
                       normals, viz. non-planar solid walls).
<float4>               :math:`n_z` , z-component of normal vector to the geometry
                       boundary. (see important note below regarding variable wall
                       normals, viz. non-planar solid walls).
<float5>               currently not used.
<float6>               :math:`\sigma` is the surface tension. This value is  
                       multiplied by the surface tension value stipulated by 
                       the surface tension material model.
<float7>               :math:`t_{relax}` is a relaxation time which can be 
                       used to smooth the imposed contact point velocity for 
                       transient problems. Set to zero for no smoothing.
<float8>               :math:`V_{old}` is an initial velocity used in velocity
                       smoothing for transient problems. Set to zero when 
                       smoothing is not used.
[integer]              *blk_id*, an optional parameter that is the element block
                       number in conjugate problems that identifies the material
                       region to which the contact angle applies (usually the 
                       liquid element block in solid/liquid conjugate problems).
====================== ==============================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = VELO_THETA_HOFFMAN NS   100   {45}   0.   1.   0.   0.   72.0   0   0   2

This condition applies a contact angle of 45 degrees between the free surface normal at
the 100 nodeset and the vector (0,1,0). The surface tension is 72 and the dynamic
contact angle is taken from element block 2. Normally, this latter vector is the normal
to the solid surface in contact with the free surface at this point.

-------------------------
**Technical Discussion**
-------------------------

* The constraint that is imposed at the nodeset node is as follows:

.. figure:: /figures/123_goma_physics.png
	:align: center
	:width: 90%

where :math:`v_{Hoffman}` is computed from the implicit solution of the Hoffman correlation;

.. figure:: /figures/124_goma_physics.png
	:align: center
	:width: 90%

or

.. figure:: /figures/125_goma_physics.png
	:align: center
	:width: 90%

where the Hoffman functions, fHoff and gHoff, which are inverses of each other are
given by;

.. figure:: /figures/126_goma_physics.png
	:align: center
	:width: 90%

.. figure:: /figures/127_goma_physics.png
	:align: center
	:width: 90%

When smoothing is not used, i.e. :math:`t_{relax}` = 0 , the imposed velocity is equal 
to that stipulated by the Hoffman correlation. Also see WETTING_SPEED_HOFFMAN
and SHARP_HOFFMAN_VELOCITY for level-set versions.

* For steady problems, the substrate velocity will be extracted from adjoining
  VELO_TANGENT, VELO_SLIP, or VELO_SLIP_ROT boundary conditions.

* Because the Hoffman functions are implicit, iteration is required in the
  determination of the wetting velocity. As a result, for very high capillary numbers,
  i.e. > :math:`10^6`, the iteration procedure in Goma may need to be modified.

* This condition was motivated by the Hoffman empirical correlation (cf. Stephan F.
  Kistler, “Hydrodynamics of Wetting,” in Wettability edited by John Berg, 1993).




.. TODO - Lines 84, 90, 97, and 101 have photos that needs to be replaced with the real equation.
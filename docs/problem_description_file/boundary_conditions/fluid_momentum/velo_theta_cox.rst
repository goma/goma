******************
**VELO_THETA_COX**
******************

::

	BC = VELO_THETA_COX NS <bc_id> <float_list> [integer]

-----------------------
**Description / Usage**
-----------------------

**(PCC/R_MOM_TANG1)**

This boundary condition card applies a dynamic contact angle dependent velocity
condition in a similiar fashion to VELO_THETA_TPL, but the functional form of the
velocity is different. The functional form stems from the hydrodynamic theory of
wetting by Cox.

The <float_list> for this boundary condition has eight values; definitions of the input
parameters are as follows:

=================== =============================================================
**VELO_THETA_COX**  Name of the boundary condition.
**NS**              Type of boundary condition (<bc_type>), where **NS** denotes
                    node set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (node set in
                    EXODUS II) in the problem domain.
<float1>            :math:`\theta`, equilibrium (static) contact angle subtended 
                    by wall normal and free surface normal, in units of degrees.
<float2>            :math:`n_x` , x-component of normal vector to the geometry
                    boundary (see important note below regarding variable wall
                    normals, viz. non-planar solid walls).
<float3>            :math:`n_y` , y-component of normal vector to the geometry
                    boundary. (see important note below regarding variable wall
                    normals, viz. non-planar solid walls).
<float4>            :math:`n_z` , z-component of normal vector to the geometry
                    boundary. (see important note below regarding variable wall
                    normals, viz. non-planar solid walls).
<float5>            :math:`\varepsilon_s` is the dimensionless slip length, i.e.
                    the  ratio of the slip
                    length to the characteristic length scale of the macroscopic
                    flow.
<float6>            :math:`\sigma` is the surface tension. This value is 
                    multiplied by the surface tension value stipulated by the 
                    surface tension material model.
<float7>            :math:`t_{relax}` is a relaxation time which can be used 
                    to smooth the imposed contact point velocity for transient
                    problems. Set to zero for no smoothing.
<float8>            :math:`V_{old}` is an initial velocity used in velocity 
                    smoothing for transient problems. Set to zero when smoothing
                    is not used.
[integer]           *blk_id*, an optional parameter that is the element block
                    number in conjugate problems that identifies the material
                    region to which the contact angle applies (usually the liquid
                    element block in solid/liquid conjugate problems).
=================== =============================================================

------------
**Examples**
------------

Following is a sample card:
::

    BC = VELO_THETA_COX NS   100   {45}   0.   1.   0. 0.1   72.0   0   0   2

This condition applies a contact angle of 45 degrees between the free surface normal at
the 100 nodeset and the vector (0,1,0). The surface tension is 72, the reciprocal of the
slip coefficient is 0.1, and the dynamic contact angle is taken from element block 2.
Normally, this latter vector is the normal to the solid surface in contact with the free
surface at this point.

-------------------------
**Technical Discussion**
-------------------------

* The constraint that is imposed at the nodeset node is as follows:

.. figure:: /figures/119_goma_physics.png
	:align: center
	:width: 90%

where :math:`v_{Cox}` is computed from

.. figure:: /figures/120_goma_physics.png
	:align: center
	:width: 90%

where the Cox functions, f and g, are given by;

.. figure:: /figures/121_goma_physics.png
	:align: center
	:width: 90%

* The parameters :math:`\lambda`, :math:`q_{inner}`, and :math:`q_{outer}` are currently not accessible from the input card
  and are hard-set to zero. :math:`\lambda` is the ratio of gas viscosity to liquid viscosity whereas
  :math:`q_{inner}` and :math:`q_{outer}` represent influences from the inner and outer flow regions

.. figure:: /figures/122_goma_physics.png
	:align: center
	:width: 90%

When smoothing is not used, i.e. :math:`t_{relax}` = 0 , the imposed velocity is equal to that
stipulated by the Hoffman correlation. Also see WETTING_SPEED_COX and
SHARP_COX_VELOCITY for level-set versions.

* For steady problems, the substrate velocity will be extracted from adjoining
  VELO_TANGENT, VELO_SLIP, or VELO_SLIP_ROT boundary conditions.

* The Cox wetting velocity requires evaluation of integrals for the function g
  ( :math:`\theta`, :math:`\lambda`)
  which is currently done numerically using 10-point Gaussian quadrature. As such
  the evaluation of the integrals is expected to become inaccurate as either 
  :math:`\theta_eq` tends
  toward zero or :math:`\theta` tends toward 180 degrees. Note that the integrand becomes
  singular as :math:`\theta` tends toward 0 or 180 degrees.

* This condition was motivated by the Cox hydrodynamic theory of wetting (cf.
  Stephan F. Kistler, “Hydrodynamics of Wetting,” in Wettability edited by John
  Berg, 1993 ).




.. TODO - Lines 81, 87, 93, and 101 have photos that needs to be replaced with the real equation.
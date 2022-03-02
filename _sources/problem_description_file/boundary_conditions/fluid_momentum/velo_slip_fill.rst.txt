******************
**VELO_SLIP_FILL**
******************

::

	BC = VELO_SLIP_FILL SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is applied only in problems involving embedded interface
tracking, that is, level set or volume of fluid. As in the case of the *VELO_SLIP* card, it
allows for slip to occur between fluid and solid substrate, but in this case slipping is
allowed only in a narrow region around the location of the interface where it intercepts
the solid boundary. Elsewhere, this boundary condition enforces a no-slip condition
between fluid and substrate.

When using the level set tracking, slip is allowed only near the intersection of the zero
level set contour and the substrate boundary, and then only in a region twice the level
set length scale wide centered on the zero level set. When using volume of fluid, the
criterion for slipping is that the absolute value of the color function should be less than
0.25.

This boundary condition is most often used in conjunction with the *FILL_CA* boundary
condition. The latter applies forces to contact lines in order to simulate wetting line
motion. These forces are applied in a weak sense to the same regions near the interface
so it is necessary to use *VELO_SLIP_FILL* with a large slipping coefficient so that
effectively no-slip is relaxed completely near the interface.

Definitions of the input parameters are as follows:

=================== ==============================================================
**VELO_SLIP_FILL**  Name of the boundary condition.
**SS**              Type of boundary condition (<bc_type>), where **SS**
                    denotes side set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (side set
                    in EXODUS II) in the problem domain.
<float1>            :math:`\beta`, the slip coefficient. The inverse of 
                    :math:`\beta` defines the
                    scaling between stress and slip. The parameter supplied
                    on the input deck is used only in the region define
                    above. Elsewhere, the slip coefficient is uniformly set to
                    :math:`10^{-6}`.
<float2>            :math:`v_{s,x}`, the x-component of surface velocity vector.
                    This would be the x-component of the fluid velocity if
                    a noslip condition were applied.
<float3>            :math:`v{s,y}`, the y-component of surface velocity vector. 
                    This would be the y-component of the fluid velocity 
                    if a noslip condition were applied.
<float4>            :math:`v{s,z}`, the z-component of surface velocity vector. 
                    This would be the z-component of the fluid velocity 
                    if a noslip condition were applied.
=================== ==============================================================

------------
**Examples**
------------

Following is a sample card without the optional parameters:
::

     BC = VELO_SLIP SS 10 100000.0 0.0 0.0 0.0

The large value of slip coefficient ensures nearly perfect slip in the region around the
interface.

-------------------------
**Technical Discussion**
-------------------------

* See the documentation under *VELO_SLIP* boundary condition for a description of
  the nature of this boundary condition.

* An important caveat when using this boundary condition to relax no-slip in the
  vicinity of the interface is that it relaxes all constraints on the velocities in the
  region. This includes the constraint to keep fluid from passing through the
  substrate boundary. For this region, it is usually also necessary to use a
  impenetrability condition, *VELO_NORMAL* for example, in conjunction with this
  boundary condition for appropriate results.



--------------
**References**
--------------

No References.
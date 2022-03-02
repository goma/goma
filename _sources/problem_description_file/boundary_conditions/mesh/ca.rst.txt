******
**CA**
******

::

	BC = CA NS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MESH)**

This boundary condition card applies a specified contact-angle on the mesh at a single
node nodeset. It is used exclusively in two dimensional computations. Its primary
application is imposing contact angles at static or dynamic contact lines. Consequently,
the nodeset is usually found where a free-surface boundary intersects a fixed,
“geometry” boundary.

The <float_list> for this boundary condition has four values; definitions of the input
parameters are as follows:

============= ==================================================================
**CA**        Name of the boundary condition.
**NS**        Type of boundary condition (<bc_type>), where **NS** denotes
              node set in the EXODUS II database.
<bc_id>       The boundary flag identifier, an integer associated with
              <bc_type> that identifies the boundary location (node set in
              EXODUS II) in the problem domain.
<float1>      θ, angle subtended by wall normal and free surface normal,
              in units of radians.
<float2>      n\ :sub:`x`, x-component of normal vector to the geometry
              boundary (see important note below regarding variable wall
              normals, viz. non-planar solid walls).
<float3>      n\ :sub:`y`, y-component of normal vector to the geometry
              boundary. (see important note below regarding variable wall
              normals, viz. non-planar solid walls).
<float4>      n\ :sub:`z`, z-component of normal vector to the geometry
              boundary. (see important note below regarding variable wall
              normals, viz. non-planar solid walls).
============= ==================================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = CA NS 100 1.4 0. 1. 0.

This condition applies a contact angle of 1.4 radians between the free surface normal at
the 100 nodeset and the vector (0,1,0). Normally, this latter vector is the normal to the
solid surface in contact with the free surface at this point.

-------------------------
**Technical Discussion**
-------------------------
* The constraint that is imposed at the nodeset node is:

	.. figure:: /figures/054_goma_physics.png
		:align: center
		:width: 90%


  where *n* is the normal to the geometry specified on the card itself, and n\ :sub:`fs` is the normal to the outward free surface computed internally by *Goma*. Also see the
  *CA_OR_FIX* card for an extension to this condition and *CA_EDGE* for its
  extension to three dimensions.

* In addition for the case in which the geometry normal components are set to zero, 
  the wall normal is allowed to vary with a geometrical boundary condition, i.e.,
  *GD_TABLE, SPLINE, PLANE*, etc. The geometry normal is found on the same or
  on a neighboring element that contains the dynamic contact angle in question. If a
  *GD*_ type boundary condition is used to describe the wall (i.e., *GD_TABLE*), one
  must specify the *R_MESH_NORMAL* equation type on that equation for the
  variable wall normal to take effect.

* **Important: Variable Wall Normals**. Situations for which the wall shape is nonplanar,
  meaning that the normal vector is not invariant as the contact line moves,
  there is an option to leave all of the normal-vector components zero. In this case
  *Goma* then seeks to determine the local wall normal vector from the geometry it is
  currently on, using the element facets. It is recommended that this option not be
  used unless the geometry is truly nonplanar, as the logic is complex and not 100%
  reliable. An example of such a case is as follows:

::

     BC = CA NS 100 1.4 0. 0. 0.

|

  Notice how all three components of the normal vector are set to zero.

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

.. figure:: /figures/055_goma_physics.png
	:align: center
	:width: 90%



--------------
**References**
--------------

No References.
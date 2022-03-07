********************
**FLOW_HYDROSTATIC**
********************

::

	BC = FLOW_HYDROSTATIC SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition allows the user to impose a pressure force that varies linearly
with position over the boundary. It functions in much the same manner as the
*FLOW_PRESSURE* boundary condition except that more variability is allowed in the
imposed pressure. As the name implies, this boundary condition is most often used to
impose hydrostatic pressure profiles in problems in which gravitational forces play a
role.

The <float_list> has four values to be specified; definitions of the input parameters 
are as follows:

===================== ==================================================================
**FLOW_HYDROSTATIC**  Boundary condition name.
**SS**                Type of boundary condition (<bc_type>), where **SS**
                      denotes side set in the EXODUS II database.
<bc_id>               The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location (side set
                      in EXODUS II) in the problem domain.
<float1>              :math:`\delta P_x`, the pressure variation in x-direction.
<float2>              :math:`\delta P_y`, the pressure variation in y-direction.
<float3>              :math:`\delta P_z`, the pressure variation in z-direction.
<float4>              :math:`P_0`, the pressure value at the coordinate point (0,0,0).
                      This serves as a means of establishing a datum and it is
                      not required that (0,0,0) lie on the sideset.
===================== ==================================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = FLOW_HYDROSTATIC SS 15   0.0 0.0 -1.5 10.0

This card will impose a pressure profile on side set 15 so that the pressure decreases by
1.5 as the z coordinate increases by one unit. At the point, (0,0,0) the pressure imposed
is 10.0. Note that (0,0,0) does not necessarily have to be on side set 15.

-------------------------
**Technical Discussion**
-------------------------

* The mathematical form of the boundary condition imposed by this card is as
  follows:

.. figure:: /figures/098_goma_physics.png
	:align: center
	:width: 90%

where *n* is the outward normal vector to the boundary, *T* is the total fluid stress
tensor, and *x, y, z* are the global coordinate positions.

* Like the *FLOW_PRESSURE* conditions, this is a weakly integrated condition and
  the comments appearing with that card apply equally well here.

* Most often this boundary condition is used in problems in which gravity is present.
  Under these circumstances, the pressure profile across a fully-developed flow inlet
  is not constant but varies according to hydrostatic head. Hence, the
  *FLOW_PRESSURE* condition cannot be used to provide the inlet pressure.
  Instead, this card is used with the variation in the pressure being imposed
  according to the direction of gravity. Thus, some if not all of 
  :math:`\delta P_x`, :math:`\delta P_y`, or :math:`\delta P_z` will
  be functions of gravity and the fluid density.

* It is true that this variation could be determined automatically by *Goma* from its
  known values for density and gravitational direction. But for a variety of reasons,
  this may not always be the best option. Instead, the user is allowed to vary the
  pressure on a boundary independently of the density and gravitational forces set
  elsewhere in the material file. If consistency is important in the problem at hand,
  then the user is cautioned to be consistent.

* The input parameter :math:`P_0` as noted above serves as a datum to the relationship. 
  In theory, it is the pressure value that would be computed at the point (0,0,0), but 
  in reality it is chosen to impose a known pressure at some point in the domain.



--------------
**References**
--------------

No References.
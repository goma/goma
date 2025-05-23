*******************
**VELO_TANGENT_3D**
*******************

::

	BC = VELO_TANGENT_3D SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MOMENTUM)**

This boundary condition is the three dimensional analog of the *VELO_TANGENT*
condition. It is used to strongly set the tangential velocity component along a side set in
a three-dimensional problem. It is not a completely general condition since it can set
only a single tangential velocity component. It can only be applied to flat surfaces or
surfaces which have only one radius of curvature such as a cylinder.

The <float_list> requires four values be specified; a description of the input parameters follows:

=================== =================================================================
**VELO_TANGENT_3D** The name of the boundary condition.
**SS**              Type of boundary condition (<bc_type>), where **SS**
                    denotes side set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (side set
                    in EXODUS II) in the problem domain.
<float1>            :math:`v_t`, the value assigned to the tangential velocity
                    component.
<float2>            :math:`t_x`, the x-component of a unit normal vector tangent to
                    the surface; this vector must be tangent at all points on
                    the surface. The direction of the imposed tangential
                    velocity component is n × t with *n* the outward-pointing
                    normal.
<float3>            :math:`t_y`, the y-component of a unit normal vector tangent to
                    the surface; this vector must be tangent at all points on
                    the surface. The direction of the imposed tangential
                    velocity component is n × t with *n* the outwardpointing
                    normal.
<float4>            :math:`t_z`, the z-component of a unit normal vector tangent to
                    the surface; this vector must be tangent at all points on
                    the surface. The direction of the imposed tangential
                    velocity component is n × t with *n* the outwardpointing
                    normal.
=================== =================================================================

------------
**Examples**
------------

The following is an example of the card:
::

     BC = VELO_TANGENT_3D SS 10   1.0   0.0 0.0 1.0

One could use this card to set the tangential velocity on a cylindrically shaped side set
10 provided that the cylinders axis was parallel to the z-axis. In this fashion, the
tangential velocity component perpendicular to the z-axis is set to 1.0.

-------------------------
**Technical Discussion**
-------------------------

* The constraint applied to the velocity vector by this condition on the side set is:

.. figure:: /figures/088_goma_physics.png
	:align: center
	:width: 90%

where :math:`\tilde{t}` = n × t with the components of *t* supplied on the card. The advantages of
introducing the normal vector is that it permits use of this card on curving surfaces
provided the curvature occurs in only one direction and a single tangent vector
exists that is perpendicular to both the surface normal and the direction of
curvature. This of course implies that the tangential component can only be
applied in the direction of the curvature.

* Such conditions are of course met by a planar surface, but also a cylindrical
  surface. In the latter case, the vector *t* should be parallel to the axis of the cylinder.
  One application for this condition is in three-dimensional eccentric roll coating in
  which the roll speed can be set using this condition. The axis vectors of both roll
  coaters are supplied on the card.




.. TODO - In line 68, the photo needs to be repalced by the proper equation.
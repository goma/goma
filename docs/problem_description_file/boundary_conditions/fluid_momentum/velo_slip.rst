*************
**VELO_SLIP**
*************

::

	BC = VELO_SLIP SS <bc_id> <float_list> [integer1] [float5]

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition allows for slip between the fluid and a boundary using an
implementation of the Navier slip relation. This relation fixes the amount of slip as a
function of the applied shear stress. The scaling between stress and slip is a user
parameter. This implementation also permits (in two dimensions only) variable scaling
dependent upon distance from a mesh node. The latter can be used in modeling
dynamic contact lines. This condition cannot currently be used on connecting surfaces.

There are four required values in <float_list> and two optional values; definitions of
the input parameters are as follows:

============== ==================================================================
**VELO_SLIP**  Name of the boundary condition (<bc_name>).
**SS**         Type of boundary condition (<bc_type>), where **SS** denotes
               side set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set in
               EXODUS II) in the problem domain.
<float1>       :math:`\beta`, the slip coefficient. The inverse of :math:`\beta` 
               defines the scaling
               between stress and slip. Hence, for small values of 
               :math:`\beta`, large
               shear stresses are needed for a given amount of slip, and
               conversely, for large values of :math:`\beta`, 
               the amount of stress
               needed for the same degree of slip decreases (see below for
               a more rigorous description).
<float2>       :math:`v_{s,x}`, the x-component of surface velocity vector. This would
               be the x-component of the fluid velocity if a no slip
               condition were applied.
<float3>       :math:`v_{s,y}`, the y-component of surface velocity vector. This would
               be the y-component of the fluid velocity if a no slip
               condition were applied.
<float4>       :math:`v_{s,z}`, the z-component of surface velocity vector. This would
               be the z-component of the fluid velocity if a no slip
               condition were applied.
[integer]      :math:`N_{cl}`, a single-node node set identification number. When the
               variable coefficient slip relation is used, distance is 
               measured relative to this node (see discussion below).
               Normally, this node set represents the location of the
               dynamic contact line. Note that this option is generally only
               used in two-dimensional simulations.
[float5]       :math:`\alpha`, the distance scale in the variable slip model (see the
               discussion below). Both :math:`N_{cl}` and :math:`\alpha` should be present to
               activate the variable slip model.
============== ==================================================================

------------
**Examples**
------------

Following is a sample card without the optional parameters:
::

     BC = VELO_SLIP SS 10 0.1 0.0 0.0 0.0

-------------------------
**Technical Discussion**
-------------------------

* The general form of this boundary condition is

.. figure:: /figures/089_goma_physics.png
	:align: center
	:width: 90%

where :math:`\tau` is the deviatoric portion of the fluid stress tensor, :math:`\beta` is the Navier slip
coefficient and :math:`v_s` is the velocity of the solid surface. The velocity of the surface
must be specified, as described in the Description/Usage subsection above. It is a
weakly integrated vector condition, as noted above, so it will be added to each of
the three momentum equation components.

This last point is important to keep in mind, especially when applying this
condition to boundaries that are not parallel to any of the principle axes. It is
possible under these circumstances that this condition will allow motion through a
boundary curve in addition to slip tangential to it. This can be avoided by including
a rotated boundary condition like *VELO_NORMAL* on the same sideset. This will
cause the momentum equations to be rotated to normal and tangential components
and also enforce no normal flow of the material. Whatever slipping that takes place
will be in the tangential direction.

* The variable slip coefficient model is quite simple: 
  **EQUATION**, where *d* is
  the absolute distance from node :math:`N_{cl}` identified on the card; the coefficients :math:`\beta` and :math:`\alpha`
  are also supplied on input. This relation is protected against overflowing as *d*
  increases. This model can be used to allow slipping to occur in a region close to the
  node set, but at points further removed, a no slip boundary (:math:`\beta` large) is reinstated on the sideset.




.. TODO - In line 76, the photo needs to be repalced by the proper equation. Also, In line 96 the equation need to be written so it read correctly where it reads "EQUATION". 
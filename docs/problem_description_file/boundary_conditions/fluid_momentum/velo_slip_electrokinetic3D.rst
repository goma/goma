******************************
**VELO_SLIP_ELECTROKINETIC3D** 
******************************

::

	BC = VELO_SLIP_ELECTROKINETIC3D SS <bc_id> [floatlist]

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MOMENTUM)**

This is a 3D generalization of the VELO_SLIP_ELECTROKINETIC boundary
condition. It is similar to VELO_TANGENT_3D except the slip velocity is calculated
based on Helmholtz-Smulkowski relation. This boundary condition allows for slip
between the fluid and a solid boundary due to electrokinetic effects on the charged
solid wall. The user provides the following parameters: zeta potential at the wall,
permittivity of the fluid and.

============================== ========================================================
**VELO_SLIP_ELECTROKINETIC3D** Name of the boundary condition (<bc_name>).
**SS**                         Type of boundary condition (<bc_type>), where **SS**
                               denotes side set in the EXODUS II database.
<bc_id>                        The boundary flag identifier, an integer associated with
                               <bc_type> that identifies the boundary location 
                               (side set in EXODUS II) in the problem domain.
<float1>                       :math:`\varepsilon`, absolute permittivity of the fluid.
<float2>                       :math:`\zeta`, the surface potential of solid boundary. 
                               It is referred to as the zeta potential.
<float3>                       :math:`t_x`, the x-component of a unit normal vector 
                               tangent to the surface; this vector must be tangent at 
                               all points on
                               the surface. The direction of the imposed tangential
                               velocity component is n × t with *n* the outwardpointing
                               normal.
<float4>                       :math:`t_y`, the y-component of a unit normal vector 
                               tangent to the surface; this vector must be tangent at all points on
                               the surface. The direction of the imposed tangential
                               velocity component is n × t with *n* the outwardpointing
                               normal.
<float5>                       :math:`t_z`, the z-component of a unit normal vector 
                               tangent to the surface; this vector must be tangent at all points on
                               the surface. The direction of the imposed tangential
                               velocity component is n × t with *n* the outwardpointing
                               normal.
============================== ========================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = VELO_SLIP_ELECTROKINETIC3D SS 10 1.e-5 1.e-2 0. 0. 1.

-------------------------
**Technical Discussion**
-------------------------

* The general form of this boundary condition is

.. figure:: /figures/091_goma_physics.png
	:align: center
	:width: 90%

where :math:`\varepsilon` is the absolute permittivity of the medium, :math:`\zeta` is the zeta potential, :math:`E_t` is the
electric field tangent to the solid surface, and :math:`v_s` is the slip velocity.




.. TODO - Line 65 contains a photo that needs to be exchanged for the equation.
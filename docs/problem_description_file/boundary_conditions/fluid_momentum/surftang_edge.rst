*****************
**SURFTANG_EDGE** 
*****************

::

	BC = SURFTANG_EDGE SS <bc_id1> <bc_id2> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition card is used to insert surface tension forces on an edge
boundary defined by the primary and secondary sidesets. The direction of application
of the surface tension ( :math:`\sigma` ) is specified by the vector (defined by 
(< :math:`m_x` >, < :math:`m_y` >, < :math:`m_z` >)).
This card is the three-dimensional analog of the *CAP_ENDFORCE* card. It is often
used at free-surface outflow boundaries if the outflow velocity is not set by a strong
condition. This condition is an unrotated, weak integrated vector condition.

There are four values to be supplied in the <float_list>; definitions of the input
parameters are as follows:

================= =================================================================
**SURFTANG_EDGE** Name of the boundary condition.
**SS**            Type of boundary condition (<bc_type>), where **SS**
                  denotes side set in the EXODUS II database.
<bc_id1>          The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the primary boundary location
                  (side set in EXODUS II) in the problem domain. This
                  side set is usually attached to a free surface.
<bc_id2>          The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the secondary boundary
                  location (side set in EXODUS II) in the problem
                  domain. The boundary condition is applied on the edge
                  defined by the intersection of this side set with the
                  primary side set.
<float1>          :math:`m_x`, the x-component of direction of application of
                  surface tension force.
<float2>          :math:`m_y`, the y-component of direction of application of
                  surface tension force.
<float3>          :math:`m_z`, the z-component of direction of application of
                  surface tension force.
<float4>          a factor multiplying the surface tension value read from
                  the material file when evaluating the surface integral
                  imposed by this boundary condition.
================= =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = SURFTANG_EDGE SS 80 60   0. -1. 0. 1.

-------------------------
**Technical Discussion**
-------------------------

* The need for this boundary condition appears out of the formulation used to apply
  capillary forces to surfaces in three-dimensions. The surface divergence theorem is 
  used to simplify the curvature term in the capillary stress jump condition. This
  produces integrals of the form:

.. figure:: /figures/108_goma_physics.png
	:align: center
	:width: 90%

where *C* is the bounding curve of the capillary free surface, :math:`\sigma` is the surface
tension, :math:`\phi_i` is a finite element shape function and **m** is a vector that is at once
normal to the capillary surface and also normal to the curve *C*. It always points
outward from the domain in question.

Most often this boundary condition appears at outflow boundaries of free-surfaces.
It is applied along the edge where the free-surface intercepts the outflow plane. In
this case, the **m** vector is normal to the outflow plane. If the outflow velocity is not
strongly set by a Dirichlet condition or other strongly enforced condition, this
boundary condition needs to be present so that a proper inclusion of all relevant
surface tension terms is performed.

* The <factor> parameter is provided to allow the user to independently vary the
  surface tension value associated with this term alone. The value for :math:`\sigma` used in the
  preceding expression is the surface tension value obtained from the model
  specified in the material file multiplied by the value of <float>. Reasons for doing
  this are somewhat obscure but important to the practitioners of this art.




.. TODO - Lines 69 has a photo that needs to be replaced with the real equation.
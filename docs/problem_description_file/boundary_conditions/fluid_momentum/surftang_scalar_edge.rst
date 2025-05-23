************************
**SURFTANG_SCALAR_EDGE**
************************

::

	BC = SURFTANG_SCALAR_EDGE SS <bc_id1> <bc_id2> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

Like the *SURFTANG_EDGE* card, this boundary condition card is used to insert
surface tension forces on an outflow edge boundary defined by the primary and
secondary sidesets. In contrast to the *SURFTANG_EDGE* card, the direction of
application of the surface tension ( :math:`\sigma` ) is predetermined automatically as the binormal
along the edge with respect to the outward facing normal of the primary sideset. This
condition is also an unrotated, weak integrated vector condition. It should be used 
only in three-dimensional applications.

Definitions of the input parameters are as follows:

======================== ===========================================================
**SURFTANG_EDGE_SCALAR** Name of the boundary condition.
**SS**                   Type of boundary condition (<bc_type>), where **SS**
                         denotes side set in the EXODUS II database. Since it is
                         an EDGE condition, it applies to a curve defined as the
                         intersection of the primary and secondary sideset.
<bc_id1>                 The boundary flag identifier, an integer associated with
                         <bc_type> that identifies the primary boundary location
                         (side set in EXODUS II) in the problem domain. This
                         side set is used in defining the edge and the local vector
                         basis (normal, tangent, binormal) and is usually attached
                         to a free surface.
<bc_id2>                 The boundary flag identifier, an integer associated with
                         <bc_type> that identifies the secondary boundary
                         location (side set in EXODUS II) in the problem
                         domain. It is used in defining the edge and the local
                         vector basis (normal, tangent, binormal).The boundary
                         condition is applied on the edge defined by the
                         intersection of this side set with the primary side set.
<float>                  A factor multiplying the surface tension value read from
                         the material file when evaluating the surface integral
                         imposed by this boundary condition.
======================== ===========================================================

------------
**Examples**
------------

The following is a sample input card:
::

    BC = SURFTANG_EDGE_SCALAR SS 5 10 1.0

applies the boundary integral (see the Technical Discussion) along the curve described
by the intersection of side sets 5 and 10. The value for surface tension in the material
file is used unmodified since the multiplying factor is 1.0.

-------------------------
**Technical Discussion**
-------------------------

* The need for this boundary condition appears out of the formulation used to apply
  capillary forces to surfaces in three-dimensions. The surface divergence theorem is
  used to simplify the curvature term in the capillary stress jump condition. This
  produces integrals of the form:

.. figure:: /figures/110_goma_physics.png
	:align: center
	:width: 90%

where *C* is the bounding curve of the capillary free surface, :math:`\sigma` is the surface
tension, :math:`\phi_i` is a finite element shape function and **m** is the outward binormal vector
to the curve *C* with respect to the normal of the primary side set.

Most often this boundary condition appears at outflow boundaries of free surfaces.
It is applied along the edge where the free surface intercepts the outflow plane. If
the outflow velocity is not strongly set by a Dirichlet condition or other strongly
enforced condition, this boundary condition needs to be present so that a proper
inclusion of all relevant surface tension terms is performed.

* The <factor> parameter is provided to allow the user to independently vary the
  surface tension value associated with this term alone. The value for :math:`\sigma` used in the
  preceding expression is the surface tension value obtained from the model
  specified in the material file multiplied by the value of <float>. Reasons for doing
  this are somewhat obscure but important to the practitioners of this art.




.. TODO - Lines 71 has a photo that needs to be replaced with the real equation.
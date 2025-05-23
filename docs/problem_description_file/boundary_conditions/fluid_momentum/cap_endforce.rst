****************
**CAP_ENDFORCE** 
****************

::

	BC = CAP_ENDFORCE NS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SPECIAL/VECTOR MOMENTUM)**

This boundary condition card adds surface tangent forces to the momentum equations
at the endpoint of a free-surface. There are four values to be input for the <
float_list>; definitions of the input parameters are as follows:

================ ==========================================================
**CAP_ENDFORCE** Name of the boundary condition (<bc_name>).
**NS**           Type of boundary condition (<bc_type>), where **NS**
                 denotes node set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (node
                 set in EXODUS II) in the problem domain.
<float1>         X-component of surface tangent vector at end point.
<float2>         Y-component of surface tangent vector at end point.
<float3>         Z-component of surface tangent vector at end point.
<float4>         Equilibrium surface tension value. (See Technical
                 Discussion.)
================ ==========================================================

This condition need only be applied at the intersection of outflow or inflow surfaces
and the free-surface. The sign on the tangent vector depends on whether the computed
tangent vector is facing inward or outward. This can be figured by t = n × k.

------------
**Examples**
------------

The following is a sample input card using several APREPRO variables:
::

     BC = CAP_ENDFORCE NS 100   {sind(th2)} {-cosd(th2)} 0.0 {surf_tens}

-------------------------
**Technical Discussion**
-------------------------

* The need for this boundary condition appears out of the formulation used to apply
  capillary forces to surfaces. The surface divergence theorem is used to simplify the
  curvature term in the capillary stress jump condition. This produces integrals of the
  form:

.. figure:: /figures/107_goma_physics.png
	:align: center
	:width: 90%

where *C* is the bounding curve of the capillary free surface, :math:`\sigma` is the surface
tension, :math:`\phi_i` is a finite element shape function and *m* is a vector that is at once
normal to the capillary surface and also normal to the curve *C*. It always points
outward from the domain in question. While this is completely general for threedimensions,
a surface can be reduced to a curve for two-dimensions and the
divergence theorem still applies (for this boundary condition).

* This card or the *CAP_ENDFORCE_SCALAR* is used in conjunction with the
  *CAPILLARY* card to complete (as indicated above) the treatment of capillarity
  conditions. It is only required when an inflow or outflow boundary intersects a free
  surface.

* The *CAP_ENDFORCE* boundary condition is applied through function
  fapply_ST (in file *mm_ns_bc.c*). The boundary term is computed as the product
  of the surface tension supplied on this card (<float4>) and the value supplied on
  the *Surface Tension* card in the material file. When the latter card is missing, 
  *Goma* defaults its value to 1.0.

* This card was previously called *SURFTANG* for the surface tangent component of
  the capillary force. Old input decks can be updated simply by changing the name
  of the boundary condition without changing the parameters.




.. TODO - Lines 55 has a photo that needs to be replaced with the real equation.
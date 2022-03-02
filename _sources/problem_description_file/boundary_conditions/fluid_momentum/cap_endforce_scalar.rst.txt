***********************
**CAP_ENDFORCE_SCALAR**
***********************

::

	BC = CAP_ENDFORCE_SCALAR NS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(SPECIAL/VECTOR MOMENTUM)**

This boundary condition card is very similar to the *CAP_ENDFORCE* card. It adds
surface tangent forces to the momentum equations at the endpoint of a free-surface, but
does not require specification of the surface tangent vector. The current free-surface
tangent vector is used as the surface tangent vector. Definitions of the input 
parameters are as follows:

======================== ======================================================
**CAP_ENDFORCE_SCALAR**  Name of the boundary condition (<bc_name>).
**NS**                   Type of boundary condition (<bc_type>), where **NS**
                         denotes node set in the EXODUS II database.
<bc_id>                  The boundary flag identifier, an integer associated with
                         <bc_type> that identifies the boundary location (node
                         set in EXODUS II) in the problem domain.
<float>                  Equilibrium surface tension value. (See Technical
                         Discussion.)
======================== ======================================================

This condition need only be applied at the intersection of outflow or inflow surfaces
and the free-surface.

------------
**Examples**
------------

The following is a sample input card:
::

    BC = CAP_ENDFORCE_SCALAR NS 100   60.0

-------------------------
**Technical Discussion**
-------------------------

* The need for this boundary condition appears out of the formulation used to apply
  capillary forces to surfaces. The surface divergence theorem is used to simplify the
  curvature term in the capillary stress jump condition. This produces integrals of the
  form:

.. figure:: /figures/109_goma_physics.png
	:align: center
	:width: 90%

where *C* is the bounding curve of the capillary free surface, :math:`\sigma` is the surface
tension, :math:`\phi_i` is a finite element shape function and **m** is a vector that is at once
normal to the capillary surface and also normal to the curve *C*. It always points
outward from the domain in question. While this is completely general for threedimensions,
a surface can be reduced to a curve for two-dimensions and the
divergence theorem still applies (for this boundary condition).

* This card or the *CAP_ENDFORCE* is used in conjunction with the *CAPILLARY*
  card to complete (as indicated above) the treatment of capillarity conditions. It is
  only required when an inflow or outflow boundary intersects a free surface.

* The *CAP_ENDFORCE_SCALAR* boundary condition is applied through function
  fapply_ST_scalar (in file *mm_ns_bc.c*). The boundary term is computed as the
  product of the surface tension supplied on this card (<float>) and the value
  supplied on the *Surface Tension* card in the material file. When the latter card is
  missing, *Goma* defaults its value to 1.0.

* This card was previously called *SURFTANG_SCALAR* for the surface tangent
  component of the capillary force. Old input decks can be updated simply by
  changing the name of the boundary condition without changing the parameters.




.. TODO - Lines 53 has a photo that needs to be replaced with the real equation.
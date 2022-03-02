***********************
**POROUS_LIQ_PRESSURE**
***********************

::

	BC = POROUS_LIQ_PRESSURE NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/POR_LIQ_PRES)**

This Dirichlet boundary condition is used to set the liquid phase pore pressure at a noden set. It can be applied to a node set on a boundary of a *POROUS_SATURATED,
POROUS_UNSATURATED* or *POROUS_TWO_PHASE* medium type (see *Media Type*
card).

======================= ========================================================
**POROUS_LIQ_PRESSURE** Boundary condition name (bc_name).
**NS**                  Type of boundary condition (<bc_type>), where **NS**
                        denotes node set in the EXODUS II database.
<bc_id>                 The boundary flag identifier, an integer associated with
                        <bc_type> that identifies the boundary location (node
                        set in EXODUS II) in the problem domain.
<float1>                Value of liquid phase pressure.
[float2]                An optional parameter (that serves as a flag to the code
                        for a Dirichlet boundary condition). If a value is present,
                        and is not -1.0, the condition is applied as a residual
                        equation. Otherwise, it is a "hard set" condition and is
                        eliminated from the matrix. *The residual method must
                        be used when this Dirichlet boundary condition is used
                        as a parameter in automatic continuation sequences*.
======================= ========================================================

------------
**Examples**
------------

The boundary condition card
::

   BC = POROUS_LIQ_PRESSURE NS 101 {pcmin}

sets the porous liquid pressure at the boundary denoted by node set 101 to the value
represented by the APREPRO variable {pcmin}.

-------------------------
**Technical Discussion**
-------------------------

Setting the porous liquid pressure to a value cannot be done independently of the
saturation as the two are related through the vapor pressure curve for simulations in
partially saturated media (see *Saturation* model card). Keep in mind that when using
this card in these situations, you are setting also the saturation level based on the
capillary pressure, defined as :math:`p_{gas}` - :math:`p_{liq}` = :math:`p_c` . The convention in *Goma* is that when
the capillary pressure :math:`p_c` is greater than zero, the saturation level is less than unity, viz.
the medium is partially saturated. When :math:`p_c` is less than zero, i.e., when the liquid
phase pressure is greater than the gas phase pressure, then the medium is saturated (in
this case the capillary pressure is poorly defined, though). Also, for *Media Type* options
of *POROUS_UNSATURATED*, the ambient gas pressure is constant within the pore
space and is set by the *Porous Gas Constants* card in the material file. This boundary
condition, when setting the liquid phase pressure, must be used with consideration of
these definitions.

For saturated media (viz. *Media Type* of *POROUS_SATURATED*), this discussion is
not relevant. In this case, one must only consider the pressure level as it may effect the
isotropic stress in poroelastic problems.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk
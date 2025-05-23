***********************
**POROUS_GAS_PRESSURE**
***********************

::

	BC = POROUS_GAS_PRESSURE NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/POR_GAS_PRES)**

This Dirichlet boundary condition is used to set the gas-phase pore pressure at the
boundary of a *POROUS_TWO_PHASE* medium type (see *Media Type* card). This
condition makes no sense on other *POROUS Media Types*; the gas pressure in those
cases is constant and set using the *Porous Gas Constants* card (*Microstructure
Properties*).

======================= =============================================================
**POROUS_GAS_PRESSURE** Boundary condition name (<bc_name>).
**NS**                  Type of boundary condition (<bc_type>), where **NS**
                        denotes node set in the EXODUS II database.
<bc_id>                 The boundary flag identifier, an integer associated with
                        <bc_type> that identifies the boundary location (node
                        set in EXODUS II) in the problem domain.
<float>                 Value of gas phase pressure.
[float2]                An optional parameter (that serves as a flag to the code
                        for a Dirichlet boundary condition). If a value is present,
                        and is not -1.0, the condition is applied as a residual
                        equation. Otherwise, it is a “hard set” condition and is
                        eliminated from the matrix. *The residual method must 
                        be used when this Dirichlet boundary condition is used
                        as a parameter in automatic continuation sequences*.
======================= =============================================================

------------
**Examples**
------------

The boundary condition card
::

   BC = POROUS_GAS_PRESSURE NS 101 {pgas}

sets the porous gas pressure at the boundary denoted by node set 101 to the value
represented by the APREPRO variable {pgas}.

-------------------------
**Technical Discussion**
-------------------------

Setting the porous liquid pressure to a value cannot be done independently of the
saturation as the two are related through the vapor pressure curve for simulations in
partially saturated media (see *Saturation* model card). Keep in mind that when using
this card in these situations, you are setting also the saturation level based on the
capillary pressure, defined as :math:`p_{gas}` - :math:`p_{liq}` = :math:`p_c` . The convention in Goma is that when
the capillary pressure :math:`p_c` is greater than zero, the saturation level is less than unity, viz.the medium is partially saturated. When :math:`p_c` is less than zero, i.e., when the liquid phase pressure is greater than the gas phase pressure, then the medium is saturated (in this case the capillary pressure is poorly defined, though). Also, this pressure sets the datum of pressure for deformable porous media and must be set in a manner compatible with the solid-stress values on the boundaries of the porous matrix.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk
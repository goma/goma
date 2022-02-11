********
**PSPG**
********

::

	BC = PSPG SS <bc_id>

-----------------------
**Description / Usage**
-----------------------

**(WIC/CONTINUITY)**

This special type of boundary condition exists for pressure-stabilized incompressible
flow simulations only. This card should be used only if the value of the *Pressure
Stabilization* card has been set to *yes*. In conjunction with this feature, equal-order
interpolation should be used for the velocity and pressure. If *PSPG* is used, a boundary
integral will be added to the continuity equation to represent the gradients of velocity in the momentum residual, which has been added onto the continuity equation for
stabilization. This term is only needed on inflow and outflow boundaries; in the rest of
the domain, it cancels out. For more details about the derivation of this term, see the
paper by Droux and Hughes (1994).

This boundary condition card requires no integer or floating point constants.
Definitions of the input parameters are as follows:

=========== =============================================================
**PSPG**    Name of the boundary condition (<bc_name>).
**SS**      Type of boundary condition (<bc_type>), where SS denotes
            side set in the EXODUS II database.
<bc_id>     The boundary flag identifier, an integer associated with
            <bc_type> that identifies the boundary location (side set in
            EXODUS II) in the problem domain.
=========== =============================================================

------------
**Examples**
------------

The following is an example of using this card on both the inflow and outflow planes of
the domain.
::

   BC = PSPG   SS   40

::

   BC = PSPG   SS   20

-------------------------
**Technical Discussion**
-------------------------

Please see Rao (1996) memo for a more detailed discussion of pressure stabilization
and its implementation in *Goma*.



--------------
**References**
--------------

GTM-001.0: Pressure Stabilization in Goma using Galerkin Least Squares, July 17,
1996, R. R. Rao

Droux, J. J. and T. J. R. Hughes, “A Boundary Integral Modification of the Galerkin
Least Squares Formulation for the Stokes Problem, ” Comput. Methods Appl. Mech.
Engrg., 113 (1994) 173-182.


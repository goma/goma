***************
**Q_VELO_SLIP**
***************

::

	BC = Q_VELO_SLIP SS <bc_id>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition card is used to calculate the surface integral for viscous
heating due to slip in the tangential velocity component on a surface. Definitions of 
the input parameters are as follows:

================ =================================================================
**Q_VELO_SLIP**  Name of the boundary condition (<bc_name>).
**SS**           Type of boundary condition (<bc_type>), where **SS** denotes
                 side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set in
                 EXODUS II) in the problem domain.
================ =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = Q_VELO_SLIP_BC SS 10

-------------------------
**Technical Discussion**
-------------------------

Use of this boundary condition requires specification of the slip velocity components
by using either the *VELO_SLIP* or *VELO_SLIP_ROT* boundary condition.




**************
**POROUS_GAS**
**************

::

	BC = POROUS_GAS SS <bc_id> <integer_list> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC/POR_LIQ_PRES)**

This boundary condition card is used to equate flux of solvent in the porous medium
and external gas. The condition is similar to the solid-liquid interface conditions that
apply to interfaces between a porous medium and an external gas (in which the energy
equation is used to solve for solvent concentration in the gas phase). This boundary
condition is still in development.

There are three values in the <integer_list> and two values in the <float_list> for which to supply values; definitions of the input parameters are as follows:

=============== ===============================================================
**POROUS_GAS**  Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS**
                denotes side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set
                in EXODUS II) in the problem domain.
<integer1>      Element block ID of solid, porous phase from the
                EXODUS II database.
<integer2>      Element block ID of gas phase from the EXODUS II
                database.
<integer3>      Species number of liquid phase in porous medium.
<float1>        Vapor density.
<float2>        Factor to allow normal velocity in gas.
=============== ===============================================================

------------
**Examples**
------------

Users are referred to the Cairncross (1999) reference for the best example of card
usage.

-------------------------
**Technical Discussion**
-------------------------

This highly specialized boundary condition is best explained in a paper by Cairncross
(1999).



--------------
**References**
--------------

GTM-028.0: Modeling Drying of Dip-Coated Films with Strongly-Coupled Gas Phase
Natural Convection, R. A. Cairncross, 1999.

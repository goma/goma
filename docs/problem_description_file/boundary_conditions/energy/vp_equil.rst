************
**VP_EQUIL**
************

::

	BC = VP_EQUIL SS <bc_id> <integer_list> <float>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ENERGY)**

This boundary condition card is used to equate solvent partial pressure in the gas
between the porous medium and the external phase. The condition is similar to the
solid-liquid interface conditions that apply to interfaces between a porous medium and
an external gas phase (in which the energy equation is used to solve for solvent
concentration in the gas phase). This boundary condition is still under development.

There are three values to be specified for the <integer_list>; definitions of the input
parameters are as follows:

============ ==================================================================
**VP_EQUIL** Name of the boundary condition (<bc_name>).
**SS**       Type of boundary condition (<bc_type>), where **SS** denotes
             side set in the EXODUS II database.
<bc_id>      The boundary flag identifier, an integer associated with
             <bc_type> that identifies the boundary location (side set in
             EXODUS II) in the problem domain.
<integer1>   Element block ID of solid, porous phase from the EXODUS
             II database.
<integer2>   Element block ID of gas phase from the EXODUS II database.
<integer3>   Species number of liquid phase in porous medium.
<float>      Ambient pressure in external gas phase. 
============ ==================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = VP_EQUIL SS 100   1 2 0 0.0

where the solid/porous phase is present in element block 1 and the gas phase is present
in element block 2. The external gas phase pressure has been set to 0.0.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.





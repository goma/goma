*****************
**CONT_TANG_VEL**
*****************

::

	BC = CONT_TANG_VEL SS <bc_id>

-----------------------
**Description / Usage**
-----------------------

**(SIC/MOMENTUM)**

This boundary condition card enforces continuity of tangential velocity between two
phases with discontinuous velocity treatment. The condition only applies to interphase
mass, heat, and momentum transfer problems with discontinuous (or multivalued)
variables at an interface, and it must be invoked on fields that employ the **Q1_D** or
**Q2_D** interpolation functions to “tie” together or constrain the extra degrees of
freedom at the interface in question.

Definitions of the input parameters are as follows:

================= ==========================================================
**CONT_TANG_VEL** Name of the boundary condition (<bc_name>).
**SS**            Type of boundary condition (<bc_type>), where **SS**
                  denotes side set in the EXODUS II database.
<bc_id>           The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set
                  in EXODUS II) in the problem domain.
================= ==========================================================

------------
**Examples**
------------

The following is a sample card:
::

     BC = CONT_TANG_VEL SS 10

-------------------------
**Technical Discussion**
-------------------------

No Discussion. 




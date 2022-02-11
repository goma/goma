****************
**CURRENT_USER**
****************

::

	BC = CURRENT_USER SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/POTENTIAL)**

This boundary condition card is used to define a routine for a user-defined electrical
current density model. Definitions of the input parameters are as follows:

================ ===========================================================
**CURRENT_USER** Name of the boundary condition (<bc_name>).
**SS**           Type of boundary condition (<bc_type>), where **SS**
                 denotes side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain.
<float_list>     A list of float values separated by spaces which will be
                 passed to the user-defined subroutines so the user can
                 vary the parameters of the boundary condition. This list
                 of float values is passed as a one-dimensional double
                 array to the appropriate C function.
================ ===========================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = CURRENT_USER SS 100   10.0 3.14159

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




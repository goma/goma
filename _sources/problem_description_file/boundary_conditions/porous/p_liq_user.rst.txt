**************
**P_LIQ_USER**
**************

::

	BC = P_LIQ_USER SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(PCC/R_POR_LIQ_PRES)**

This boundary condition card is used to call a routine for a user-defined liquid-phase
pressure for porous flow problems at an external boundary of a material of one of the
following media types: *POROUS_SATURATED, POROUS_UNSATURATED,
POROUS_TWO_PHASE*.. Specification is made via the function
p_liq_user_surf in file “user_bc.c.” Definitions of the input parameters are as
follows:

============== ==============================================================
**P_LIQ_USER** Name of the boundary condition (<bc_name>).
**SS**         Type of boundary condition (<bc_type>), where **SS** denotes
               side set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set in
               EXODUS II) in the problem domain.
<float_list>   A list of float values separated by spaces which will be
               passed to the user-defined subroutine so the user can vary
               the parameters of the boundary condition. This list of float
               values is passed as a one-dimensional double array to the
               appropriate C function in file user_bc.c.
============== ==============================================================

------------
**Examples**
------------

The following is a sample input card with two parameters passed to function tuser:
::

   BC =P_LIQ_USER SS 100   273.13 100.0

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
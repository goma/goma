*********
**YUSER**
*********

::

	BC = YUSER SS <bc_id> <integer> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC/MASS)**

This is a user-defined mass concentration boundary. The user must supply the
relationship in function *yuser_surf* within *user_bc.c*.

============ ==============================================================
**YUSER**    Name of the boundary condition.
**SS**       Type of boundary condition (<bc_type>), where **SS** denotes
             node set in the EXODUS II database.
<bc_id>      The boundary flag identifier, an integer associated with
             <bc_type> that identifies the boundary location (node set in
             EXODUS II) in the problem domain.
<integer>    Species number.
<float_list> A list of float values separated by spaces which will be
             passed to the user-defined subroutines so the user can vary
             the parameters of the boundary condition. This list of float
             values is passed as a one-dimensional double array to the
             appropriate C function.
============ ==============================================================

------------
**Examples**
------------

The following sample input card applies a user flux condition to side set 100 for species 0 that requires two input parameters:
::

   BC = YUSER SS 100 0 .5 .5

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
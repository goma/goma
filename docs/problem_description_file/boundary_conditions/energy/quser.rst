*********
**QUSER** 
*********

::

	BC = QUSER SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition card is used to call a routine for a user-defined heat flux
model. Definitions of the input parameters are as follows:

============= ================================================================
**QUSER**     Name of the boundary condition (<bc_name>).
**SS**        Type of boundary condition (<bc_type>), where **SS** denotes
              side set in the EXODUS II database.
<bc_id>       The boundary flag identifier, an integer associated with
              <bc_type> that identifies the boundary location (side set in
              EXODUS II) in the problem domain.
<float_list>  A list of float values separated by spaces which will be
              passed to the user-defined subroutines so the user can vary
              the parameters of the boundary condition. This list of float
              values is passed as a one-dimensional double array to the
              quser_surf C function in file user_bc.c.
============= ================================================================

------------
**Examples**
------------

The following is a sample input card for a heat flux model requiring two parameters:
::

   BC = QUSER SS 100   10.0 3.14159

-------------------------
**Technical Discussion**
-------------------------

No Discussion.







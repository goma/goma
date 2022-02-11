**************
**YFLUX_USER**
**************

::

	BC = YFLUX_USER SS <bc_id> <integer> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary condition card is used to set mass flux to a user-prescribed function and
integrate by parts again. The user should provide detailed flux conditions in the
*mass_flux_user_surf* routine in *user_bc.c*. The flux quantity is specified on a per
mass basis so the heat and mass transfer coefficients are in units of L/t.

Definitions of the input parameters are as follows:

============== =================================================================
**YFLUX_USER** Name of the boundary condition (<bc_name>).
**SS**         Type of boundary condition (<bc_type>), where **SS** denotes
               side set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set in
               EXODUS II) in the problem domain.
<integer>      Species number of concentration.
<float_list>   A list of float values separated by spaces which will be
               passed to the user-defined subroutines so the user can vary
               the parameters of the boundary condition. This list of float
               values is passed as a one-dimensional double array to the
               appropriate C function.
============== =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = YFLUX_USER SS 2   0 .5 .5

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
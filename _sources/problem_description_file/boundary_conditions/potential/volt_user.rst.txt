*************
**VOLT_USER**
*************

::

	BC = VOLT_USER SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/POTENTIAL)**

This boundary condition card is used to specify a voltage or potential computed via a
user-defined function. Definitions of the input parameters are as follows:

============= =========================================================
**VOLT_USER** Name of the boundary condition (<bc_name>).
**SS**        Type of boundary condition (<bc_type>), where **SS**
              denotes side set in the EXODUS II database.
<bc_id>       The boundary flag identifier, an integer associated with
              <bc_type> that identifies the boundary location (side set
              in EXODUS II) in the problem domain.
<float_list>  A list of float values separated by spaces which will be
              passed to the user-defined subroutines so the user can
              vary the parameters of the boundary condition. This list
              of float values is passed as a one-dimensional double
              array to the appropriate C function.
============= =========================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = VOLT_USER SS 14 0.33 1000. 0.001 4e-5 1. 1. 353. 0.

-------------------------
**Technical Discussion**
-------------------------

In the **VOLT_USER** model currently implemented in GOMA, the electrolyte potential
is given by the linearized Butler-Volmer kinetic model as in the hydrogen-oxidation
reaction of a hydrogen-fueled polymer-electrolyte-membrane fuel cell. See the
*user_bc.c* routine for details.



--------------
**References**
--------------

No References.
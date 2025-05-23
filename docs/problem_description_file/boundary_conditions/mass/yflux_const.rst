***************
**YFLUX_CONST**
***************

::

	BC = YFLUX_CONST SS <bc_id> <integer> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary condition card is used to specify a constant diffusive mass flux of a
given species. This flux quantity can be specified on a per mass basis (e.g. with units of
g/ :math:`cm^2` /s) or on a per mole basis (e.g. with units of moles/ :math:`cm^2` /s), depending on the
user’s choice of units in the species unknown.

Definitions of the input parameters are as follows:

================ ===========================================================
**YFLUX_CONST**  Name of the boundary condition (<bc_name>).
**SS**           Type of boundary condition (<bc_type>), where SS
                 denotes side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain.
<integer>        Species number.
<float>          Value of diffusive mass flux; the units of this quantity
                 depends on the user’s choice of units for species
                 concentration.
================ ===========================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = YFLUX_CONST SS 1 0 10000.2

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
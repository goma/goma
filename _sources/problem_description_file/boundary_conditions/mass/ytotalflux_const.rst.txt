********************
**YTOTALFLUX_CONST**
********************

::

	BC = YTOTALFLUX_CONST SS <bc_id> <integer> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary condition card is used to specify a constant total mass flux (including
contribution from diffusion, migration, and convection) of a given species. This card
enables the treatment of the situation in which diffusion, migration and convection
fluxes cancel each other such that the total flux vanishes (e.g. is equal to zero). This
flux quantity can be specified on a per mass basis (i.e., with units of g/ 
:math:`cm^2`/s) or on a
per mole basis (e.g. with units of moles/:math:`cm^2`/s), depending on the user’s choice of units
in the species concentration unknown.

Definitions of the input parameters are as follows:

==================== =======================================================
**YTOTALFLUX_CONST** Name of the boundary condition (<bc_name>).
**SS**               Type of boundary condition (<bc_type>), where **SS**
                     denotes side set in the EXODUS II database.
<bc_id>              The boundary flag identifier, an integer associated with
                     <bc_type> that identifies the boundary location (side set
                     in EXODUS II) in the problem domain.
<integer>            Species number of concentration.
<float>              Value of total mass flux - the units of this quantity
                     depends on the user’s choice of units for species
                     concentration.
==================== =======================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = YTOTALFLUX_CONST SS 5   0   0.0

-------------------------
**Technical Discussion**
-------------------------

No Discussion.





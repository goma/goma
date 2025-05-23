************
**KIN_CHEM**
************

::

	BC = KIN_CHEM SS <bc_id> <float1> ... <floatn>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MESH)**

This boundary condition card is used to establish the sign of flux contributions to the
overall mass balance on boundaries so that movements are appropriately advancing or
receding depending on whether a species is a reactant or product in a surface reaction.

Definitions of the input parameters are as follows:

============ =======================================================================
**KIN_CHEM** Name of the boundary condition (<bc_name>).
**SS**       Type of boundary condition (<bc_type>), where **SS** denotes side set.
<bc_id>      The boundary flag identifier, an integer associated with
             <bc_type> that identifies the boundary location (set in
             EXODUS II) in the problem domain.
<float1>     Stoichiometric coefficient for species 0.
<floatn>     Stoichiometric coefficient for species *n* +1.
============ =======================================================================

The input function will read as many stoichiometric coefficients as specified by the
user for this card; the number of coefficients read is counted and saved. The
stoichiometric coefficient is +1 for products or -1 for reactants. When a species is a
product, the surface will advance corresponding to production/creation of mass of that
species, versus recession of that interface when a reaction leads to consumption of that
species.

------------
**Examples**
------------

Following is a sample card for two reactant and one product species:
::

     BC = KIN_CHEM SS 25   -1.0   -1.0   1.0

-------------------------
**Technical Discussion**
-------------------------

This function is built from the same function as boundary condition *KIN_LEAK*, i.e.,
*kin_bc_leak*, so the user is referred to discussions for this boundary condition for
appropriate details. The stoichiometric coefficients are read from the *KIN_CHEM* card
or set equal to 1.0 in the absence of *KIN_CHEM*.



--------------
**References**
--------------

No References.
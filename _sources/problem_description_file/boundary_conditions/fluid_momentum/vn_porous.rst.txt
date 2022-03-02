*************
**VN_POROUS**
*************

::

	BC = VN_POROUS SS <bc_id> <integer_list> <float>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MOMENTUM)**

This boundary condition is used to calculate the normal component of gas phase
velocity at the interface between a continuous gas phase and porous phase. The
condition is basically the unsaturated equivalent to *DARCY_CONTINUOUS*, and hence
is a condition on the normal component of fluid (gas) velocity on the continuous side of
the interface (see below). The flux on the porous medium side includes Darcy flux and
Fickian diffusive flux in the porous phase. The vapor flux into gas is used to determine
gas velocity. The condition is similar to the solid-liquid interface conditions that apply
to interfaces between a porous medium and an external gas (in which the energy
equation is used to solve for solvent concentration in the gas phase). This boundary
condition is still under development and has not been recently tested. Its last use was
for evaporation from a porous unsaturated film in a sol-gel application (see references
below).

There are three values to be supplied for the <integer_list>; definitions of the input
parameters are as follows:

============= ==================================================================
**VN_POROUS** Name of the boundary condition (<bc_name>).
**SS**        Type of boundary condition (<bc_type>), where **SS** denotes
              side set in the EXODUS II database.
<bc_id>       The boundary flag identifier, an integer associated with
              <bc_type> that identifies the boundary location (side set in
              EXODUS II) in the problem domain.
<integer1>    Element block ID of solid, porous phase from the EXODUS
              II database.
<integer2>    Element block ID of gas phase from the EXODUS II
              database.
<integer3>    **Set to zero for now**, indicating that this condition pertains
              to the primary liquid solvent vapor. At some point this will
              be generalized to include all vapor components.
<float>       Density of pure solvent vapor.
============= ==================================================================

------------
**Examples**
------------

The following is a sample input card:
::

    BC = VN_POROUS SS 5   1 2 0 1.e-3

This condition applies to internal side set 5, which defines the interface between
element block 1 (the solid porous phase which has *Media Type* of
*POROUS_PART_SAT* or *POROUS_TWO_PHASE*) and element block 2 (the fluid
phase which has Media Type *CONTINUOUS*). It is based on the flux of liquid solvent
in the porous phase (denoted by the integer 0), the vapor form of which has a density of
1.e-3. The condition results in a blowing or sucking velocity at the interface in the fluid
(gas) continuous phase.

-------------------------
**Technical Discussion**
-------------------------

The functional form of this boundary condition is

.. figure:: /figures/116_goma_physics.png
	:align: center
	:width: 90%

Here, the left hand side is the total flux of liquid solvent, in both gas and liquid phases.
The first two terms are the Darcy pressure driven contributions, and the second two
terms are the Fickian flux contributions.

This condition would be useful for predicting the gas-flow pattern above a drying
porous matrix, in which the vapor flux being driven out of the porous skeleton were a
mass source to drive flow in the surrounding gas. The condition has not been tested
since 1995.



--------------
**References**
--------------

GTM-028.0: Modeling Drying of Dip-Coated Films with Strongly-Coupled Gas Phase
Natural Convection, R. A. Cairncross, 1999.

.. TODO - Line 71 have photos that needs to be replaced with the real equation.
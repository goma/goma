*********************
**KINEMATIC_SPECIES**
*********************

::

	BC = KINEMATIC_SPECIES SS <bc_id> <integer>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary condition card is used to impose an interphase species flux continuity
constraint on species components undergoing phase change between two materials.
The species conservation equation (see *EQ* card and *species_bulk*) for a single gas or
liquid phase component requires two boundary conditions because of the multivalued,
discontinuous concentration at the interface. This condition should be used in
conjunction with *VL_EQUIL* tie condition for each species. Definitions of the input
parameters are as follows:

===================== ==========================================================
**KINEMATIC_SPECIES** Name of the boundary condition (<bc_name>).
**SS**                Type of boundary condition (<bc_type>), where **SS**
                      denotes side set in the EXODUS II database.
<bc_id>               The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location (side set
                      in EXODUS II) in the problem domain.
<integer>             Species number.
<float1>              Unused floating point number.
===================== ==========================================================

This boundary condition is typically applied to multicomponent two-phase flows that
have rapid mass exchange between phases, rapid enough to induce a diffusion velocity
at the interface, and to thermal contact resistance type problems. The best example of
this is rapid evaporation of a liquid component into a gas.

------------
**Examples**
------------

Following is a sample card:
::

     BC = KINEMATIC_SPECIES SS 10 2 0.0

This card invokes the species flux balance condition on species 2 at shared side set 10
to be applied to the liquid phase convective diffusion equation. It should be used in
conjunction with a *VL_EQUIL* type condition on the same species, but from the
bounding phase. Note: side set 10 must be a double-sided side set between two
materials (i.e., must be attached to both materials), each deploying basis function
interpolation of type **Q1_D** or **Q2_D**.

-------------------------
**Technical Discussion**
-------------------------

The condition only applies to interphase mass transfer problems with discontinuous (or
multivalued) variables at an interface, and it must be invoked on fields that employ the
**Q1_D** or **Q2_D** interpolation functions to “tie” together or constrain the extra degrees
of freedom at the interface in question. The mathematical form is

.. math::

   \underline{n} \cdot \left[ \left( \underline{v}_s - {\underline{v}}^l \right) {y_i}^l \rho^l - {\underline{j}_i}^l \right] = \underline{n} \cdot \left[ \left( \underline{v}_s - {\underline{v}}^g \right) {y_i}^g \rho^g - {\underline{j}_i}^g \right]

Here :math:`{\underline{v}}^l` and :math:`{\underline{v}}^g` are the gas and liquid velocity vectors at the free surface, respectively; :math:`\underline{v}_s` is the mesh velocity at the same location; :math:`\rho^l` and :math:`\rho^g` are the liquid and gas phase
densities, respectively; :math:`{y_i}^l` and :math:`{y_i}^g` are the liquid and gas phase volume fractions of
component :math:`i`; and :math:`{\underline{j_i}}^l` and :math:`{\underline{j_i}}^g` the mass fluxes of component :math:`i`. This condition constrains
only one of two phase concentrations at the discontinuous interface. The other needs to
come from a Dirichlet boundary condition like (BC =) Y, or an equilibrium boundary
condition like *VL_EQUIL*.



--------------
**References**
--------------

Schunk, P. R. and Rao, R. R. 1994. “Finite element analysis of multicomponent twophase
flows with interphase mass and momentum transport”, Int. J. Numer. Meth.
Fluids, 18, 821-842.

GTM-007.1: New Multicomponent Vapor-Liquid Equilibrium Capabilities in GOMA,
December 10, 1998, A. C. Sun

..
	 TODO - The picture in line 64 needs to be exchanged with the equation. In lines 68-73, where it says "**EQUATION**" everytime there is supposed to be something from the equation that needs to be written. 
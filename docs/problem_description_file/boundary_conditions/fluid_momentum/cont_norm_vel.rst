*****************
**CONT_NORM_VEL**
*****************

::

	BC = CONT_NORM_VEL SS <bc_id>

-----------------------
**Description / Usage**
-----------------------

**(SIC/MOMENTUM)**

This boundary condition card is similar to the *VELO_NORM_DISC* card except that it
enforces a continuous normal velocity component in a discontinuous boundary field.
The condition only applies to interphase mass, heat, and momentum transfer problems
with discontinuous (or multivalued) variables at an interface, and it must be invoked on
fields that employ the **Q1_D** or **Q2_D** interpolation functions to “tie” together or
constrain the extra degrees of freedom at the interface in question.

Definitions of the input parameters are as follows:

================= ============================================================
**CONT_NORM_VEL** Name of the boundary condition (<bc_name>).
**SS**            Type of boundary condition (<bc_type>), where **SS**
                  denotes side set in the EXODUS II database.
<bc_id>           The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set
                  in EXODUS II) in the problem domain.
================= ============================================================

This boundary condition is typically applied to multicomponent two-phase flows that
have rapid mass exchange between phases, rapid enough to induce a diffusion velocity
at the interface, and to thermal contact resistance type problems. The best example of
this is rapid evaporation of a liquid component into a gas.

------------
**Examples**
------------

The following is a sample card:
::

     BC = CONT_NORM_VEL SS 10

-------------------------
**Technical Discussion**
-------------------------

No Discussion. 



--------------
**References**
--------------

No References.
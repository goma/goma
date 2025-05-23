********************
**VELO_NORMAL_DISC**
********************

::

	BC = VELO_NORMAL_DISC SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MOMENTUM)**

This boundary condition card balances mass loss from one phase to the gain from an
adjacent phase. It is the same as the *KINEMATIC_DISC* card but is applied to the 
fluid
momentum equation. The condition only applies to interphase mass, heat, and
momentum transfer problems with discontinuous (or multivalued) variables at an
interface, and it must be invoked on fields that employ the **Q1_D** or **Q2_D**
interpolation functions to “tie” together or constrain the extra degrees of freedom 
at the
interface in question.

Definitions of the input parameters are as follows:

===================== ==========================================================
**VELO_NORMAL_DISC**  Name of the boundary condition.
**SS**                Type of boundary condition (<bc_type>), where **SS**
                      denotes side set in the EXODUS II database.
<bc_id>               The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location (side set
                      in EXODUS II) in the problem domain. It is important
                      to note that this side set should be shared by both
                      element blocks for internal boundaries.
<float>               Set to zero for internal interfaces; otherwise used to
                      specify the mass average velocity across the interface
                      for external boundaries.
===================== ==========================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = VELO_NORMAL_DISC SS 66 0.0

is used at internal side set 10 (note, it is important that this side set include 
elements
from both abutting materials) to enforce the overall conservation of mass exchange.

-------------------------
**Technical Discussion**
-------------------------

* This boundary condition card applies the following constraint to nodes on the side
  set: 

.. math::

  \rho_1 n \cdot \left(v - v_s\right) |_1 = \rho_2 n \cdot \left(v - v_s\right) |_2

  

where 1 denotes evaluation in phase 1 and 2 denotes evaluation in phase 2. This
constraint replaces only one of the momentum equations present at an internal
discontinuous boundary between materials. There usually must be another
momentum boundary condition applied to this side set. In addition, there must also
be a distinguishing condition applied to the mesh equations if mesh motion is part
of the problem.

* This boundary condition is typically applied to multicomponent two-phase flows
  that have rapid mass exchange between phases, rapid enough to induce a diffusion
  velocity at the interface, and to thermal contact resistance type problems. The best
  example of this is rapid evaporation of a liquid component into a gas.



--------------
**References**
--------------

No References.
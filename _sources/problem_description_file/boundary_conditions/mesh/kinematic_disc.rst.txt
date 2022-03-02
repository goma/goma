******************
**KINEMATIC_DISC**
******************

::

	BC = KINEMATIC_DISC SS <bc_id> <float1>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MESH)**

This boundary condition card is used as a distinguishing condition on the mesh motion
equations (viz. *mesh1, mesh2*, and *mesh3* under the *EQ* card) in the special case of an
interface between two fluids of different density (e.g. a gas and a liquid, both meshed
up as *Goma* materials) through which a phase transition is occurring and there is a
discontinuous velocity (see the mathematical form in the technical discussion below).
Like the *KINEMATIC* boundary condition, it is used to distinguish a material surface
between two phases exchanging mass. In two dimensions, this condition is
automatically applied to the normal component of the vector mesh equations which is
rotated into normal-tangential form. In three dimensions, the application of this
boundary condition needs to be further directed with the *ROT* cards (see *Rotation
Specifications*). The application of this condition should be compared with
*KINEMATIC_PETROV* and *KINEMATIC_COLLOC*.

This condition must be applied to problem description regions using the Q1_D or
Q2_D interpolation type, indicating a discontinuous variable treatment at the interface
(see *EQ* card).

Definitions of the input parameters are as follows:

=================== ===========================================================
**KINEMATIC_DISC**  Name of the boundary condition (<bc_name>).
**SS**              Type of boundary condition (<bc_type>), where **SS**
                    denotes side set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (side set
                    in EXODUS II) in the problem domain.
<float1>            Set to zero for internal interfaces; otherwise used to
                    specify the mass average velocity across the interface
                    for external boundaries.
=================== ===========================================================

------------
**Examples**
------------

The following sample card
::

     BC = KINEMATIC_DISC SS 10 0.0

is used at internal side set 10 (note, it is important that this side set include elements
from both abutting materials) to enforce the overall conservation of mass exchange.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is typically applied to multicomponent two-phase flows that
have rapid mass exchange between phases, rapid enough to induce a diffusion velocity
at the interface. The best example of this is rapid evaporation of a liquid component
into a gas.

This boundary condition card is used for a distinguishing condition and its functional
form is:

.. math::

   \rho_1 \underline{n} \cdot \left( \underline{v} - \underline{v}_s \right) \rvert_1 = \rho_2 \underline{n} \cdot \left( \underline{v} - \underline{v}_s \right) \rvert_2

where 1 denotes evaluation in phase 1 and 2 denotes evaluation in phase 2.

This condition is applied to the rotated form of the mesh equations. The condition only
applies to interphase mass, heat, and momentum transfer problems with discontinuous
(or multivalued) variables at an interface, and it must be invoked on fields that employ
the **Q1_D** or **Q2_D** interpolation functions to “tie” together or constrain the extra
degrees of freedom at the interface in question (see for example boundary condition
*VL_EQUIL_PSEUDORXN*).



--------------
**References**
--------------

GTM-015.1: Implementation Plan for Upgrading Boundary Conditions at
Discontinuous-Variable Interfaces, January 8, 2001, H. K. Moffat
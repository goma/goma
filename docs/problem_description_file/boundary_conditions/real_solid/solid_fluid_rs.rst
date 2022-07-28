******************
**SOLID_FLUID_RS**
******************

::

	BC = SOLID_FLUID_RS SS <bc_id> <integer1> <integer2> [float]

-----------------------
**Description / Usage**
-----------------------

**(PCC/VECTOR REALSOLID)**

Used for fluid-structure interaction problems, the *SOLID_FLUID_RS* condition
equates the normal traction between adjacent fluid and solid materials. (By “normal
traction” we mean the tangential and normal force components, per unit area.) This
condition is only to be used on boundaries between regions of *ARBITRARY* mesh
motion with fluid-momentum equations and of *TOTAL_ALE* mesh motion (cf.
*SOLID_FLUID* boundary condition card for *LAGRANGIAN* mesh motion regions),
with solid momentum equations (or mesh equations) - see *Mesh Motion* and *EQ* cards.
All elements on both sides of the interface must have the same element type (the same
order of interpolation and basis functions) e.g., Q1 or Q2. Also, such interfaces must
include element sides from both sides of the interface in the defining side set.

Definitions of the input parameters are as follows:

================== ============================================================
**SOLID_FLUID_RS** Name of the boundary condition (<bc_name>).
**SS**             Type of boundary condition (<bc_type>), where **SS**
                   denotes side set in the EXODUS II database.
<bc_id>            The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (side set
                   in EXODUS II) in the problem domain.
<integer1>         Element block ID of the solid phase (of TOTAL_ALE
                   motion type) from the EXODUS II database.
<integer2>         Element block ID of the liquid phase from the
                   EXODUS II database.
[float]            Scale factor for stress balance for nondimensionalization.
                   This parameter, a multiplier on the
                   liquid phase contribution, is optional; the default is 1.0.
================== ============================================================

------------
**Examples**
------------

The following set of input cards is a sample specification for a fluid-structure
interaction problem:
::

     BC = SOLID_FLUID_RS SS 5   2   1

::

     BC = NO_SLIP_RS SS 5   2   1

::

     BC = KIN_DISPLACEMENT SS   5   2

In this example, side set 5 is a boundary between a solid rubber blade and a liquid; the
material in element block 2 is the blade, and the material in element block 1 is the fluid.
Along the blade, a companion boundary condition is applied to ensure no slip on the
same side set. Also, because this condition involves a *TOTAL_ALE* mesh region, a
*KIN_DISPLACEMENT* boundary condition is needed on the same side set to force the
solid boundary to follow the side set.

-------------------------
**Technical Discussion**
-------------------------

The functional form of the *SOLID_FLUID_RS* boundary condition is:

.. math::

     \lambda (\underline n \cdot \underline T) = \underline n \cdot \underline \sigma

     

where :math:`\underline{T}` is the fluid phase stress tensor given by any one of the specified fluid-phase
constitutive equations, and :math:`\underline{\sigma}` is the real-solid solid phase stress tensor, also given by
any one of the solid-phase constitutive equation (see *Mat* file specifications). :math:`\lambda` is a
scaling factor that defaults to unity (and is usually best taken as such unless some
scaling is invoked). With this boundary condition, the local residual and Jacobian
contributions from the fluid mechanics momentum equations (on the *ARBITRARY* side
of the boundary) are added into the weak form of the residual and Jacobian entries for
the real-solid solid mechanics equations (viz. the *EQ = mom_solid* options on the
real-solid *TOTAL_ALE* side of the boundary).

*TOTAL_ALE* mesh motion regions cannot be porous and deformable (as of 11/19/
2001).



--------------
**References**
--------------

GT-005.3: THE NEW TOTAL-ARBITRARY-LAGRANGIAN-EULERIAN (TALE)
CAPABILITY and its applicability to coating with/on deformable media, August 6,
1999, P. R. Schunk

SAND2000-0807: TALE: An Arbitrary Lagrangian-Eulerian Approach to Fluid-
Structure Interaction Problems, P. R. Schunk (May 2000)

.. 
	TODO - Line 72 contains a picture that needs to be repalced with the equation. 
***************
**Mesh Motion**
***************

::

	Mesh Motion = {char_string}

-----------------------
**Description / Usage**
-----------------------

This card is required for each material section in the *Problem Description File* even if a moving mesh problem is not being solved. It is used to specify the method which prescribes the movement of nodes within the mesh. Valid options for {char_string} are:

====================== ========================================================
**ARBITRARY**          This option triggers the implicit pseudo-solid
                       domainmapping technique using the constitutive equation
                       designated in the corresponding file.mat (see *Material
                       File* description); with this technique, the boundaries of the domain are controlled by distinguishing conditions coupled with the problem physics, and the interior nodes move independently of the problem physics.
**LAGRANGIAN**         This option triggers coupling the motion of nodes on the
                       interior of the domain to the deformation of an elastic solid. If the solid is incompressible, this technique uses a pressure (Lagrange multiplier) to couple the solid deformation and the local solvent concentration.
**DYNAMIC_LAGRANGIAN** This option triggers coupling the motion of nodes on the
                       interior of the domain to the deformation of an elastic solid, *including solid inertia*. If the solid is incompressible, this technique uses a pressure (Lagrange multiplier) to couple the solid deformation and the local solvent concentration. Together with the equation term multiplier on the mass matrix (see *EQ* card) and a “transient” specification on the *Time Integration* Card, this option will invoke a Newmark- Beta time integration scheme for the inertia term in the *R_MESH* equations.
**TOTAL_ALE**          This option allows motion of nodes on the interior of the
                       domain of a solid region to be independent of the material motion. TALE is an acronym for “Total Arbitrary
                       Lagrangian Eulerian” mesh motion. This is typically used in elastic solids in which large scale deformation makes motions under the LAGRANGIAN option unmanageable. If the solid is incompressible, this technique uses a pressure (Lagrange multiplier) to couple the solid deformation and the local solvent concentration. Invoking this option requires mesh equations and real solid equations, as described on the EQ card. Other relevant cards that are often used with this option are *KINEMATIC_DISPLACEMENT* boundary condition, *DX_RS, DY_RS, DZ_RS* boundary conditions, *FORCE_RS, FLUID_SOLID_RS*, and others. See references for more detailed usage procedures.
====================== ========================================================

------------
**Examples**
------------

The following is a sample card that sets the mesh motion scheme to be arbitrary:
::

   Mesh Motion = ARBITRARY

-------------------------
**Technical Discussion**
-------------------------

For the *TOTAL_ALE* mesh motion option we must supply elastic properties and solid constitutive equations for both the mesh and the real solid. It is best to consult the example tutorials cited below for details.



--------------
**References**
--------------

GT-005.3: THE NEW TOTAL-ARBITRARY-LAGRANGIAN-EULERIAN (TALE)
CAPABILITY and its applicability to coating with/on deformable media, August 6,
1999, P. R. Schunk

GT-006.3: Slot and Roll coating with remeshing templates and tutorial for GOMA and CUBIT/MAPVAR, August 3, 1999, R. R. Lober and P. R. Schunk
*****************
Polymer Viscosity
*****************

::

   Pseudo-Solid Constitutive Equation = {model_name}

-----------------------
**Description / Usage**
-----------------------

This card specifies the constitutive equation used to control mesh motion for arbitrary
Lagrangian Eulerian solid mechanics and is required for use with the *TOTAL_ALE*
mesh motion type (see *Mesh Motion* card). Details are discussed in references provided
below.

The single input parameter is the type of model for the constitutive equation:

+-----------------+------------------------------------------------------------------------------------------+
|{model_name}     |The name of the constitutive equation; {model_name} can be one of the following:          |
|                 |                                                                                          |
|                 | * **LINEAR** - the mesh deformations are assumed to be small and thus simplifies the     |
|                 |   analysis of strain and stress.                                                         |
|                 | * **NONLINEAR** - a nonlinear neo-Hookean elastic model for which the deformations can be|
|                 |   large without loss of frame invariance. This is the recommended model (and all         |
|                 |   materials currently default to NONLINEAR if the mesh is arbitrary).                    |
|                 |                                                                                          |
|                 |The following models are allowed but not recommended.                                     |
|                 |                                                                                          |
|                 | * **HOOKEAN_PSTRAIN** - a nonlinear neo-Hookean model with plane strain assumption       |
|                 |   2D only).                                                                              |
|                 | * **INCOMP_PSTRAIN** - an incompressible nonlinear neo-Hookean model with plane strain   |
|                 |   and a Lagrangian pressure constraint on the volume.                                    |
|                 | * **INCOMP_3D** - Incompressible version of the neo-Hookean solid in a special Segalman  |
|                 |   formulation that removes the volume-change from the strain tensor (like the            |
|                 |   *INCOMP_PSTRAIN* model above), and is specifically designed for 3D applications (not a |
|                 |   widely used option).                                                                   |
+-----------------+------------------------------------------------------------------------------------------+

Note again the requirement that the *Mesh Motion* type for the material in which this
constitutive equation applies must be *TOTAL_ALE*.

------------
**Examples**
------------

::

   Pseudo-Solid Constitutive Equation = NONLINEAR

This card specifies the mesh motion in the ALE solid region is to conform to the
nonlinear elastic model, as described on the *Solid Constitutive Equation* card. This card
is required together with *Pseudo-Solid Lame Mu and Pseudo-Solid Lame Lambda* cards.

-------------------------
**Technical Discussion**
-------------------------

The Pseudo-Solid mesh motion, like the *ARBITRARY* mesh motion, is governed by the
equations of elasticity. These cards, together with the other cards required by the real
solid constitutive behavior, are required for ALE solid mechanics. The theory is
explained in detail in the provided references. Throughout the boundary condition
options, the user will notice an appended *_RS*. This signifies that the boundary
conditions apply to the real-solid elasticity in *TOTAL_ALE* problems. All other
boundary conditions on force and displacement, viz. those without the *_RS*, are applied
to the mesh motion.



--------------
**References**
--------------

GT-005.3: THE NEW TOTAL-ARBITRARY-LAGRANGIAN-EULERIAN (TALE)
CAPABILITY and its applicability to coating with/on deformable media, August 6,
1999, P. R. Schunk

SAND2000-0807: TALE: An Arbitrary Lagrangian-Eulerian Approach to Fluid-
Structure Interaction Problems, P. R. Schunk (May 2000)

Sackinger, P. A., Schunk, P. R. and Rao, R. R. 1995. "A Newton-Raphson Pseudo-Solid
Domain Mapping Technique for Free and Moving Boundary Problems: A Finite
Element Implementation", J. Comp. Phys., 125 (1996) 83-103.


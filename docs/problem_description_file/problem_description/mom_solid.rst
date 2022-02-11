*************
**mom_solid**
*************

::

	EQ = mom_solid{1|2|3} {Galerkin_wt} {D1|D2|D3}_RS {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for one component of
a solid momentum equation. Definitions of the input parameters are defined below.
Note that <floatlist> contains five constants for the Solid Momentum equation defining
the constant multipliers for each type of term in the equation. The Galerkin weight and
the interpolation function must be the same for the code to work properly.

+----------------------------------------+---------------------------------------------------+
|**mom_solid1 | mom_solid2 | mom_solid3**|Name of the equation to be solved, where the 1, 2  |
|                                        |and 3 components correspond to one of the principal|
|                                        |coordinate directions, e.g. X, Y and Z for         |
|                                        |Cartesian geometry.                                |
+----------------------------------------+---------------------------------------------------+
|{Galerkin_wt}                           |Two-character value that defines the type of       |
|                                        |weighting function for this equation, where:       |
|                                        |                                                   |
|                                        | * **Q1**-Linear                                   |
|                                        | * **Q2**-Quadratic                                |
+----------------------------------------+---------------------------------------------------+
|**D1_RS | D2_RS | D3_RS**               |Name of the variable associated with the 1, 2 or 3 |
|                                        |principal coordinate direction for this component  |
|                                        |equation.                                          |
+----------------------------------------+---------------------------------------------------+
|{Interpol_fnc}                          |Two-character value that defines the interpolation |
|                                        |function used to represent the variable            |
|                                        |**D1_RS, D2_RS** or **D3_RS** where:               |
|                                        |                                                   |
|                                        | * **Q1**-Linear                                   |
|                                        | * **Q2**-Quadratic                                |
+----------------------------------------+---------------------------------------------------+
|<float1>                                |Multiplier on mass matrix term ( d ⁄dt ).          |
+----------------------------------------+---------------------------------------------------+
|<float2>                                |Multiplier on advective term.                      |
+----------------------------------------+---------------------------------------------------+
|<float3>                                |Multiplier on boundary term                        |
|                                        |( :math:`\underline{n}` • flux  ).                 |
+----------------------------------------+---------------------------------------------------+
|<float4>                                |Multiplier on diffusion term.                      |
+----------------------------------------+---------------------------------------------------+
|<float5>                                |Multiplier on source term.                         |
+----------------------------------------+---------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card for the first solid mesh equation that uses linear
continuous interpolation and turns on all term multipliers except for the mass matrix:
::

   EQ = mom_solid1 Q1 D1_RS Q1 0. 1. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

The solid momentum equations are used as a second set of displacement equations
when the ALE (arbitrary-Lagrangian Eulerian) technique is used in the solid phase as
well as the liquid phase. We have termed this capability TALE for “total arbitraryd
Lagrangian Eulerian” and details of implementation and usage for Goma can be found
in Schunk (2000).



--------------
**References**
--------------

SAND2000-0807: TALE: An Arbitrary Lagrangian-Eulerian Approach to Fluid-
Structure Interaction Problems, P. R. Schunk (May 2000)
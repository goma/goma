*****************************************
**lagr_mult_1, lagr_mult_2, lagr_mult_3**
*****************************************

::

	EQ = lagr_mult_{1|2|3} {Galerkin_wt} LM{1|2|3} {Interpol_fnc}

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a Langrange multiplier vector equation for
imposition of the kinematic boundary condition at a fluid/solid interface. It is used
soley for the overset grid capability in Goma (cf. GT-026.2). Definitions of the input
parameters are defined below. The Galerkin weight and the interpolation function must
be the same for the code to work properly.

+-------------------------------------------+-----------------------------------------+
|**lagr_mult_1 | lagr_mult_2 | lagr_mult_3**|Name of the equation to be solved. The   |
|                                           |appended number indexes with the         |
|                                           |dimension of the problem, viz.           |
|                                           |lagr_mult_1 and lagr_mult_2 equations    |
|                                           |are required for a                       |
|                                           |two dimensional problem.                 |
+-------------------------------------------+-----------------------------------------+
|{Galerkin_wt}                              |Two-character value that defines the type|
|                                           |of weighting                             |
|                                           |function for this equation, where:       |
|                                           |                                         |
|                                           | * **P0**-Constant Discontinuous         |
|                                           | * **P1**-Linear Discontinuous           |
|                                           | * **Q1**-Linear Continuous              |
|                                           | * **Q2**-Quadratic Continuous           |
+-------------------------------------------+-----------------------------------------+
|**LM1 | LM2 | LM3**                        |Name of the variable associated with     |
|                                           |this equation.                           |
+-------------------------------------------+-----------------------------------------+
|{Interpol_fnc}                             |Two-character value that defines the     |
|                                           |interpolation function used to represent |
|                                           |the variable **P**, where:               |
|                                           |                                         |
|                                           | * **P0**-Constant Discontinuous         |
|                                           | * **P1**-Linear Discontinuous           |
|                                           | * **Q1**-Linear Continuous              |
|                                           | * **Q2**-Quadratic Continuous           |
+-------------------------------------------+-----------------------------------------+

Basically when the level-set field (actually phase field 1, cf. F1 equation) that
corresponds to solid/fluid boundary defined by an overset grid (using the Slave Level
Set Card) intersects an element, the equations associated with that element will get the
kinematic boundary condition for the fluid-structure interaction, which basically
equates the fluid velocity to the solid velocity. In elements that don’t contain the solid/
fluid boundary, the equations are trivialized so that they are condensed out of the
system to be solved.

------------
**Examples**
------------

The following is a sample cards are required for the overset grid capability for twodimensional problems. It is recommended that P0 (element constant) interpolation
functions be used.
::

   EQ = lagr_mult_1 P0 LM1 P0 1

::

   EQ = lagr_mult_2 P0 LM2 P0 1

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

GT-026.3 GOMA’s Overset Mesh Method: User Tutorial. November 19, 2003. P. R.
Schunk and E. D. Wilkes
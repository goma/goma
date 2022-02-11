***************
**shell_angle**
***************

::

	EQ =shell_angle{1 | 2} {Galerkin_wt} {SH_ANG1 | SH_ANG2} {Interpol_fnc}

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a definition equation for the surface
orientation angle in a 2-dimensional bar element. It applies only to shell element
blocks. Note that this equation is available in three-dimensional problems but is in fact
untested at this time.. The shell angle equation(s) determine the components of the
normal vector to the shell surface; since its magnitude is 1 by definition, one less
degree of freedom is required tha the number of coordinates. The Galerkin weight and
the interpolation function must be the same for the code to work properly.

+--------------------+----------------------------------------------------------+
|**shell_angle{1|2}**|Name of the equation to be solved.                        |
+--------------------+----------------------------------------------------------+
|{Galerkin_wt}       |Two- or four-character value that defines the type of     |
|                    |weighting function for this equation, where:              |
|                    |                                                          |
|                    | * **Q1**-Linear                                          |
|                    | * **Q2**-Quadratic                                       |
+--------------------+----------------------------------------------------------+
|**SH_ANG{1|2}**     |SH_ANG{1|2} Name of the variable associated with the shell|
|                    |angle equation.                                           |
+--------------------+----------------------------------------------------------+
|{Interpol_fnc}      |Two- or four-character value that defines the             |
|                    |interpolation function used to represent the variable     |
|                    |**SH_ANG** where:                                         |
|                    |                                                          |
|                    | * **Q1**-Linear Continuous                               |
|                    | * **Q2**-Quadratic Continuous                            |
+--------------------+----------------------------------------------------------+

This equation requires no equation term multiplier entries.

------------
**Examples**
------------

The following are sample cards that use linear continuous curvature interpolation and
weight function:
::

   EQ = shell_angle1 Q1 SH_ANG1 Q1

::

   EQ = shell_angle2 Q1 SH_ANG2 Q2

The second card applies only to 3D problems.

-------------------------
**Technical Discussion**
-------------------------

For 2D problems, the defining equation is: :math:`\Theta` = atan[ :math:`n_x`, :math:`n_y`]  where Q is shell_angle1 and :math:`n_x` and :math:`n_y` are the components of the normal vector to the shell surface. There is an analogous definition for shell_angle2.



--------------
**References**
--------------

No References.
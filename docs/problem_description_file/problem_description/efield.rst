**********
**efield**
**********

::

	EQ = efield{1 | 2 | 3} {Galerkin_wt} {E1 | E2 | E3} {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a definition equation for the vector electric
field, which is the gradient of the voltage or potential field (see voltage equation).
Hence, these equations (two components in two dimensions, and three components in
three dimensions) must be solved together with the voltage equation.

+-------------------------------+--------------------------------------------------------+
|**efield1 | efield2 | efield3**|Name of the equation to be solved.                      |
+-------------------------------+--------------------------------------------------------+
|{Galerkin_wt}                  |Two-character value that defines the type of weighting  |
|                               |function for this equation, where:                      |
|                               |                                                        |
|                               | * **Q1**-Linear                                        |
|                               | * **Q2**-Quadratic                                     |
+-------------------------------+--------------------------------------------------------+
|**E1 | E2 | E3**               |Name of the variable associated with this equation.     |
+-------------------------------+--------------------------------------------------------+
|{Interpol_fnc}                 |Two-character value that defines the interpolation      |
|                               |function used to represent the variable                 |
|                               |**E1, E2,** or **E3,** where:                           |
|                               |                                                        |
|                               | * **Q1**-Linear                                        |
|                               | * **Q2**-Quadratic                                     |
+-------------------------------+--------------------------------------------------------+
|<float1>                       |Multiplier on advective term.                           |
+-------------------------------+--------------------------------------------------------+
|<float2>                       |Multiplier on source term.                              |
+-------------------------------+--------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is an example of efield-equation specification using linear elements in
two dimensions. Notice the companion voltage equation.
::

   EQ = efield1 Q1 E1 Q1 1. 1.

::

   EQ = efield2 Q1 E1 Q1 1. 1.

::

   EQ = voltage Q1 V Q1 1. 1. 1. 1. 1. 1.

This set of equations is required for applying an electrohydrodynamic force to the fluid
momentum equations (see *Navier-Stokes Source* card. )

-------------------------
**Technical Discussion**
-------------------------

The electric field is defined by :math:`\underline{E}` = –:math:`\Delta` :math:`\phi`. In some cases it may be more convenient to
solve equations for the potential field and the electric field.




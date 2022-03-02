*****************
**shell_tension**
*****************

::

	EQ = shell_tension {Galerkin_wt} {TENS} {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving an equation for tension using the structural
shell capability in Goma. The capability is based on inextensible cylindrical shells.
One material property is associated with this equation and that is the bending stiffness.
Note that <floatlist> contains one constant and it should always be set to one. The
Galerkin weight and the interpolation function must be the same for the code to work
properly.

+-----------------+----------------------------------------------------------+
|**shell_tens**   |Name of the equation to be solved.                        |
+-----------------+----------------------------------------------------------+
|{Galerkin_wt}    |Two- or four-character value that defines the type of     |
|                 |weighting function for this equation, where:              |
|                 |                                                          |
|                 | * **Q1**-Linear                                          |
|                 | * **Q2**-Quadratic                                       |
+-----------------+----------------------------------------------------------+
|**TENS**         |Name of the variable associated with the shell tension    |
|                 |equation.                                                 |
+-----------------+----------------------------------------------------------+
|{Interpol_fnc}   |Two- or four-character value that defines the             |
|                 |interpolation function used to represent the variable     |
|                 |**TENS** where:                                           |
|                 |                                                          |
|                 | * **Q1**-Linear Continuous                               |
|                 | * **Q2**-Quadratic Continuous                            |
+-----------------+----------------------------------------------------------+
|<float1>         |Multiplier on whole equation. Set to 1.0.                 |
+-----------------+----------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that uses linear continuous tension interpolation and
weight function:
::

   EQ = momentum1 Q1 TENS Q1 1.0

-------------------------
**Technical Discussion**
-------------------------

Complete tutorial on the use of this equation exists. See GT-27.1.

----------
**Theory**
----------

The structural shell equation capability in Goma builds on the shell-element capability
built by Pat Notz and Ed Wilkes in FY03. Basically we are solving the following
equations for the shell tension (this card) and shell curvature (see shell_curvature
equation):

.. figure:: /figures/303_goma_physics.png
	:align: center
	:width: 90%


--------------
**References**
--------------

GT-27

..
	TODO - Line 67 contains a photo that needs to be written as an equation.
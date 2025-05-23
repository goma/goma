*******************
**shell_curvature**
*******************

::

	EQ = shell_curvature {Galerkin_wt} {K} {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a definition equation for total curvature of a
two-dimensional shell element. Note that this equation is not yet available in three
dimensions. The curvature is required by the inextensible cylindrical shell capability in
Goma. See references cited below. Note that <floatlist> contains one constant and it
should always be set to one. The Galerkin weight and the interpolation function must
be the same for the code to work properly.

+-------------------+----------------------------------------------------------+
|**shell_curvature**|Name of the equation to be solved.                        |
+-------------------+----------------------------------------------------------+
|{Galerkin_wt}      |Two- or four-character value that defines the type of     |
|                   |weighting function for this equation, where:              |
|                   |                                                          |
|                   | * **Q1**-Linear                                          |
|                   | * **Q2**-Quadratic                                       |
+-------------------+----------------------------------------------------------+
|**K**              |Name of the variable associated with the shell curvature  |
|                   |equation.                                                 |
+-------------------+----------------------------------------------------------+
|{Interpol_fnc}     |Two- or four-character value that defines the             |
|                   |interpolation function used to represent the variable     |
|                   |**K** where:                                              |
|                   |                                                          |
|                   | * **Q1**-Linear Continuous                               |
|                   | * **Q2**-Quadratic Continuous                            |
+-------------------+----------------------------------------------------------+
|<float1>           |Multiplier on whole equation. Set to 1.0.                 |
+-------------------+----------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that uses linear continuous curvature interpolation and
weight function:
::

   EQ = momentum1 Q1 K Q1 1.0

-------------------------
**Technical Discussion**
-------------------------

Complete tutorial on the use of this equation exists. See GT-27.

----------
**Theory**
----------

The structural shell equation capability in Goma builds on the shell-element capability
built by Pat Notz and Ed Wilkes in FY03. Basically we are solving the following
equations for the shell tension and shell curvature (this card):

.. figure:: /figures/304_goma_physics.png
	:align: center
	:width: 90%


--------------
**References**
--------------

GT-27

..
	TODO - Line 66 contains a photo that needs to be written as an equation.
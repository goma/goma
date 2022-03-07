*******************
**shell_surf_curv**
*******************

::

	EQ =shell_surf_curv {Galerkin_wt} gamma2 {Interpol_fnc} <float1>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a definition equation for the total surface
curvature in a 2-dimensional bar element. Note that this equation is not yet available in
three dimensions and is in fact untested at this time. These building blocks are required
by the non-Newtonian surface rheology capability in Goma. Note that <floatlist>
contains one constant and it should always be set to one. The Galerkin weight and the
interpolation function must be the same for the code to work properly.

+-------------------+----------------------------------------------------------+
|**shell_surf_curv**|Name of the equation to be solved.                        |
+-------------------+----------------------------------------------------------+
|{Galerkin_wt}      |Two- or four-character value that defines the type of     |
|                   |weighting function for this equation, where:              |
|                   |                                                          |
|                   | * **Q1**-Linear                                          |
|                   | * **Q2**-Quadratic                                       |
+-------------------+----------------------------------------------------------+
|**gamma2**         |Name of the variable associated with the shell curvature  |
|                   |equation.                                                 |
+-------------------+----------------------------------------------------------+
|{Interpol_fnc}     |{Interpol_fnc} Two- or four-character value that defines  |
|                   |the interpolation function used to represent the variable |
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

   EQ = shell_surf_curv Q1 gamma2 Q1 1.0

-------------------------
**Technical Discussion**
-------------------------

See discussion for EQ = shell_surf_div_v



--------------
**References**
--------------

Edwards, D. A., Brenner, H., Wasan, D. T., 1991. Interfacial Transport Processes and
Rheology. Butterworth-Heinemann, Boston.
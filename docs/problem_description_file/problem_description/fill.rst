********
**fill**
********

::

	EQ = fill {Galerkin_wt} F {Interpol_fnc} <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for the fill equation.
Definitions of the input parameters are defined below. Note that <float1> through
<float3> define the constant multipliers for each term in the equation. The Galerkin
weight and the interpolation function must be the same for the code to work properly.

+----------------+--------------------------------------------------------------------+
|**fill**        |Name of the equation to be solved.                                  |
+----------------+--------------------------------------------------------------------+
|{Galerkin_wt}   |Two-character value that defines the type of weighting              |
|                |function for this equation, where:                                  |
|                |                                                                    |
|                | * **P0**-Constant Discontinuous                                    |
|                | * **P1**-Linear Discontinuous                                      |
|                | * **Q1**-Bilinear/Trilinear Continuous                             |
|                | * **Q2**-Biquadratic/Triquadratic Continuous                       |
|                | * **Q1_D**-Standard linear interpolation with special              |
|                |   allowance for discontinuous degrees of freedom at interfaces.    |
|                | * **Q2_D**-Standard quadratic interpolation with special           |
|                |   allowance for discontinuous degrees of freedom at interfaces.    |
|                | * **PQ1**-Q1 Discontinuous                                         |
|                | * **PQ2**-Q2 Discontinuous                                         |
+----------------+--------------------------------------------------------------------+
|**F**           |Name of the variable associated with this equation.                 |
+----------------+--------------------------------------------------------------------+
|{Interpol_fnc}  |Two-character value that defines the interpolation function         |
|                |used to represent the variable **F**, where:                        |
|                |                                                                    |
|                | * **P0**-Constant Discontinuous                                    |
|                | * **P1**-Linear Discontinuous                                      |
|                | * **Q1**-Bilinear/Trilinear Continuous                             |
|                | * **Q2**-Biquadratic/Triquadratic Continuous                       |
|                | * **Q1_D**-Standard linear interpolation with special              |
|                |   allowance for discontinuous degrees of freedom at interfaces.    |
|                | * **Q2_D**-Standard quadratic interpolation with special           |
|                |   allowance for discontinuous degrees of freedom                   |
|                |   at interfaces.                                                   |
|                | * **PQ1**-Q1 Discontinuous                                         |
|                | * **PQ2**-Q2 Discontinuous                                         |
+----------------+--------------------------------------------------------------------+
|<float1>        |Multiplier on mass matrix term ( d ⁄dt ).                           |
+----------------+--------------------------------------------------------------------+
|<float2>        |Multiplier on advective term.                                       |
+----------------+--------------------------------------------------------------------+
|<float3>        |Multiplier on source term.                                          |
+----------------+--------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that uses continuous linear interpolation for the fill
equation:
::

   EQ = fill Q1 F Q1 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

The fill equation is used in the calculation of volume of fluid interface tracking. It
solves an advection equation of a color function that takes on a different integer value
depending on which fluid phase you are in. For most applications this capability has
been superseded by the level set method of interface tracking.

The interpolation/weight functions that are discontinuous, e.g. have the prefix “P”
invoke the discontinuous Galerkin (DG) method for solving the fill equations where the
interpolation is discontinuous and flux continuity is maintained by evaluating surface
integrals. For details of the implementation of the DG method in *Goma* please see the
viscoelastic tutorial memo (Rao, 2000).



--------------
**References**
--------------

GT-014.1: Tutorial for Running Viscoelastic Flow Problems with GOMA, June 21,
2000, R. R. Rao
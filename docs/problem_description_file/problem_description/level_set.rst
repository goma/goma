*************
**level set**
*************

::

	EQ = level set {Galerkin_wt} F {Interpol_fnc} <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for the level set
equation. Definitions of the input parameters are defined below. Note that <float1>
through <float3> define the constant multipliers for each term in the equation. The
Galerkin weight and the interpolation function must be the same for the code to work
properly. If upwinding is desired, we can set this through a Petrov-Galerkin weight
function in the level set section of the input file (*Time Integration Specifications*).

+----------------+--------------------------------------------------------------------+
|**level set**   |Name of the equation to be solved.                                  |
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
|                |   allowance for discontinuous degrees of freedom at interfaces.    |
|                | * **PQ1**-Q1 Discontinuous                                         |
|                | * **PQ2**-Q2 Discontinuous                                         |
+----------------+--------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that uses continuous linear interpolation for the level set
equation and turns on all term multipliers:
::

   EQ = level_set Q1 F Q1 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

The interpolation/weight functions that are discontinuous, e.g. have the prefix “P”,
invoke the discontinuous Galerkin (DG) method for solving the level set equations
where the interpolation is discontinuous and flux continuity is maintained by
evaluating surface integrals. For details of the implementation of the DG method in
*Goma* please see the viscoelastic tutorial memo (Rao, 2000). Note that DG methods are
not necessarily recommended for the level set equation since it is inherently smooth.



--------------
**References**
--------------

GT-014.1: Tutorial for Running Viscoelastic Flow Problems with GOMA, June 21,
2000, R. R. Rao
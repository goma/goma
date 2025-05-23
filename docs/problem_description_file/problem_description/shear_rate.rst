**************
**shear_rate**
**************

::

	EQ = shear_rate {Galerkin_wt} SH {Interpol_fnc} <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for the scalar shear
rate invariant. Definitions of the input parameters are defined below. Note that <float1>
through <float3> define the constant multipliers in front of each type of term in the
equation. The Galerkin weight and the interpolation function must be the same for the
code to work properly.

+----------------+--------------------------------------------------------------------+
|**shear_rate**  |Name of the equation to be solved.                                  |
+----------------+--------------------------------------------------------------------+
|{Galerkin_wt}   |Two-character value that defines the type of weighting              |
|                |function for this equation, where:                                  |
|                |                                                                    |
|                | * **Q1**-Linear                                                    |
|                | * **Q2**-Quadratic                                                 |
+----------------+--------------------------------------------------------------------+
|**SH**          |Name of the variable associated with this equation.                 |
+----------------+--------------------------------------------------------------------+
|{Interpol_fnc}  |Two-character value that defines the interpolation function         |
|                |used to represent the variable **SH**, where:                       |
|                |                                                                    |
|                | * **Q1**-Linear                                                    |
|                | * **Q2**-Quadratic                                                 |
+----------------+--------------------------------------------------------------------+
|<float1>        |Multiplier on advective term.                                       |
+----------------+--------------------------------------------------------------------+
|<float2>        |Multiplier on diffusion term.                                       |
+----------------+--------------------------------------------------------------------+
|<float3>        |Multiplier on source term.                                          |
+----------------+--------------------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses quadratic continuous interpolation for the
species equation and turns on all the term multipliers:
::

   EQ = shear_rate Q2 SH Q2 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
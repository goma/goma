**********
**stress**
**********

::

	EQ = {eqname} {Galerkin_wt} {varname} {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation. Definitions of the
input parameters are defined below. Note that <floatlist> contains five constants for the
Stress equation defining the constant multipliers for each type of term in the equation.
The Galerkin weight and the interpolation function must be the same for the code to
work properly. If upwinding is desired for advection dominated problems, we can set
this through a Petrov-Galerkin weight function in the material file.

+---------------+--------------------------------------------------------------------+
|{eqname}       |The name of the component of the stress equation to be solved, one  |
|               |of the following: **stress11, stress12, stress13, stress22,         |
|               |stress23, stress33**.                                               |
+---------------+--------------------------------------------------------------------+
|{Galerkin_wt}  |Two-character or three-character value that defines the type        |
|               |of weighting function for this equation, where:                     |
|               |                                                                    |
|               | * **P0**-Constant Discontinuous                                    |
|               | * **P1**-Linear Discontinuous                                      |
|               | * **Q1**-Bilinear/Trilinear Continuous                             |
|               | * **Q2**-Biquadratic/Triquadratic Continuous                       |
|               | * **PQ**-Q1 Discontinuous                                          |
|               | * **PQ**-Q2 Discontinuous                                          |
+---------------+--------------------------------------------------------------------+
|{varname}      |The name of the variable associated with the respective             |
|               |components (11, 12, 13, 22, 23, and 33) of the symmetric            |
|               |Stress tensor, which are **S11, S12, S13, S22, S23, S33**.          |
+---------------+--------------------------------------------------------------------+
|{Interpol_fnc} |Two-character or three-character value that defines the             |
|               |interpolation function used to represent the variable               |
|               |**S11, S12, S13, S22, S23** or **S33**, where:                      |
|               |                                                                    |
|               | * **P0**-Constant Discontinuous                                    |
|               | * **P1**-Linear Discontinuous                                      |
|               | * **Q1**-Bilinear/Trilinear Continuous                             |
|               | * **Q2**-Biquadratic/Triquadratic Continuous                       |
|               | * **PQ**-Q1 Discontinuous                                          |
|               | * **PQ2**-Q2 Discontinuous                                         |
+---------------+--------------------------------------------------------------------+
|<float1>       |Multiplier on mass matrix term ( d ⁄dt ).                           |
+---------------+--------------------------------------------------------------------+
|<float2>       |Multiplier on advective term.                                       |
+---------------+--------------------------------------------------------------------+
|<float3>       |Multiplier on boundary term                                         |
|               |( :math:`\underline{n}` • flux  ).                                  |
+---------------+--------------------------------------------------------------------+
|<float4>       |Multiplier on diffusion term.                                       |
+---------------+--------------------------------------------------------------------+
|<float5>       |Multiplier on source term.                                          |
+---------------+--------------------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses a linear continuous interpolation for stress and
turns on all the term multipliers:
::

   EQ = stress11 Q1 S11 Q1 1. 1. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

The interpolation/weight functions that are discontinuous, e.g. have the prefix “P”,
invoke the discontinuous Galerkin method for solving the stress equations where the
interpolation is discontinuous and flux continuity is maintained by performing surface
integrals. For details of the implementation of the discontinuous Galerkin method in
*Goma* please see the viscoelastic tutorial memo (Rao, 2000).



--------------
**References**
--------------

GT-014.1: Tutorial for Running Viscoelastic Flow Problems with GOMA, June 21,
2000, R. R. Rao
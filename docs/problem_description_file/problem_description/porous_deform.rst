*****************
**porous_deform**
*****************

::

	EQ = porous_deform {Galerkin_wt} P_POR {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for porous solid
phase porosity. Definitions of the input parameters are defined below. Note that
<floatlist> has five parameters to define the constant multipliers in front of each type of
term in the equation. The Galerkin weight and the interpolation function must be the
same for the code to work properly. If upwinding is desired for advection dominated
problems, we can set this through a Petrov-Galerkin weight function in the material
file.

+-----------------+----------------------------------------------------------+
|**porous_deform**|Name of the equation to be solved.                        |
+-----------------+----------------------------------------------------------+
|{Galerkin_wt}    |Two-character value that defines the type of weighting    |
|                 |function for this equation, where:                        |
|                 |                                                          |
|                 | * **Q1**-Linear                                          |
|                 | * **Q2**-Quadratic                                       |
+-----------------+----------------------------------------------------------+
|**P_POR**        |Name of the variable associated with this equation.       |
+-----------------+----------------------------------------------------------+
|{Interpol_fnc}   |Two-character value that defines the interpolation        |
|                 |function used to represent the variable **P_POR**, where: |
|                 |                                                          |
|                 | * **Q1**-Linear Continuous                               |
|                 | * **Q2**-Quadratic Continuous                            |
+-----------------+----------------------------------------------------------+
|<float1>         |Multiplier on mass matrix term ( d ⁄dt ).                 |
+-----------------+----------------------------------------------------------+
|<float2>         |Multiplier on advective term.                             |
+-----------------+----------------------------------------------------------+
|<float3>         |Multiplier on boundary term                               |
|                 |( :math:`\underline{n}` • flux ).                         |
+-----------------+----------------------------------------------------------+
|<float4>         |Multiplier on diffusion term.                             |
+-----------------+----------------------------------------------------------+
|<float5>         |Multiplier on source term.                                |
+-----------------+----------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses a linear continuous interpolation and weight
function for the deforming porous porosity equation and has all the term multipliers on
except for the mass matrix for time derivatives:
::

   EQ = porous_deform Q1 P_POR Q1 0. 1. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk
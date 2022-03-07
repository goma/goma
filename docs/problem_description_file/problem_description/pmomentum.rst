*************
**pmomentum**
*************

::

	EQ = pmomentum{1|2|3} {Galerkin_wt} {PU1|PU2|PU3} {Interpol_fnc}<floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for one component of
a vector particle momentum equation. Definitions of the input parameters are defined
below. Note that <floatlist> contains six constants for the Pmomentum equation
defining the constant multipliers for each type of term in the equation. The Galerkin
weight and the interpolation function must be the same for the code to work properly.

+----------------------------------------+-------------------------------------------------+
|**pmomentum1 | pmomentum2 | pmomentum3**|Name of the equation to be solved, where the 1, 2|
|                                        |and 3 components correspond to one of the        |
|                                        |principal coordinate directions, e.g. X, Y and Z |
|                                        |for Cartesian geometry.                          |
+----------------------------------------+-------------------------------------------------+
|{Galerkin_wt}                           |Two- or four-character value that defines the    |
|                                        |type of weighting function for this equation,    |
|                                        |where:                                           |
|                                        |                                                 |
|                                        | * **Q1**-Linear                                 |
|                                        | * **Q2**-Quadratic                              |
|                                        | * **Q1_D**-Standard linear interpolation with   |
|                                        |   special allowance for discontinuous degrees of|
|                                        |   freedom at interfaces.                        |
|                                        | * **Q2_D**-Standard quadratic interpolation with|
|                                        |   special allowance for discontinuous degrees of|
|                                        |   freedom at interfaces.                        |
+----------------------------------------+-------------------------------------------------+
|**PU1 | PU2 | PU3**                     |Name of the variable associated with the 1, 2 or |
|                                        |3 principal coordinate direction for this        |
|                                        |component equation.                              |
+----------------------------------------+-------------------------------------------------+
|{Interpol_fnc}                          |Two- or four-character value that defines the    |
|                                        |interpolation function used to represent the     |
|                                        |variable **PU1**, **PU2** or **PU3** where:      |
|                                        |                                                 |
|                                        | * **Q1**-Linear Continuous                      |
|                                        | * **Q2**-Quadratic Continuous                   |
|                                        | * **Q1_D**-Standard linear interpolation with   |
|                                        |   special allowance for discontinuous degrees   |
|                                        |   of freedom at interfaces.                     |
|                                        | * **Q2_D**-Standard quadratic interpolation with|
|                                        |   special allowance for discontinuous degrees   |
|                                        |   of freedom at interfaces.                     |
+----------------------------------------+-------------------------------------------------+
|<float1>                                |Multiplier on mass matrix term ( d ⁄dt ).        |
+----------------------------------------+-------------------------------------------------+
|<float2>                                |Multiplier on advective term.                    |
+----------------------------------------+-------------------------------------------------+
|<float3>                                |Multiplier on boundary term                      |
|                                        |( :math:`\underline{n}` • flux  ).               |
+----------------------------------------+-------------------------------------------------+
|<float4>                                |Multiplier on diffusion term.                    |
+----------------------------------------+-------------------------------------------------+
|<float5>                                |Multiplier on source term.                       |
+----------------------------------------+-------------------------------------------------+
|<float6>                                |Multiplier on porous term (linear source).       |
+----------------------------------------+-------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses quadratic continuous velocity interpolation
and weight function and turns on all equation term multipliers except for the mass
matrix and the porous term:
::

   EQ = momentum1 Q2 PU1 Q2 0. 1. 1. 1. 1. 0.

-------------------------
**Technical Discussion**
-------------------------

The particle momentum equations have been added to *Goma* as part of a research
project and are not currently in use for production computing.



--------------
**References**
--------------

No References.
********
**mesh**
********

::

	EQ = mesh{1|2|3} {Galerkin_wt} {D1|D2|D3} {Interpol_fnc} <floatlist}

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for one component of
mesh motion. Definitions of the input parameters are defined below. Note that
<floatlist> contains five constants for the Mesh equation defining the constant
multipliers for each type of term in the equation. The Galerkin weight and the
interpolation function must be the same for the code to work properly.

+-------------------------+----------------------------------------------------------+
|**mesh1 | mesh2 | mesh3**|Name of the equation to be solved, where the 1, 2 and 3   |
|                         |components correspond to one of the principal coordinate  |
|                         |directions, e.g. X, Y and Z for Cartesian geometry.       |
+-------------------------+----------------------------------------------------------+
|{Galerkin_wt}            |Two-character value that defines the weighting function   |
|                         |type for this equation, where:                            |
|                         |                                                          |
|                         | * **Q1**-Linear                                          |
|                         | * **Q2**-Quadratic                                       |
+-------------------------+----------------------------------------------------------+
|**D1 | D2 | D3**         |Name of the variable associated with the 1, 2 or 3        |
|                         |principal coordinate direction for this component         |
|                         |equation.                                                 |
+-------------------------+----------------------------------------------------------+
|{Interpol_fnc}           |Two-character value that defines the interpolation        |
|                         |function used to represent the variable                   |
|                         |**D1, D2** or **D3** where:                               |
|                         |                                                          |
|                         | * **Q1**-Linear                                          |
|                         | * **Q2**-Quadratic                                       |
+-------------------------+----------------------------------------------------------+
|<float1>                 |Multiplier on mass matrix term ( d ⁄dt ).                 |
+-------------------------+----------------------------------------------------------+
|<float2>                 |Multiplier on advective term.                             |
+-------------------------+----------------------------------------------------------+
|<float3>                 |Multiplier on boundary term                               |
|                         |( :math:`\underline{n}` • flux  ).                        |
+-------------------------+----------------------------------------------------------+
|<float4>                 |Multiplier on diffusion term.                             |
+-------------------------+----------------------------------------------------------+
|<float5>                 |Multiplier on source term.                                |
+-------------------------+----------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card for the first mesh equation that uses linear continuous
interpolation and turns on all term multipliers except for the mass matrix:
::

   EQ = mesh1 Q1 D1 Q1 0. 1. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
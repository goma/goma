**************
**potential2**
**************

::

	EQ = potential2 {Galerkin_wt} PHI2 {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for the liquid-phase
electrolyte potential. This electrolyte-potential equation is solved together with the
solid-phase electrode-potential equation (see the potential1 card) for simulating
electrochemical processes (such as thermal batteries and proton-exchange-membrane
fuel cells) involving simultaneous charge transport in both the liquid-electrolyte and
solid-electrode phases (as in the porous anode and cathode). Definitions of the input
parameters are defined below. Note that <floatlist> has five parameters to define the
constant multipliers in front of each type of term in the equation. The Galerkin weight
and the interpolation function must be the same for the code to work properly.

+--------------------+----------------------------------------------------------+
|**potential2**      |Name of the equation to be solved.                        |
+--------------------+----------------------------------------------------------+
|{Galerkin_wt}       |Two-character value that defines the type of weighting    |
|                    |function for this equation, where:                        |
|                    |                                                          |
|                    | * **Q1**-Linear                                          |
|                    | * **Q2**-Quadratic                                       |
+--------------------+----------------------------------------------------------+
|**PHI2**            |Name of the variable associated with this equation.       |
+--------------------+----------------------------------------------------------+
|{Interpol_fnc}      |Two-character value that defines the interpolation        |
|                    |function used to represent the variable **PHI2**, where:  |
|                    |                                                          |
|                    | * **Q1**-Linear                                          |
|                    | * **Q2**-Quadratic                                       |
+--------------------+----------------------------------------------------------+
|<float1>            |Multiplier on mass matrix term ( d ⁄dt ).                 |
+--------------------+----------------------------------------------------------+
|<float2>            |Multiplier on advective term.                             |
+--------------------+----------------------------------------------------------+
|<float3>            |Multiplier on boundary term                               |
|                    |( :math:`\underline{n}` • flux ).                         |
+--------------------+----------------------------------------------------------+
|<float4>            |Multiplier on diffusion term.                             |
+--------------------+----------------------------------------------------------+
|<float5>            |Multiplier on source term.                                |
+--------------------+----------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses a quadratic interpolation and weight function
and turns off the mass (or transient) and advection terms but turns on the boundary,
diffusion, and source terms:
::

   EQ = potential2 Q2 PHI2 Q2 0. 0. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
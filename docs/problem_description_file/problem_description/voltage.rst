***********
**voltage**
***********

::

	EQ = voltage {Galerkin_wt} V {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for the voltage.
Definitions of the input parameters are defined below. Note that <floatlist> has five
parameters to define the constant multipliers in front of each type of term in the
equation. The Galerkin weight and the interpolation function must be the same for the
code to work properly.

+----------------+--------------------------------------------------------------------+
|**voltage**     |Name of the equation to be solved.                                  |
+----------------+--------------------------------------------------------------------+
|{Galerkin_wt}   |Two-character value that defines the type of weighting              |
|                |function for this equation, where:                                  |
|                |                                                                    |
|                | * **Q1**-Linear                                                    |
|                | * **Q2**-Quadratic                                                 |
+----------------+--------------------------------------------------------------------+
|**V**           |Name of the variable associated with this equation.                 |
+----------------+--------------------------------------------------------------------+
|{Interpol_fnc}  |Two-character value that defines the interpolation function         |
|                |used to represent the variable **V**, where:                        |
|                |                                                                    |
|                | * **Q1**-Linear                                                    |
|                | * **Q2**-Quadratic                                                 |
+----------------+--------------------------------------------------------------------+
|<float1>        |Multiplier on mass matrix term ( d ⁄dt ).                           |
+----------------+--------------------------------------------------------------------+
|<float2>        |Multiplier on advective term.                                       |
+----------------+--------------------------------------------------------------------+
|<float3>        |Multiplier on boundary term ( :math:`\underline{n}` • flux ).       |
+----------------+--------------------------------------------------------------------+
|<float4>        |Multiplier on diffusion term.                                       |
+----------------+--------------------------------------------------------------------+
|<float5>        |Multiplier on source term.                                          |
+----------------+--------------------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses a linear interpolation function for voltage:
::

   EQ = voltage Q1 V Q1 0. 1. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

The voltage equation has no mass term, viz. it is quasistatic. So it won’t matter
whether that multiplier is 1 or 0.



--------------
**References**
--------------

No References.
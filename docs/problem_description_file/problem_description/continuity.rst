**************
**continuity**
**************

::

	EQ = continuity {Galerkin_wt} P {Interpol_fnc} <float1> <float2>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation. Definitions of the
input parameters are defined below. Note that <float1> and <float2> define the
constant multipliers for each type of term in the Continuity equation. The Galerkin
weight and the interpolation function must be the same for the code to work properly.

+----------------+--------------------------------------------------------------------+
|**continuity**  |Name of the equation to be solved.                                  |
+----------------+--------------------------------------------------------------------+
|{Galerkin_wt}   |Two-character value that defines the type of weighting              |
|                |function for this equation, where:                                  |
|                |                                                                    |
|                | * **P0**-Constant Discontinuous                                    |
|                | * **P1**-Linear Discontinuous                                      |
|                | * **Q1**-Linear Continuous                                         |
|                | * **Q2**-Quadratic Continuous                                      |
|                | * **P0_XV**-Constant, discontinuous, enriched (level-set only)     |
|                | * **P1_XV**-Linear, discontinuous, enriched (level-set only)       |
|                | * **Q1_XV**-Linear, continuous, enriched (level-set only)          |
|                | * **Q2_XV**-Linear, continuous, enriched (level-set only)          |
+----------------+--------------------------------------------------------------------+
|**P**           |Name of the variable (pressure) associated with this equation.      |
+----------------+--------------------------------------------------------------------+
|{Interpol_fnc}  |Two-character value that defines the interpolation function         |
|                |used to represent the variable **P**, where:                        |
|                |                                                                    |
|                | * **P0**-Constant Discontinuous                                    |
|                | * **P1**-Linear Discontinuous                                      |
|                | * **Q1**-Linear Continuous                                         |
|                | * **Q2**-Quadratic Continuous                                      |
|                | * **P0_XV**-Constant, discontinuous, enriched (level-set only)     |
|                | * **P1_XV**-Linear, discontinuous, enriched (level-set only)       |
|                | * **Q1_XV**-Linear, continuous, enriched (level-set only)          |
|                | * **Q2_XV**-Linear, continuous, enriched (level-set only)          |
+----------------+--------------------------------------------------------------------+
|<float1>        |Multiplier on divergence term.                                      |
+----------------+--------------------------------------------------------------------+
|<float2>        |Multiplier on source term. This multiplier is equal to the          |
|                |initial volume fraction of solvents for Lagrangian mesh             |
|                |motion with swelling.                                               |
+----------------+--------------------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses a constant discontinuous pressure interpolation
and weight function and turns on both the divergence and source terms.
::

   EQ = continuity P0 P P0 1. 1.

-------------------------
**Technical Discussion**
-------------------------

Please see the EQ=energy equation card for a more detailed description of **P0_XV,
P1_XV, Q1_XV, Q2_XV** interpolations. These are MOST COMMONLY used for the
continuity equation for better accuracy of representing pressure across level-set
interfaces with surface tension.



--------------
**References**
--------------

No References.
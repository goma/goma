***************
**acous_pimag**
***************

::

	EQ = acous_pimag {Galerkin_wt} API {Interpol_fnc} <float list>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for the imaginary part
of the harmonic acoustic wave equation. Definitions of the input parameters are
defined below. Note that <float1> through <float5> define the constant multipliers in
front of each type of term in the equation. The Galerkin weight and the interpolation
function must be the same for the code to work properly.

+--------------------+----------------------------------------------------------+
|**acous_pimag**     |Name of the equation to be solved.                        |
+--------------------+----------------------------------------------------------+
|{Galerkin_wt}       |Two-character value that defines the type of weighting    |
|                    |function for this equation, where:                        |
|                    |                                                          |
|                    | * **Q1**-Linear                                          |
|                    | * **Q2**-Quadratic                                       |
+--------------------+----------------------------------------------------------+
|**API**             |Name of the variable associated with this equation.       |
+--------------------+----------------------------------------------------------+
|{Interpol_fnc}      |Two-character value that defines the interpolation        |
|                    |function used to represent the variable **SH**, where:    |
|                    |                                                          |
|                    | * **Q1**-Linear                                          |
|                    | * **Q2**-Quadratic                                       |
+--------------------+----------------------------------------------------------+
|<float1>            |currently not used.                                       |
+--------------------+----------------------------------------------------------+
|<float2>            |Multiplier on acoustic absorption term.                   |
+--------------------+----------------------------------------------------------+
|<float3>            |Multiplier on boundary terms.                             |
+--------------------+----------------------------------------------------------+
|<float4>            |Multiplier on Laplacian term.                             |
+--------------------+----------------------------------------------------------+
|<float5>            |Multiplier on pressure term.                              |
+--------------------+----------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses quadratic continuous interpolation for the
species equation and turns on all the term multipliers:
::

   EQ = acous_pimag Q2 API Q2 0. 1. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

Harmonic form of the wave equation with absorption (attenuation) included. P is the
amplitude of the acoustic pressure (complex), k is the wavenumber, :math:`\alpha` is the absorption
coefficient, and :math:`\omega` is the frequency (rad/sec).

.. figure:: /figures/308_goma_physics.png
	:align: center
	:width: 90%




..
	TODO - Line 69 contains a photo that needs to be written as an equation.
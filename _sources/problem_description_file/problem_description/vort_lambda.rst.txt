***************
**vort_lambda**
***************

::

	EQ = vort_lambda {Galerkin_wt} VLAMBDA {Interpol_fnc}

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for the vorticity
direction. Definitions of the input parameters are defined below; there are no <float>
input parameters for this equation. The Galerkin weight and the interpolation function
must be the same for the code to work properly.

+---------------+----------------------------------------------------------+
|**vort_lamda** |Name of the equation to be solved.                        |
+---------------+----------------------------------------------------------+
|{Galerkin_wt}  |Two-character value that defines the type of weighting    |
|               |function for this equation, where:                        |
|               |                                                          |
|               | * **Q1**-Linear                                          |
|               | * **Q2**-Quadratic                                       |
+---------------+----------------------------------------------------------+
|**VLAMBDA**    |Name of the variable associated with this equation.       |
+---------------+----------------------------------------------------------+
|{Interpol_fnc} |Two-character value that defines the interpolation        |
|               |function used to represent the variable **VLAMBDA**,      |
|               |where:                                                    |
|               |                                                          |
|               | * **Q1**-Linear Continuous                               |
|               | * **Q2**-Quadratic Continuous                            |
+---------------+----------------------------------------------------------+

------------
**Examples**
------------

No Examples.

-------------------------
**Technical Discussion**
-------------------------

This equation type is used for a research capability involving the flows of suspensions
in curvilinear coordinates and is not currently being used for production computations.




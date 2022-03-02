***************
**shell_filmh**
***************

::

	EQ = shell_filmh {Galerkin_wt} SHELL_FILMH {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving the film lubrication equation for free surface
flow. Definitions of the input parameters are defined below. The Galerkin weight and
the interpolation function must be set the same for the code to work properly.

Counterparts to this equation for lubrication flow of capillary films (film-equations) are
lup_p equation.

+--------------------+----------------------------------------------------------+
|**shell_filmh**     |Name of equation to be solved.                            |
+--------------------+----------------------------------------------------------+
|{Galerkin_wt}       |Two-or four-character value that defines the type of      |
|                    |weighting function for this equation, where:              |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|**SHELL_FILMH**     |Name of the variable associated with this equation.       |
+--------------------+----------------------------------------------------------+
|{Interpol_fnc}      |Two-or four-character value that defines the              |
|                    |interpolation function for the variable **SHELL_FILMH**,  |
|                    |where:                                                    |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|<float1>            |Multiplier for the mass matrix term. It is not activated. |
+--------------------+----------------------------------------------------------+
|<float2>            |Multiplier for the advection term. It is not activated.   |
+--------------------+----------------------------------------------------------+
|<float3>            |Multiplier for the boundary term. It is not activated.    |
+--------------------+----------------------------------------------------------+
|<float4>            |Multiplier for the diffusion term.                        |
+--------------------+----------------------------------------------------------+
|<float5>            |Multiplier for the source term.                           |
+--------------------+----------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:
::

   EQ = shell_filmh Q1 SHELL_FILMH Q1 0. 0. 0. 1. 1.

This applies the film flow equation with all terms activated.

-------------------------
**Technical Discussion**
-------------------------

The equation solved is as follows:

.. figure:: /figures/314_goma_physics.png
	:align: center
	:width: 90%

* The diffusion multiplier activates the capillary pressure term.

* The source activates the first term.

* This equation does not fit the general prototype of conservation equation where the
  diffusion and source terms really apply. In all cases, both diffusion and source
  terms need to be activated.

* This equation has to be used with the equation describing SHELL_FILMP.




..
	TODO - Line 66 contains a photo that needs to be written as an equation.
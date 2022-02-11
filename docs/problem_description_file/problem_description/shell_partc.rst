***************
**shell_partc**
***************

::

	EQ = shell_partc {Galerkin_wt} SHELL_PARTC {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving the z-averaged concentration of particles
inside film flow. Definitions of the input parameters are defined below. The Galerkin
weight and the interpolation function must be set the same for the code to work
properly.

+--------------------+----------------------------------------------------------+
|**shell_partc**     |Name of equation to be solved.                            |
+--------------------+----------------------------------------------------------+
|{Galerkin_wt}       |Two-or four-character value that defines the type of      |
|                    |weighting function for this equation, where:              |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|**SHELL_PARTC**     |Name of the variable associated with this equation.       |
+--------------------+----------------------------------------------------------+
|{Interpol_fnc}      |Two-or four-character value that defines the              |
|                    |interpolation function for the variable                   |
|                    |**SHELL_PARTC**, where:                                   |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|<float1>            |Multiplier for the mass matrix term.                      |
+--------------------+----------------------------------------------------------+
|<float2>            |Multiplier for the advection term.                        |
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

   EQ = shell_partc Q1 SHELL_PARTC Q1 1. 0. 0. 1. 0.

This applies the film flow equation with all terms activated.

-------------------------
**Technical Discussion**
-------------------------

The equation solved is as follows:

.. figure:: /figures/315_goma_physics.png
	:align: center
	:width: 90%

* The mass matrix multiplier activates the time-derivative term.

* The advection multiplier activates the second term, where the flow rate is dotted
  onto the gradient of particles concentration and it represents advection of particles
  due to the liquid film flow.

* The diffusion multiplier activates the terms inside the divergence operator and
  represents the Fickian diffusion of particles.

* The source activates the last term, rate of evaporation of liquid that contributes to
  the increase of the particles conentration.

* This equation has to be used with the film profile equation describing
  SHELL_FILMP and SHELL_FILMH.




..
	TODO - Line 64 contains a photo that needs to be written as an equation.
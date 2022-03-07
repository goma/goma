********************
**shell_sat_open_2**
********************

::

	EQ = shell_sat_open_2 {Galerkin_wt} SH_P_OPEN_2 {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides the capability to solve a second porous shell equation for open
(interconnected) structured pores. The use of this equation requires that the shell
material share the same nodes but be a distinct material from that which shell_sat_open
resides. Please see the associated tutorials. The Galerkin weight and the interpolation
function must be set the same for the code to work properly. The counterpart to this
equation is shell_sat_closed, which solves for non-interconnected pores.

+--------------------+----------------------------------------------------------+
|**shell_sat_open_2**|Name of equation to be solved.                            |
+--------------------+----------------------------------------------------------+
|{Galerkin_wt}       |Two-or four-character value that defines the type of      |
|                    |weighting function for this equation, where:              |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|**SH_SAT_OPEN_2**   |Name of the variable associated with this equation.       |
+--------------------+----------------------------------------------------------+
|{Interpol_fnc}      |Two-or four-character value that defines the              |
|                    |interpolation function for the variable **SH_P_OPEN_2**,  |
|                    |where:                                                    |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|<float1>            |Multiplier for the mass matrix term.                      |
+--------------------+----------------------------------------------------------+
|<float2>            |Multiplier for the source term.                           |
+--------------------+----------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:
::

   EQ = shell_sat_open_2 Q1 SH_P_OPEN_2 Q1   1.0 1.0

This applies the equation with all terms activated.

-------------------------
**Technical Discussion**
-------------------------

The equation solved is as follows:

.. figure:: /figures/318_goma_physics.png
	:align: center
	:width: 90%

* The mass matrix multiplier activates the time-derivative term.

* The source matrix multiplier activates the remaining term.

* This equation is required to couple with LUBP to solve for the lubrication forces.

* Currently, this equation assumes that the porous shell is located in the +z direction
  of the lubrication shell, and the coupling is set up to draw liquid from the
  lubrication layer by adding a sink term into the lubrication equations.

* NOT FULLY IMPLEMENTED.




..
	TODO - Line 60 contains a photo that needs to be written as an equation.
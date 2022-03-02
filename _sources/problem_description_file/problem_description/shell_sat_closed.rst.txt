********************
**shell_sat_closed**
********************

::

	EQ = shell_sat_closed {Galerkin_wt} SH_SAT_CLOSED {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides the capability to solve the porous shell equations for closed (noninterconnected) structured pores. The Galerkin weight and the interpolation function
must be set the same for the code to work properly. The counterpart to this equation is
porous_sat_open, which solves for interconnected pores.

+--------------------+----------------------------------------------------------+
|**shell_sat_closed**|Name of equation to be solved.                            |
+--------------------+----------------------------------------------------------+
|{Galerkin_wt}       |Two-or four-character value that defines the type of      |
|                    |weighting function for this equation, where:              |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|**SH_SAT_CLOSED**   |Name of the variable associated with this equation.       |
+--------------------+----------------------------------------------------------+
|{Interpol_fnc}      |Two-or four-character value that defines the              |
|                    |interpolation function for the variable                   |
|                    |**SH_SAT_CLOSED**, where:                                 |
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

   EQ = shell_sat_closed   Q1 SH_SAT_CLOSED Q1   1.0 1.0

This applies the equation with all terms activated.

-------------------------
**Technical Discussion**
-------------------------

The equation solved is as follows:

.. figure:: /figures/316_goma_physics.png
	:align: center
	:width: 90%

* The mass matrix multiplier activates the time-derivative term.

* The source matrix multiplier activates the remaining term.

* This equation is required to couple with LUBP to solve for the lubrication forces.

* Currently, this equation assumes that the porous shell is located in the +z direction
  of the lubrication shell, and the coupling is set up to draw liquid from the
  lubrication layer by adding a sink term into the lubrication equations.

Beyond the standard porous media material cards for continuum element regions, one
needs in the thin-shell material inputs in the following section:

::

   Porous Shell Closed Porosity = CONSTANT 0.1

::

   Porous Shell Height = CONSTANT 1.0

::

   Porous Shell Radius = CONSTANT 0.01

::

   Porous Shell Atmospheric Pressure = CONSTANT 1.e6

::

   Porous Shell Reference Pressure = CONSTANT 0.

::

   Porous Shell Cross Permeability = CONSTANT 0.2

::

   Porous Shell Initial Pore Pressure = CONSTANT 0.

Please read the associated material property cards sections for details.




..
	TODO - Line 57 contains a photo that needs to be written as an equation.
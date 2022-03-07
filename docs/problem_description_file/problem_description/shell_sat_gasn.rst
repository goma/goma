******************
**shell_sat_gasn**
******************

::

	EQ = shell_sat_gasn {Galerkin_wt} SH_SAT_GASN {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides the capability to solve the porous shell equation for the inventory of
trapped gas in a closed pore shell simulation, viz. the EQ=shell_sat_closed. The
equation tracks the inventory of trapped gas and accounts for the compression (ideal
gas law) and dissolution into the invading liquid. Two terms are required in this
equation:

+--------------------+----------------------------------------------------------+
|**shell_sat_gasn**  |Name of equation to be solved.                            |
+--------------------+----------------------------------------------------------+
|{Galerkin_wt}       |Two-or four-character value that defines the type of      |
|                    |weighting function for this equation, where:              |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|**SH_SAT_GASN**     |Name of the variable associated with this equation.       |
+--------------------+----------------------------------------------------------+
|{Interpol_fnc}      |Two-or four-character value that defines the              |
|                    |interpolation function for the variable                   |
|                    |**SH_SAT_GASN**, where:                                   |
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

   EQ = shell_sat_gasn Q1 SH_SAT_GASN Q1   1.0 1.0

This applies the equation with all terms activated.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

S. A. Roberts and P. R. Schunk 2012. In preparation.
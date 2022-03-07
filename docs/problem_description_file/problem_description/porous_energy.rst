*****************
**porous_energy**
*****************

::

	EQ = porous_energy {Galerkin_wt} P_TEMP {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a conservation of energy differential
equation for porous media, deploying a multiphase formulation.. Definitions of the
input parameters are defined below. Note that <floatlist> contains six constants for the
porous energy equation defining the constant multipliers for each term in the equation.
The Galerkin weight and the interpolation function must be the same for the code to
work properly.

+-----------------+----------------------------------------------------------+
|**porous_energy**|Name of the equation to be solved.                        |
+-----------------+----------------------------------------------------------+
|{Galerkin_wt}    |Two- or four-character value that defines the type of     |
|                 |weighting function for this equation, where:              |
|                 |                                                          |
|                 | * **Q1**-Linear                                          |
|                 | * **Q2**-Quadratic                                       |
+-----------------+----------------------------------------------------------+
|**T**            |Name of the variable associated with this equation.       |
+-----------------+----------------------------------------------------------+
|{Interpol_fnc}   |Two- or four-character value that defines the             |
|                 |interpolation function used to represent the variable     |
|                 |**T**, where:                                             |
|                 |                                                          |
|                 | * **Q1**-Linear Continuous                               |
|                 | * **Q2**-Quadratic Continuous                            |
+-----------------+----------------------------------------------------------+
|<float1>         |Multiplier on mass matrix term ( d ⁄dt ).                 |
+-----------------+----------------------------------------------------------+
|<float2>         |Multiplier on advective term.                             |
+-----------------+----------------------------------------------------------+
|<float3>         |Multiplier on boundary term                               |
|                 |( :math:`\underline{n}` • flux ).                         |
+-----------------+----------------------------------------------------------+
|<float4>         |Multiplier on diffusion term.                             |
+-----------------+----------------------------------------------------------+
|<float5>         |Multiplier on source term.                                |
+-----------------+----------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses a linear continuous interpolation and weight
function and has all the term multipliers on:
::

   porous_energy Q1 P_TEMP Q1     1.   1.   1.   1.   1.

-------------------------
**Technical Discussion**
-------------------------

Usage of this equation is discussed extensively in GT-009.3 Output variables in the
ExodusII database are POR_TEMP



--------------
**References**
--------------

GT-009
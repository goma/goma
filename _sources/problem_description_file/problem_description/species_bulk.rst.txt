****************
**species_bulk**
****************

::

	EQ = species_bulk {Galerkin_wt} Y {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation. Definitions of the
input parameters are defined below. Note that <floatlist> contains five parameters to
define the constant multipliers in front of each type of term in the equation. The
Galerkin weight and the interpolation function must be the same for the code to work
properly. If upwinding is desired for advection dominated problems, we can set this
through a Petrov-Galerkin weight function in the material file.

+----------------+--------------------------------------------------------------------+
|**species_bulk**|Name of the equation to be solved. This equation type               |
|                |should only be listed once regardless of the number of              |
|                |species (the *Number of bulk species* card specifies the            |
|                |number of species_bulk equations to be solved).                     |
|                |Differences in diffusion coefficients between species should        |
|                |be accounted for in the materials properties section of *Goma*.     |
+----------------+--------------------------------------------------------------------+
|{Galerkin_wt}   |Two- to four-character value that defines the type of               |
|                |weighting function for this equation, where:                        |
|                |                                                                    |
|                | * **P0**-Constant Discontinuous                                    |
|                | * **P1**-Linear Discontinuous                                      |
|                | * **Q1**-Bilinear/Trilinear Continuous                             |
|                | * **Q2**-Biquadratic/Triquadratic Continuous                       |
|                | * **Q1_D**-Standard linear interpolation with special              |
|                |   allowance for discontinuous degrees of freedom at interfaces.    |
|                | * **Q2_D**-Standard quadratic interpolation with special           |
|                |   allowance for discontinuous degrees of freedom at interfaces.    |
|                | * **PQ1**-Q1 Discontinuous                                         |
|                | * **PQ2**-Q2 Discontinuous                                         |
|                | * **Q1_XV, Q1_GN, Q1_GP**-Linear interpolation with enrichment in  |
|                |   elements of material interfaces. This enrichment function        |
|                |   allows discontinuity in value and gradient along                 |
|                |   interface but maintains continuity at element edges/faces.       |
|                | * **Q2_XV, Q2_GN, Q1_GP**-Quadratic interpolation with enrichment  |
|                |   in elements of material interfaces. This enrichment              |
|                |   function allows discontinuity in value and gradient              |
|                |   along interface but maintains continuity at element edges/faces. |
+----------------+--------------------------------------------------------------------+
|**Y**           |Name of the variable associated with this equation.                 |
+----------------+--------------------------------------------------------------------+
|{Interpol_fnc}  |Two- to four-character value that defines the interpolation         |
|                |function used to represent the variable **Y**, where:               |
|                |                                                                    |
|                | * **P0**-Constant Discontinuous                                    |
|                | * **P1**-Linear Discontinuous                                      |
|                | * **Q1**-Bilinear/Trilinear Continuous                             |
|                | * **Q2**-Biquadratic/Triquadratic Continuous                       |
|                | * **Q1_D**-Standard linear interpolation with special              |
|                |   allowance for discontinuous degrees of freedom at interfaces.    |
|                | * **Q2_D**-Standard quadratic interpolation with special           |
|                |   allowance for discontinuous degrees of freedom at interfaces.    |
|                | * **PQ1**-Q1 Discontinuous                                         |
|                | * **PQ2**-Q2 Discontinuous                                         |
|                | * **Q1_XV, Q1_GN, Q1_GP**-Linear interpolation with enrichment in  |
|                |   elements of material interfaces. This enrichment                 |
|                |   function allows discontinuity in value and                       |
|                |   gradient along interface but maintains continuity                |
|                |   at element edges/faces. *See energy equation for                 |
|                |   more discussion*.                                                |
|                | * **Q2_XV, Q2_GN, Q1_GP**-Quadratic interpolation with enrichment  |
|                |   in elements of material interfaces. This enrichment              |
|                |   function allows discontinuity in value and                       |
|                |   gradient along interface but maintains continuity                |
|                |   at element edges/faces. *See energy equation for                 |
|                |   more discussion.*                                                |
+----------------+--------------------------------------------------------------------+
|<float1>        |Multiplier on mass matrix term ( d ⁄dt ).                           |
+----------------+--------------------------------------------------------------------+
|<float2>        |Multiplier on advective term.                                       |
+----------------+--------------------------------------------------------------------+
|<float3>        |Multiplier on boundary term                                         |
|                |( :math:`\underline{n}` • flux  ).                                  |
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

The following is a sample card that uses quadratic continuous interpolation for the
species equation and turns on all the term multipliers:
::

   EQ = species_bulk Q2 Y Q2 1. 1. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

The interpolation/weight functions that are discontinuous, e.g. have the prefix “P”,
invoke the discontinuous Galerkin (DG) method for solving the species equations
where the interpolation is discontinuous and flux continuity is maintained by
performing surface integrals. For details of the implementation of the DG method in
*Goma* please see the viscoelastic tutorial memo. Note, the DG implementation for the
species equation is only for advection dominated problems; DG methods have not yet
been completely developed for diffusion operators.

Also, please see EQ=energy input for more detailed description of the Q1_GN,
Q2_GN, Q1_GP, Q2_GP, Q1_XV and Q2_XV enriched basis functions.



--------------
**References**
--------------

GT-014.1: Tutorial for Running Viscoelastic Flow Problems with GOMA, June 21,
2000, R. R. Rao
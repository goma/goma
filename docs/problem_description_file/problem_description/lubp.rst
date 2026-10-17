********
**lubp**
********

::

	EQ = lubp {Galerkin_wt} LUBP {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving the Reynold’s lubrication equation for
confined flow. Definitions of the input parameters are defined below. The Galerkin
weight and the interpolation function must be set the same for the code to work
properly. Counterparts to this equation for lubrication flow of capillary films (filmequations) are shell_filmp and shell_filmh equations.

+--------------------+----------------------------------------------------------+
|**lubp**            |Name of equation to be solved.                            |
+--------------------+----------------------------------------------------------+
|{Galerkin_wt}       |Two-or four-character value that defines the type of      |
|                    |weighting function for this equation, where:              |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|**LUBP**            |Name of the variable associated with this equation.       |
+--------------------+----------------------------------------------------------+
|{Interpol_fnc}      |Two-or four-character value that defines the              |
|                    |interpolation function for the variable **LUBP**, where:  |
|                    |                                                          |
|                    | * **Q1**–Linear                                          |
|                    | * **Q2**–Quadratic (not recommended at this time)        |
+--------------------+----------------------------------------------------------+
|<float1>            |Multiplier for the mass matrix term (not used yet as of   |
|                    |3/4/2010).                                                |
+--------------------+----------------------------------------------------------+
|<float2>            |Multiplier for the diffusion term.                        |
+--------------------+----------------------------------------------------------+
|<float3>            |Multiplier for the source term.                           |
+--------------------+----------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:
::

   EQ = lubp Q1 LUBP Q1 1. 1. 1.

This applies the confined flow lubrication equation with all terms activated.

-------------------------
**Technical Discussion**
-------------------------

The equation solved is as follows:

.. math::

   \frac{\partial(\rho h)}{\partial t}+\nabla_{II}\cdot\left(\frac{\rho h}{2}(U_A+U_B)-\frac{\rho h^3}{12\mu}\nabla_{II}p-\frac{\rho h^3}{12\mu}\sigma\kappa\delta(\phi)\boldsymbol{n}-\rho(\phi)\boldsymbol{g}\right)-(j_A+j_B)=0

* The first term multiplier, activating the mass (time-derivative) term is not currently
  activated as the gap-height is user-prescribed.

* The second term multiplier affects the third and fourth terms (grad_p and surface
  tension terms).

The third term multiplier activates the Couette flow terms.

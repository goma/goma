****************
**shell_normal**
****************

::

	EQ =shell_normal{1|2} {Galerkin_wt} {SH_N1|SH_N2} {Interpol_fnc} <float1>

-----------------------
**Description / Usage**
-----------------------

This card specifies a vector of shell normal vector component unknowns in a 2-
dimensional bar element, intended for use with shell diffusive flux problems. Note that
this equation is not yet available in three dimensions and is in fact untested at this time.
Note that <floatlist> contains one constant and it should always be set to one. The
Galerkin weight and the interpolation function must be the same for the code to work
properly.

+---------------------------------+--------------------------------------------+
|**shell_normal1 | shell_normal2**|Name of the equation to be solved.          |
+---------------------------------+--------------------------------------------+
|{Galerkin_wt}                    |Two- or four-character value that defines   |
|                                 |the type of weighting function for this     |
|                                 |equation, where:                            |
|                                 |                                            |
|                                 | * **Q1**-Linear                            |
|                                 | * **Q2**-Quadratic                         |
+---------------------------------+--------------------------------------------+
|**SH_N1 | SH_N2**                |Name of the variable associated with the    |
|                                 |shell curvature equation.                   |
+---------------------------------+--------------------------------------------+
|{Interpol_fnc}                   |Two- or four-character value that defines   |
|                                 |the interpolation function used to represent|
|                                 |the variable **SH_N1(2)** where:            |
|                                 |                                            |
|                                 | * **Q1**-Linear Continuous                 |
|                                 | * **Q2**-Quadratic Continuous              |
+---------------------------------+--------------------------------------------+
|<float1>                         |Multiplier on diffusion term (i.e. the whole|
|                                 |equation). Set to 1.0.                      |
+---------------------------------+--------------------------------------------+

------------
**Examples**
------------

The following is a pair of sample cards that use linear continuous normal interpolation
and weight function:
::

   EQ = shell_normal1 Q1 SH_N1 Q1 1.0

::

   EQ = shell_normal2 Q1 SH_N2 Q1 1.0

Note that since this equation applies only to 2D problem domains at this time, two
cards are needed as shown above (one for each component).

-------------------------
**Technical Discussion**
-------------------------

This equation merely sets the components of the shell normal vector equal to those in
fv->snormal, which are calculated rigorously in surface_determinant_and_normal().
Consideration is being given to replacing these with a single unknown for shell normal
angle, which contains the same information in a single scalar unknown.



--------------
**References**
--------------

No References.
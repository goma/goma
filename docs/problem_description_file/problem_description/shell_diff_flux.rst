*******************
**shell_diff_flux**
*******************

::

	EQ =shell_diff_flux {Galerkin_wt} SH_J {Interpol_fnc} <float1>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a conservation equation for the total surface
diffusive flux in a 2-dimensional bar (or shell) element. Note that this equation is not
yet available in three dimensions and is in fact untested at this time. The card entries
are as follows:

+-------------------+----------------------------------------------------------+
|**shell_diff_flux**|Name of the equation to be solved.                        |
+-------------------+----------------------------------------------------------+
|{Galerkin_wt}      |Two- or four-character value that defines the type of     |
|                   |weighting function for this equation, where:              |
|                   |                                                          |
|                   | * **Q1**-Linear                                          |
|                   | * **Q2**-Quadratic                                       |
+-------------------+----------------------------------------------------------+
|**SH_J**           |Name of the variable associated with the shell curvature  |
|                   |equation.                                                 |
+-------------------+----------------------------------------------------------+
|{Interpol_fnc}     |Two- or four-character value that defines the             |
|                   |interpolation function used to represent the variable     |
|                   |**SH_J** where:                                           |
|                   |                                                          |
|                   | * **Q1**-Linear Continuous                               |
|                   | * **Q2**-Quadratic Continuous                            |
+-------------------+----------------------------------------------------------+
|<float1>           |Multiplier for diffusion terms (in this case, the whole   |
|                   |equation).                                                |
+-------------------+----------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that uses bilinear shell diffusive flux interpolation and
weight function:
::

   EQ = shell_diff_flux Q2 SH_J Q2 1.0

-------------------------
**Technical Discussion**
-------------------------

This is only a preliminary implementation of a shell quantity conservation equation. It
is not currently operational. When it is fully implemented, the number of required
equation term multiplier entries will be adjusted acordingly.



--------------
**References**
--------------

No References.
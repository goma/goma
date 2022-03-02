*******************
**shell_diff_curv**
*******************

::

	EQ =shell_diff_curv {Galerkin_wt} SH_KD {Interpol_fnc} <float1>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a definition equation for the total surface
curvature in a 2-dimensional bar element, intended for use with shell diffusive flux
problems. Note that this equation is not yet available in three dimensions and is in fact
untested at this time. Note that <floatlist> contains one constant and it should always be
set to one. The Galerkin weight and the interpolation function must be the same for the
code to work properly.

+-------------------+----------------------------------------------------------+
|**shell_diff_curv**|Name of the equation to be solved.                        |
+-------------------+----------------------------------------------------------+
|{Galerkin_wt}      |Two- or four-character value that defines the type of     |
|                   |weighting function for this equation, where:              |
|                   |                                                          |
|                   | * **Q1**-Linear                                          |
|                   | * **Q2**-Quadratic                                       |
+-------------------+----------------------------------------------------------+
|**SH_KD**          |Name of the variable associated with the shell curvature  |
|                   |equation.                                                 |
+-------------------+----------------------------------------------------------+
|{Interpol_fnc}     |Two- or four-character value that defines the             |
|                   |interpolation function used to represent the variable     |
|                   |**SH_KD** where:                                          |
|                   |                                                          |
|                   | * **Q1**-Linear Continuous                               |
|                   | * **Q2**-Quadratic Continuous                            |
+-------------------+----------------------------------------------------------+
|<float1>           |Multiplier on diffusion term (i.e. the whole equation).   |
|                   |Set to 1.0.                                               |
+-------------------+----------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that uses linear continuous curvature interpolation and
weight function:
::

   EQ = shell_diff_curv Q1 SH_KD Q1 1.0

-------------------------
**Technical Discussion**
-------------------------

The equation solved is the surface curvature definition :math:`\kappa` = :math:`\Delta_s` :math:`\underline{n}` = (I – :math:`\underline{n}` :math:`\underline{n}`) • :math:`\Delta` n. See
discussion for EQ = shell_surf_div_v.



--------------
**References**
--------------

Edwards, D. A., Brenner, H., Wasan, D. T., 1991. Interfacial Transport Processes and
Rheology. Butterworth-Heinemann, Boston.
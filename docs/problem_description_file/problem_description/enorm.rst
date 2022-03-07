*********
**enorm**
*********

::

	EQ = enorm {Galerkin_wt} ENORM {Interpol_fnc} <float1> <float2>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a “dependency” equation for the norm of the
electric field. Definitions of the input parameters are defined below. Note that <float1>
and <float2> define the constant multipliers in front of each type of term in the
equation. The Galerkin weight and the interpolation function must be the same for the
code to work properly.

+----------------+--------------------------------------------------------------------+
|**enorm**       |Name of the equation to be solved.                                  |
+----------------+--------------------------------------------------------------------+
|{Galerkin_wt}   |Two-character value that defines the type of weighting              |
|                |function for this equation, where:                                  |
|                |                                                                    |
|                | * **P0**-Piecewise constant                                        |
|                | * **P1**-Piecewise linear                                          |
|                | * **Q1**-Linear                                                    |
|                | * **Q2**-Quadratic                                                 |
+----------------+--------------------------------------------------------------------+
|**ENORM**       |Name of the variable associated with this equation.                 |
+----------------+--------------------------------------------------------------------+
|{Interpol_fnc}  |Two-character value that defines the interpolation function         |
|                |used to represent the variable **ENORM**, where:                    |
|                |                                                                    |
|                | * **P0**-Piecewise constant                                        |
|                | * **P1**-Piecewise linear                                          |
|                | * **Q1**-Linear                                                    |
|                | * **Q2**-Quadratic                                                 |
+----------------+--------------------------------------------------------------------+
|<float1>        |Multiplier on advection term.                                       |
+----------------+--------------------------------------------------------------------+
|<float2>        |Multiplier on source term.                                          |
+----------------+--------------------------------------------------------------------+

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped. See below for
important information regarding this.

------------
**Examples**
------------

The following is a sample card that uses quadratic continuous interpolation for the
enorm equation and turns on all the term multipliers (the usual usage):
::

   EQ = enorm Q2 ENORM Q2 1.0 1.0

-------------------------
**Technical Discussion**
-------------------------

This equation allows the user to use the variable ENORM, the norm of the electric
field, which is equal to :math:`\mid` :math:`\underline{E}` :math:`\mid`, or :math:`\mid` :math:`\underline \Delta` V :math:`\mid`, with V being the voltage potential. As such, the
VOLTAGE equation must be present. We refer to this as a “dependent” equation or
“auxiliary” equation because, although it’s value can technically be derived from the V
variable directly, we would lose derivative information by doing so. This equation is
introduced solely so one can access higher derivatives of V than its interpolation would
normally allow. For example, V if were interpolated with a linear basis, then 
:math:`\underline \Delta` V would
have a constant interpolant. If we wanted access to :math:`\underline \Delta` 
(:math:`\underline \Delta` V) , it would be zero! (In reality,
we would use bilinear or trilinear basis functions, so this isn’t precisely true but it
expresses the essential problem.) By introducing this primitive variable, we can
retrieve useful values for :math:`\underline \Delta`_{enorm}`.

The two term multipliers refer to the multiple on the assembled value of enorm (stored
in the “advection” term--it has nothing to do with advection), and the multiple on the
assembled value derived from the voltage equation (stored in the “source” term--again
the name of the term is somewhat artificial).



--------------
**References**
--------------

No References.
**************************
Polymer Stress Formulation
**************************

::

   Polymer Constitutive Equation = {model_name}

-----------------------
**Description / Usage**
-----------------------

This card specifies which formulation of the polymer constitutive equation should be
used. Valid options are

+-----------------+------------------------------------------------------------------------------------------------------------+
|**EVSS_G**       |Uses the classic elastic-viscous stress splitting of Rajagopalan (1990) where the stress is the elastic     |
|                 |stress only without a Newtonian component. This option is the default if this *Polymer Stress Formulation*  |
|                 |card is not supplied. This formulation is almost never used.                                                |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**EVSS_F**       |Uses the EVSS formulation of Guenette and Fortin (1995) that solves the standard stress equation with the   |
|                 |addition of a new term to the momentum equation. This formulation is used most often.                       |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**EVSS_L**       |Uses a research formulation for viscoelasticity that includes a level set discretization that switches the  |
|                 |equations from solid to fluid. This option is not currently in production usage.                            |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the polymer stress formulation to EVSS_F:

::

   Polymer Stress Formulation = EVSS_F

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

Guenette, R. and M. Fortin, “A New Mixed Finite Element Method for Computing
Viscoelastic Flow,” J. Non-Newtonian Fluid Mech., 60 (1995) 27-52.

Rajagopalan, D., R. C. Armstrong and R. A. Brown, “Finite Element Methods for
Calculation of Viscoelastic Fluids with a Newtonian Viscosity”, J. Non-Newtonian
Fluid Mech., 36 (1990) 159-192.
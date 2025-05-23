***********************************
Disable Viscosity Sensitivities
***********************************

::

	Disable Viscosity Sensitivities = {yes | no}

-----------------------
Description / Usage
-----------------------

This optional card permits the analyst to omit the sensitivities of a shear-thinning
viscosity model with respect to shear rate from the Jacobian.Valid options for this card
are

yes
    Omit the sensitivities of a shear-thinning viscosity model with respect to
    shear rate from the Jacobian
no
    Form the complete Jacobian.

Currently, this card will have an effect only when using the following viscosity models:
**POWER_LAW, CARREAU, BINGHAM** (see the *Liquid Constitutive Equation*
card).

The default value is **no**.

------------
Examples
------------

The following is a sample card:
::

	Disable Viscosity Sensitivities = yes

-------------------------
Technical Discussion
-------------------------

It has been observed that when these terms are included for very highly
shear-thinning models the result can be non-convergence. In such situations,
disabling these terms can often result in a convergent answer but at
a convergence rate far less than the usual quadratic.


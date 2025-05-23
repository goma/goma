**************
Yield Exponent
**************

::

   Yield Exponent = CONSTANT <float> [t]

-------------------
Description / Usage
-------------------

This card is used to specify the model for the yield exponent parameter, *F*,
for the **BINGHAM** model option of the *Liquid Constitutive Equation* card, or
when the *Polymer Constitutive Equation* card is **SARAMITO_OLDROYDB**,
**SARAMITO_GIESEKUS**, or **SARAMITO_PTT**. Definitions of the input parameters
are as follows:

CONSTANT <float>     
    Name of the model for the yield exponent.                                                                   
    <float> the value of the yield exponent, *F*, which has the dimensions
    of inverse shear-rate in whatever units are consistent with the
    problem of interest and which connotes the steepness of the
    transition from solid to fluid behavior for the Bingham-Carreau-Yasuda
    model or the Saramito yield stress model.                                                                                            

For the **BINGHAM** model, if *F* is large, the material has an abrupt transition from solid-like to 
fluid-like behavior, whereas for a small *F*, the transition is more gradual.

For the **SARAMITO_OLDROYDB**, **SARAMITO_GIESEKUS**, and **SARAMITO_PTT** models, the material has
and abrupt transition when *F* equals zero. This Transistion becomes smooth for nonzero when *F* 
is greater than zero, with the transition becoming more gradual as *F* increases.

--------
Examples
--------

The following is a sample card that sets the yield exponent to 10.0

::

   Yield Exponent = CONSTANT 10.0.

--------------------
Technical Discussion
--------------------

See Description/Usage for *Liquid Constitutive Equation* and *Polymer Constitutive Equation* cards.


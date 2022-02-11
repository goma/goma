********************
Polymer Yield Stress
********************

::

   Polymer Yield Stress = CONSTANT <float> []

-------------------
Description / Usage
-------------------

This card is required when using the Saramito yield model. The card should be
included in the input when the option selected for the *Polymer Constitutive
Equation* card is **SARAMITO_OLDROYDB**, **SARAMITO_GIESEKUS**, or
**SARAMITO_PTT**. Definitions of the input parameters are as follows:

CONSTANT <float>    
    Name of the model for the yield stress.
    * <float> - the value of the yield stress.

This card does not have to be present for constitutive equations other than
**SARAMITO_OLDROYDB**, **SARAMITO_GIESEKUS**, and **SARAMITO_PTT**

--------
Examples
--------

The following is a sample card that sets the polymer yield stress to 12:

::

   Polymer Yield Stress = CONSTANT 12

**********************
Epsilon Regularization
**********************

::

   Epsilon Regularization = CONSTANT <float>

-------------------
Description / Usage
-------------------

This card is used to specify the model for the regularization used in the generalized 
Newtonian viscosity.

Currently used for *HERSCHEL_BULKLEY* models

CONSTANT <float>
   We set the epsilon value to the float specified



--------
Examples
--------

The following is a sample card that sets Papanastasiou regularization to a HERSCHEL_BULKLEY fluid
with an epsilon on the power law term

::

   Liquid Constitutive Equation = HERSCHEL_BULKLEY
   ...
   Regularization Model         = PAPANASTASIOU_EPSILON
   Epsilon Regularization       = CONSTANT 1e-12

Similarly applying only to the Yield term for epsilon regularization

::

   Liquid Constitutive Equation = HERSCHEL_BULKLEY
   ...
   Regularization Model         = EPSILON_YIELD_ONLY
   Epsilon Regularization       = CONSTANT 1e-5

--------------------
Technical Discussion
--------------------

See Description/Usage for *Liquid Constitutive Equation*


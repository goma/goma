********************
Regularization Model
********************

::

   Regularization Model = <MODEL>

-------------------
Description / Usage
-------------------

This card is used to specify the model for the regularization used in the generalized 
Newtonian viscosity.

Currently used for *HERSCHEL_BULKLEY* models

Choices for <MODEL> in *HERSCHEL_BULKLEY* fluids:

EPSILON
   Epsilon regularization applied to both the power law term and the yield term (default for HERSCHEL_BULKLEY)

PAPANASTASIOU
   Papanastasiou regularization applied to the yield parameter (default for HERSCHEL_BULKLEY_PAPANASTASIOU)

EPSILON_YIELD_ONLY
   Epsilon regularization applied only to the yield term

PAPANASTASIOU_EPSILON
   Papanastasiou regularization applied to the yield term and an epsilon on the power law term



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


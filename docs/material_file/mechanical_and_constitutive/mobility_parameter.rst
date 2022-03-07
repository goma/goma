**********************
**Mobility Parameter**
**********************

::

   Mobility Parameter = CONSTANT <float> []

-----------------------
**Description / Usage**
-----------------------

This card is used in the Giesekus model in the nonlinear stress terms. The card should
be included in the input when the option selected for the *Polymer Constitutive Equation*
card is **GIESEKUS**. Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the mobility parameter.                                                               |
|                 |                                                                                                            |
|                 | * <float> - the value of the mobility parameter.                                                           |
+-----------------+------------------------------------------------------------------------------------------------------------+

This card does not have to be present for constitutive equations other than
**GIESEKUS**.

------------
**Examples**
------------

The following is a sample card that sets the mobility parameter to 0.2:

::

   Mobility Parameter = CONSTANT 0.2

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




********************
**PTT Xi parameter**
********************

::

   PTT Xi parameter = {model_name} <float>

-----------------------
**Description / Usage**
-----------------------

This card is used in the Phan-Thien Tanner model in the nonlinear stress terms. The
card should be included in the input when the option selected for the *Polymer
Constitutive Equation* card is **PTT**. Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for PTT Xi parameter.                                                                     |
|                 |                                                                                                            |
|                 | * <float> - the value of the PTT Xi parameter.                                                             |
+-----------------+------------------------------------------------------------------------------------------------------------+

This card does not have to be present for constitutive equations other than **PTT**.

------------
**Examples**
------------

The following is a sample card that sets the PTT Xi parameter to 0.1:

::

   PTT Xi parameter = CONSTANT 0.10

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




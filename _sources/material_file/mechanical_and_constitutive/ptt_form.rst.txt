********
PTT Form
********

::

   PTT Form = {model_name}

-----------------------
**Description / Usage**
-----------------------

This card is used in the Phan-Thien Tanner model in the nonlinear stress terms.
This is an optional card, and the default behavior is preserved as *EXPONENTIAL* for
backwards compatibility.
Definitions of the input parameters are as follows:

EXPONENTIAL
   The default model, the exponential PTT form
LINEAR
   Linear form of the PTT constitutive equation.

------------
**Examples**
------------

The following is a sample card that sets the PTT to the Linear form

::

   PTT Form = LINEAR

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




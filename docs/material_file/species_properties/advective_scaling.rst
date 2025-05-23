*********************
**Advective Scaling**
*********************

::

   Advective Scaling = {model_name} <species> <float>

-----------------------
**Description / Usage**
-----------------------

This material property card permits the user to scale only the advective terms in one or
more of the species transport equations by a fixed constant. This may be useful when
solving problems with non-standard concentrations or for stability reasons.

A single {model_name} is available; it and its parameters are described below:

+-----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**           |Model used to specify the advective scaling.                                         |
|                       |                                                                                     |
|                       | * <species> - the index of the species equation to which the advective scaling will |
|                       |   occur.                                                                            |
|                       | * <float> - scaling, the actual value for the multiplicative scaling factor.        |
+-----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Here is an example of the card:

::

   Advective Scaling = CONSTANT 0 0.0

In this case, the card is being used to eliminate the advective terms in the conservation
equation for species 0.

-------------------------
**Technical Discussion**
-------------------------

The advective terms in the species conservation equations take the form, u ⋅ ∇c where
c is the species concentration and u the fluid velocity.



--------------
**References**
--------------

No References.
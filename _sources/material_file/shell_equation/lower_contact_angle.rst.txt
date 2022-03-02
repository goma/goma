***********************
**Lower Contact Angle**
***********************

::

   Lower Contact Angle = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card sets contact angle of the liquid phase on the lower-wall for the two-phase
capability in the lub_p equation (viz. when using the level-set equation to model the
motion of a meniscus in a thin gap, where the in-plan curvature is neglected. Currently
one model {model_name} is permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model is used to set a constant contact able of the the free surface at the     |
|                          |lower wall. Contact angle of less than 90 degrees is considered as nonwetting        |
|                          |with respect to the heavier level-set phase. Only one floating point value is        |
|                          |required.                                                                            |
|                          |                                                                                     |
|                          | * <float1> is the contact angle in degrees.                                         |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Lower Contact Angle = CONSTANT 180.

This card results in an lower wall contact able to 180 degrees, which is perfectly
wetting. If the lower wall is given the same angle, then the capillary pressure jump will
go as 2/h, where h is the gap.





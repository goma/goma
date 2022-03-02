************************************
**Second Level Set Momentum Source**
************************************

::

   Second Level Set Momentum Source = {model_name} {float_list} {char_string} [M/Lt]

-----------------------
**Description / Usage**
-----------------------

This card allows to the user to specify a second thermal Navier-Stokes volumetric
momentum source model that will be applied to one side of a level set interfacial curve:

+-----------------+------------------------------------------------------------------------------------------------------------+
|{model_name}     |The name of the momentum source model can only be **CONSTANT** at the current time.                         |
+-----------------+------------------------------------------------------------------------------------------------------------+
|{float1}         |This is a single float parameter which is the value of volumetric momentum source term applied to the       |
|                 |second level set phase fluid.[F/L3]                                                                         |
+-----------------+------------------------------------------------------------------------------------------------------------+
|{char_string}    |This string may take the values POSITIVE or NEGATIVE. It identifies which side of the interface the         |
|                 |preceding momentum source model is applied to.                                                              |
+-----------------+------------------------------------------------------------------------------------------------------------+

This card allows the user to apply one of the several momentum source models
implement in Goma to one side of the interface while the other side recieves the
constant momentum source value listed on this card. The side of the interface that
corresponds to char_string appearing on this card receives the constant momentum
source value. The opposite side’s momentum source is determined from the other,
(possibly) more complex model. Transition between them is accomplished using
smooth Heaviside functions whose width is given on the Level Set Length Scale card.
Note that it is the prescence of the this card in the material file that actually activates
this selection process.

------------
**Examples**
------------

The following is a usage example for this card:

::

   Navier-Stokes Source = SUSPEND 0. 0. -980.0 1.34e3

::

   Second Level Set Momentum Source = CONSTANT. 1.0e-4 POSITIVE

This setup will cause the negattive side of the interface to receive momentum source
values obtained from the USER model with the parameters listed above . The positive
side of the interface will show a constant momentum source of 1.0e-4.

-------------------------
**Technical Discussion**
-------------------------

An important thing to note is that the units of the quantity specified on this card are
units of force per volume in exact correspondence to the units used with the preceding
momentum source model. Note also that this card should not be used when using the
LEVEL_SET momentum source model. For one thing, it makes no sense and for
another thing the values specified on the latter model are simply the gravitational
acceleration and therefore are inconsistent with this card.




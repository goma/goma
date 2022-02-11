**********************************
**Second Level Set Heat Capacity**
**********************************

::

   Second Level Set Heat Capacity = {model_name} {float_list} {char_string} [M/Lt]

-----------------------
**Description / Usage**
-----------------------

This card allows to the user to specify a second thermal heat capacity model that will be
applied to one side of a level set interfacial curve:

+-----------------+------------------------------------------------------------------------------------------------------------+
|{model_name}     |The name of the heat capacity model can only be **CONSTANT** at the current time.                           |
+-----------------+------------------------------------------------------------------------------------------------------------+
|{float1}         |This is a single float parameter which is the value of heat capacity applied to the second level set phase  |
|                 |fluid.                                                                                                      |
+-----------------+------------------------------------------------------------------------------------------------------------+
|{char_string}    |This string may take the values POSITIVE or NEGATIVE. It identifies which side of the interface the         |
|                 |preceding heat capacity model is applied to.                                                                |
+-----------------+------------------------------------------------------------------------------------------------------------+

This card allows the user to apply a CONSTANT or USER model to one side of the
interface while the other side recieves the constant heat capacity value listed on this
card. The side of the interface that corresponds to char_string appearing on this card
receives the constant heat capacity value. The opposite side’s heat capacity is
determined from the other, (possibly) more complex model. Transition between them
is accomplished using smooth Heaviside functions whose width is given on the Level
Set Length Scale card. Note that it is the prescence of the this card in the material file
that actually activates this selection process.

------------
**Examples**
------------

The following is a usage example for this card:

::

   Heat Capacity = ENTHALPY 1.e4 0.1

::

   Second Level Set Heat Capacity = CONSTANT. 1.0 POSITIVE

This setup will cause the negattive side of the interface to receive heat capacity values
obtained from the USER model with the parameters listed above . The positive side of
the interface will show a constant heat capacity of 1.0e-4.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




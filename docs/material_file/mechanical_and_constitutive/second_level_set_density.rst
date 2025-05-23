****************************
**Second Level Set Density**
****************************

::

   Second Level Set Density = {model_name} {float_list} {char_string} [M/Lt]

-----------------------
**Description / Usage**
-----------------------

This card allows to the user to specify a second density model that will be applied to
one side of a level set interfacial curve:

+-----------------+------------------------------------------------------------------------------------------------------------+
|{model_name}     |The name of the density model can only be **CONSTANT** at the current time.                                 |
+-----------------+------------------------------------------------------------------------------------------------------------+
|{float1}         |This is a single float parameter which is the value of density applied to the second level set phase fluid. |
+-----------------+------------------------------------------------------------------------------------------------------------+
|{char_string}    |This string may take the values POSITIVE or NEGATIVE. It identifies which side of the interface the         |
|                 |preceding density model is applied to.                                                                      |
+-----------------+------------------------------------------------------------------------------------------------------------+

This card allows the user to apply one of the several complex density models currently
available in Goma to one side of the interface while the other side recieves the constant
density value listed on this card. The side of the interface that corresponds to
char_string appearing on this card recieves the constant density value. The opposite
sides density is determined from the other, more complex model. Transition between
them is accomplished using smooth Heaviside functions whose width is given on the
Level Set Length Scale card. Note that it is the prescence of the this card in the
material file that actually activates this selection process.

------------
**Examples**
------------

The following is a usage example for this card:

::

   Density = SUSPENSION 1.0 1.0 1.0

::

   Second Level Set Density = CONSTANT. 1.0 POSITIVE

This setup will cause the negattive side of the interface to receive density values
obtained from the SUSPENSION model with the parameters listed above . The
positive side of the interface will show a constant density of 1.0.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




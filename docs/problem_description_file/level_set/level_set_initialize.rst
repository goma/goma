************************
Level Set Initialize
************************

::

	Level Set Initialize = <char_string> <float1> <float2>

-----------------------
Description / Usage
-----------------------

This card is used to initialize fields around the zero level set.

<Char_string>
    A character string which identifies dependent variable to be initialized.
    It is taken from the list of names on the Initialize card.

<float1>
    Value of the variable on the negative side of the zero level set.

<float2>
    Value of the field on the positive side of the zero level set.

------------
Examples
------------

Two examples of initialization methods are provide below:
::

	Level Set Initialize = TEMPERATURE 0. 100.

-------------------------
Technical Discussion
-------------------------

Not clear whether this capability has been used and tested much. (12/3/2012)


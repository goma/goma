**********
END OF AC
**********

::

    END OF AC

-----------------------
Description / Usage
-----------------------

This card is the companion card to the Number of augmenting conditions card. The 
END OF AC card signals the end of the Augmenting Conditions Specifications section 
of the Goma input. When the Number of augmenting conditions (= -1) is set to negative 
one, Goma will read all the augmenting condition cards until this card is encountered in 
the input file. This card may omitted if the integer N on the Number of augmenting 
conditions card is not -1.

------------
Examples
------------

Typical usage of this card is illustrated below:

::

    Number of augmenting conditions = -1
    .
    .
    .
    END OF AC

-------------------------
Technical Discussion
-------------------------

See companion card Number of augmenting conditions.

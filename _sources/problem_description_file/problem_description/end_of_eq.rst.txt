*************
**END OF EQ**
*************

::

	END OF EQ

-----------------------
**Description / Usage**
-----------------------

This card specifies the end of the list of equations in a material section of the *Problem
Description File*. It is only used when automatic equation counting is used, as
described and activated in the *Number of EQ* card. If the value of <integer> in that card
is set to -1, all EQ cards below this card are ignored, and *Goma* counts the number of
EQ cards between the *Number of EQ* card and the *END OF EQ* card.

Note that the *END of EQ* card should appear in every material section for which
automatic equation counting is being used.

------------
**Examples**
------------

There are no input parameters for this card, which always appears as follows:
::

   END OF EQ

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




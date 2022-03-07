**************
**END OF MAT**
**************

::

	END OF MAT

-----------------------
**Description / Usage**
-----------------------

This card specifies the end of the list of materials. It is only used when automatic
material counting is used, as described and activated in the *Number of Materials* card.
If the value of <integer> in the *Number of Materials* card is set to -1, all MAT cards
below the *END OF MAT* card are ignored, and *Goma* counts the number of MAT cards
between these two cards.

------------
**Examples**
------------

There are no input parameters for this card, which always appears as follows:
::

   END OF MAT

-------------------------
**Technical Discussion**
-------------------------

No Discussion.





~~~~~~~~~~~
END OF BC
~~~~~~~~~~~

::

	END OF BC

***********************
Description / Usage
***********************


This card specifies the end of the list of boundary conditions (BCs), and is only used when automatic BC counting is used, as described in the *Number of BC* card. If the value of <integer> in that card is set to -1, all BC cards below the *END of BC* card are ignored, and *Goma* counts the number of BC cards between the *Number of BC* card and the *END of BC* card.

***********************
Examples
***********************

There are no input parameters for this card, which always appears as follows:
::

   END OF BC

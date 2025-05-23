~~~~~~~~~~~~~~~~
Number of BC
~~~~~~~~~~~~~~~~

::

	Number of BC = <integer>

***********************
Description / Usage
***********************

This required card indicates how many boundary condition (BC) cards are contained in
the *Problem Description File*. The single input parameter is defined as

<integer>
    The number of BC cards that follow.

If <integer> is set to -1, *Goma* will automatically count the number of BC cards
between the *Number of BC* card and the *End of BC* card. This latter usage is generally
preferred if a large number of BCs are to be specified.

***********************
Examples
***********************

Following is a sample card, indicating that there are two BC cards that follow this card.
::

	Number of BC = 2

***********************
Technical Discussion
***********************

If there are more BC cards listed in an input deck than specified on this card, *Goma*
ignores the extras; in other words, only the first <integer> cards are read by *Goma*. If the number of BCs is fewer than the amount specified by <integer>, *Goma* will stop
with an error.

Also note, that if more than one BC on the same variable is specified, only the last one
is applied.


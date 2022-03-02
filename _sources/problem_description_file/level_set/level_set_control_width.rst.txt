***************************
Level Set Control Width
***************************

::

	Level Set Control Width = <float>

-----------------------
Description / Usage
-----------------------

This card is a multiplier on the *Level Set Length Scale* to determine the size of the
region around the zero level set contour over which the level set gradient is averaged.
The value of this parameter defaults to 1.0 if this card is not included.

------------
Examples
------------

This sample card sets the control width to be equivalent to the length scale:
::

	Level Set Control Width = 0.5

-------------------------
Technical Discussion
-------------------------

As noted in the description of the *Level Set Renormalization Tolerance* card,
renormalization is triggered when the average of the level set gradient magnitude has
departed sufficiently from unity. The region over which this average is obtained is
approximately a narrow fixed-width strip on either side of the zero level set contour.
The width of this strip is twice the *Level Set Length Scale* multiplied by the float value
supplied on this card.


--------
FAQs
--------

Usually it is best practice to leave this parameter at its default setting and control the
frequency of renormalization with the renormalization tolerance.


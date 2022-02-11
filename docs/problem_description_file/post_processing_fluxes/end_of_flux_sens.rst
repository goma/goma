********************
**END OF FLUX_SENS**
********************

::

   END OF FLUX_SENS

-----------------------
**Description / Usage**
-----------------------

This card is used to denote the end of a set of *FLUX_SENS* cards and is only used when
the *Post Processing Flux Sensitivities* card is present and one or more *FLUX_SENS*
cards are specified.

------------
**Examples**
------------

A simple example of using this card in context is shown below.
::

   Post Processing Flux Sensitivities =
   FLUX_SENS = VOLUME_FLUX 1 1 0 BC 5 3 flux_sens.out
   END OF FLUX_SENS

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
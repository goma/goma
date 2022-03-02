*******************
**Heat Flux Model**
*******************

::

   Heat Flux Model = USER

-----------------------
**Description / Usage**
-----------------------

**NOT TESTED**. Use this optional card to specify a user-defined model for the
conductive heat flux. The routine “usr_heat_flux” in file user_mp.c must
appropriately define the heat flux/temperature gradient model. The single input
parameter has only one possible value:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**USER**         |the user-defined model for the conductive heat flux.                                                        |
+-----------------+------------------------------------------------------------------------------------------------------------+

If this card is missing or has a different keyword, the Fourier conductive heat flux
model will be used.

------------
**Examples**
------------

Following is the only permissible specification for the card:

::

   Heat Flux Model = USER

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
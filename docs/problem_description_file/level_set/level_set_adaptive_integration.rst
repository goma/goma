**********************************
Level Set Adaptive Integration
**********************************

::

	Level Set Adaptive Integration = {ON | YES | OFF | NO}

-----------------------
Description / Usage
-----------------------

To be used with Subelement integration to improve integration accuracy. Does not
work with subgrid integration or basic level-set. Requires a sharp interface, viz. levelset
length scale of zero. Please see usage nodes below.

{ON | YES}
    Use adaptive integration on surface level set capillary term.
    
{OFF | NO}
    Don’t use adaptive integration.

--------------
**References**
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer

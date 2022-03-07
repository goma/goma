************************************
Level Set Subelement Integration
************************************

::

	Level Set Subelement Integration = {ON | YES | OFF | NO}

-----------------------
Description / Usage
-----------------------

Subelement integration is used to improve integration accuracy for all functions which
invoke a sharp level-set interface. Note here that the Level Set Length Scale
option must be zero. This is possible because the subelement integration scheme
actually produces a geometric representation of the zero level set surface on which
exact line integrals of the surface tension source term term can be peformed. Please
see usage nodes below.

{ON | YES}
    Use subelement integration on surface level set capillary term.

{OFF | NO}
    Don’t use subelement integration.

------------
Examples
------------

This example invokes the subelement integration
::

	Level Set Subelement Integration = ON

-------------------------
Technical Discussion
-------------------------

* **NOTE**: Level Set Length Scale must be set to zero.

* Because of the construction of an in-element interface meshing to find this
  representation, subelement integration cannot be used currently for three
  dimensional problems. Subgrid integration can be, however, but it is inefficient.

* Best to use this integration approach with the property specification method of
  “Second Level-Set “property_name”, e.g. Second Level Set Density, etc.

* Typically this capability greatly improves mass conservation and avoids parasitics
  for surface tension dominated problems.

* **NOTE** that the Level Set Renormalization method must be set to Huygens.

--------------
References
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer

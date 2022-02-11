*************************
**Material is Nondilute**
*************************

::

	Material is nondilute = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This card is a optional for each material section in the *Problem Description File*. It is
used to specify the number of species equations in a phase. The single string parameter
is a boolean, yes or no.

========= ====================================================================
**yes**   the number of species equations is set equal to one less than
          the number of species.
**no**    the number of species equations is set equal to the number
          of species.
========= ====================================================================

When the number of species is equal to the number of species equations, there is an
implied additional species, i.e., the solute, which is not part of species loops, which fills out the specification of the phase.

------------
**Examples**
------------

Following is a sample card:
::

   Material is nondilute = yes

-------------------------
**Technical Discussion**
-------------------------

See the discussion for the “*Number of bulk species*” card.



--------------
**References**
--------------

No References.
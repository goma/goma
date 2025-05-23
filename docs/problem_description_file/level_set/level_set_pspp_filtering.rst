****************************
Level Set PSPP Filtering
****************************

::

	Level Set PSPP filtering = <YES | NO>

-----------------------
Description / Usage
-----------------------

On this card, the user specifies a single char_string.

<YES | ON>
    This string turns on level set PSPP filtering if it is “yes” or “on”.

------------
Examples
------------

A typical PSPP filtering input card looks like:
::

	Level Set PSPP filtering = yes

-------------------------
Technical Discussion
-------------------------

Not entirely clear what this card does, but in the vicinity of the level-set interface, the
Bochev PSPP stabilization scheme is altered. This is recommended when this pressure
stabilization scheme is deployed. See the Pressure Stabilization card.


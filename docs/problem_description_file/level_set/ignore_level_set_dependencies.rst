*********************************
Ignore Level Set Dependencies
*********************************

::

	Ignore Level Set Dependencies = {yes | no}

-----------------------
Description / Usage
-----------------------

Including this card in your input deck with the string parameter set to “yes” instructs
*Goma* to discard the sensitivities of all equations to the level set variable when
constructing the Jacobian matrix. This may have benefits when it comes to stability
and convergence; although, the effectiveness of this card is very much case by case.
Note also that use of this card is consistent only with **Fill Weight Function = Explicit.**
Any other choice will result in an error.

------------
Examples
------------

A sample input card is:
::

	Ignore Level Set Dependencies = yes


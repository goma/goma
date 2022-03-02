*******************************************
Force Initial Level Set Renormalization
*******************************************

::

	Force Initial Level Set Renormalization = <char_string>

-----------------------
Description / Usage
-----------------------

This card is used to invoke a renormalization step prior to the first time step of any
transient computation.

<char_string>
    ``YES|ON`` (not case sensitive) will cause the renormalization procedure to
    occur on the first step. If this card is not included or some other string
    is used here a renormalization will automatically occur on the first time
    step.

------------
Examples
------------

A typical length scale input card looks like:
::

	Force Initial Level Set Renormalization = yes

-------------------------
Technical Discussion
-------------------------

Restarts occur fairly frequently during level set computations. It has been discovered
that the robustness of the subsequent computation can be improved by quite a bit if the
level set field is renormalized at the start of the restart, regardless of the current average
gradient norm error. This card is employed to invoke a renormalization at the start of
any computation, that is, a renormalization procedure is conducted prior to the initial
time step if this card is present in the input deck. It has become standard operating
procedure that when a level set computation runs into computational difficulty the first
step in recovery should be to restart with a forced initial renormalization using this
card.


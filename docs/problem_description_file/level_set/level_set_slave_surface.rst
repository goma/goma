***************************
Level Set Slave Surface
***************************

::

	Level Set Slave Surface = {yes|no}

-----------------------
Description / Usage
-----------------------

This card specifies whether the level set distance function is constrained during the
calculation or evolves with the typical advection equation. Permissible values for this
option are:

yes|on
    The surface is constrained to remain on the initial surfaces throughout the
    calculation (moving with these surfaces if they are moving).

no|off
    The surface evolves normally according to the local velocity field; this is
    the default.

------------
Examples
------------

This is a sample card:
::

	Level Set Slave Surface = on

-------------------------
Technical Discussion
-------------------------

In a typical level set simulation, the surface is first initialized with the *Level Set
Initialization Method* card, and then the surface evolves in time according to the local
velocity field. Using this card, however, the surface is constrained to remain on the
initial surfaces. If the initial surfaces are static, then the level set surface remains stationary. For moving interfaces such as those defined by an isosurface or a side set,
the level set function is reinitialized at each Newton iteration to match the moving
surface.


--------------
References
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer

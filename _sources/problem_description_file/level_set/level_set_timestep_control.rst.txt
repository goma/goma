******************************
Level Set Timestep Control
******************************

::

	Level Set Timestep Control = <YES | NO>

-----------------------
Description / Usage
-----------------------

On this card, the user specifies a single char_string.

<YES | ON>
    This string turns on level set timestep control if it is “yes” or “on”.

------------
Examples
------------

A typical length scale input card looks like:
::

	Level Set Timestep Control = yes

-------------------------
Technical Discussion
-------------------------

In normal operations, the error norm of the level set function is not included in
controlling the size of the time step decided upon by the variable timestep size
integrator. Inclusion of this card will add the level set unknown to the list of update
error norms used to decide the time step size. In other words, use this card when you
want the changes of the level set function to affect the timestep size. If this card is not
used, the default behavior is to ignore the level set degrees of freedom in controlling
the timestep size.


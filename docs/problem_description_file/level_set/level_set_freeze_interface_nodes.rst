********************************
Level Set Freeze Interface Nodes
********************************

::

	Level Set Freeze Interface Nodes = {yes|no}

-----------------------
Description / Usage
-----------------------

This card specifies whether the interface nodes of the level set are affected by
the renormalization method. Permissible values for this option are:

yes|on
    The interface nodes are frozen during the renormalization process, meaning
    that their level set values are not changed.

no|off
    The interface nodes are not frozen, and their level set values may be
    changed during the renormalization process; this is the default.

------------
Examples
------------

This is a sample card:
::

	Level Set Freeze Interface Nodes = yes

-------------------------
Technical Discussion
-------------------------

This option is useful to avoid unintended changes to the level set values at the
interface nodes during the renormalization process. 

--------------
References
--------------

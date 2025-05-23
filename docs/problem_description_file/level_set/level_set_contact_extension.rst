*******************************
Level Set Contact Extension
*******************************

::

	Level Set Contact Extension = {yes|no}

-----------------------
Description / Usage
-----------------------

This card specifies whether the level set surface is considered to extend into boundaries
when performing renormalization of the level set distance function. This card applies
only if *Level Set Renormalization Method* = **Huygens_Constrained.** Permissible
values for this option are:

yes|on
    The level set interface is considered to extend smoothly into the
    boundaries.

no|off
    The level set interface ends at the boundaries; this is the default.

------------
Examples
------------

This is a sample input card:
::

	Level Set Contact Extension = no

-------------------------
Technical Discussion
-------------------------

When renormalizing the level set distance function, the behavior of the interface near
boundaries is important. When the interface is considered to end at the boundary, a
large number of grid points may be closest to this boundary point. This appears as a
cusp in the interface and can make it difficult to achieve sharp contact angles because
of the very large capillary force that results. One method to alleviate this is to extend
the interface smoothly into the boundaries to eliminate the cusp in the interface. The
current algorithm, however, can cause errors when employed near corners of the domain. Until this is resolved, this option can only be recommended for domains
without interior corners.


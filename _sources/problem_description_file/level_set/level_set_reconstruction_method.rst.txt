***********************************
Level Set Reconstruction Method
***********************************

::

	Level Set Reconstruction Method = {char_string}

-----------------------
Description / Usage
-----------------------

This card indicates the method used to perform the Huygens renormalization of the
level set function. This card applies only if *Level Set Renormalization Method* is set to
**Huygens** or **Huygens_Constrained.** Permissible values of {char_string} are:

POINTS
    A list of points on the interface is formed and the renormalized distance
    is computed as the distance to the nearest point in this list; this is the
    default method.

FACETS
    A list of connected facets on the interface is formed and the renormalized
    distance is computed as the distance to the nearest point on the nearest
    facet in this list. Currently this option is not supported for
    3-dimensional calculations.

------------
Examples
------------

This is a sample input card:
::

	Level Set Reconstruction Method = FACETS

-------------------------
Technical Discussion
-------------------------

As described for the *Level Set Renormalization Method* card, Huygens based
renormalization is performed by reconstructing the level set surface and computing the
distance to the nearest point on this surface. Here, the method of reconstructing the
level set surface is addressed. Either a set of points on the interface is formed or a
connected set of facets is formed. The advantage to using connected facets is that the
interface is better described between the points on the interface. However, the
calculation of the faceted geometry is slightly more expensive computationally. Also,
the current implementation is limited to 2-dimensional simulations.


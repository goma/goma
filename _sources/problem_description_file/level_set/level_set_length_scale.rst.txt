**************************
Level Set Length Scale
**************************

::

	Level Set Length Scale = <float>

-----------------------
Description / Usage
-----------------------

On this card, the user specifies a single char_string.

<float>
    This value represents the size of the region around the zero level set
    function contour in which interfacial physical quantities, for example,
    surface tension, will be present.

Stability and conservation of phase volume are dependent upon this value to a
significant degree. Experimentation has revealed that this float value should be
between two and three times the average linear dimension of the elements in the mesh.

------------
Examples
------------

A typical length scale input card looks like:
::

	Level Set Length Scale = 0.3

-------------------------
Technical Discussion
-------------------------

The level set method is an *embedded* interface method. That is, the location of the
interface is not known explicitly as a geometric parameter of the problem, but rather it
is abstracted as a level contour of a higher dimensional function. This is convenient in
many ways, but it does mean that phenomena associated with the interface, for
example, surface tension, must enter the problem spread over a region near the zero
level set contour. The *Level Set Length Scale* sets the size of this region.

A good example of the application of the *Level Set Length Scale* parameter is in how
surface tension is included in problems using level set interface tracking. The following
tensor is added to the fluid momentum equation:

.. math:: 

   \underline{\underline{T}} = \sigma \delta_{\alpha} \left( F \right) \left( \underline{\underline{I}} - \underline{n} \,\underline{n} \right)

where :math:`F` is the level set function itself, :math:`\underline{n}` is the unit normal to the level
set contour, :math:`\underline{\underline{I}}` is the unit tensor, :math:`\sigma` the surface tension, and :math:`\delta_{\alpha} \left( F \right)` is a “smooth” Dirac
function given by:

.. math::

   \delta_{\alpha} \left( F \right) = \lvert \nabla F \rvert \frac{1 + \cos \left( \pi F/\alpha \right)}{2 \alpha}, \quad \lvert F \rvert \leq \alpha

In this example, the parameter :math:`\alpha` would be equal to one-half the *Level Set Length Scale*
value specified on this card.


--------
FAQs
--------

How should the Length Scale value be chosen? Trial and error is often the best method
to determine an appropriate value for this parameter. However, experience has shown
that values for Level Set Length Scale that are between two and three times the average
element linear dimension seem to work best.

--------------
References
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer

..
	 TODO - There are two equation pictures that I put in for someone to write. There is also a equation in line 53 between the comas that needs to be inserted. 

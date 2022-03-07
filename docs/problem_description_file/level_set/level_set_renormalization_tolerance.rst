***************************************
Level Set Renormalization Tolerance
***************************************

::

	Level Set Renormalization Tolerance = <float>

-----------------------
Description / Usage
-----------------------

This parameter provides a means for controlling how often renormalization
(redistancing) operations are performed on the level set function as it evolves by fixing
the size of the deviation allowed between the average absolute magnitude of the level
set function gradient near the level set interface and unity, the theoretical value
observed for a pure distance function.

<float>
    Value of the tolerance, the allowable deviation.

The range of this parameter is any positive real number, however, it is rare to use values
smaller than 0.1 or larger than 5.0. The value of the tolerance defaults to 0.5 if this card
is not specified.

------------
Examples
------------

This is a sample renormalization card:
::

	Level Set Renormalization Tolerance = 0.05

-------------------------
Technical Discussion
-------------------------

One of the key properties of the level set function is that it is a smooth function near to
the interface. In particular, if the level set function is a distance function then the
magnitude of its gradient on the zero level contour should always be unity. This fact is
used to provide a criterion for invoking a renormalization procedure. The gradient of
the level set function is found within a fixed region around the zero level set contour
(see *Level Set Control Width*). The integrated average of the magnitude of this vector is
determined and compare to unity. Should this difference differ by greater than the value
for Renormalization Tolerance identified on this card, a renormalization procedure will
presently be initiated.


--------
FAQs
--------

What is a proper value for this parameter? Values on the order of unity should work
well. Renormalization based on gradient can be disabled completely by choosing a very large value for this parameter. Conversely, a very small value will always result in
a renormalization step.

Is it possible to renormalize too often? Yes. Renormalization is an extraphysical
procedure designed solely to improve the numerical performance of the interface
tracker. As such, it can add or subtract volume to or from the phases represented by the
interface contour. Renormalizing too often, therefore, can result in errors being
introduced. The renormalization procedure, Huygens_Constrained, attempts to
mitigate this effect.


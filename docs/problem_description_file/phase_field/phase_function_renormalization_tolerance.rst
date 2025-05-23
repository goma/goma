********************************************
**Phase Function Renormalization Tolerance**
********************************************

::

	Phase Funtion Renormalization Tolerance = <float>

-----------------------
**Description / Usage**
-----------------------

This parameter provides a means for controlling how often renormalization
(redistancing) operations are performed on the phase function fields as they evolve by
fixing the size of the deviation allowed between the average absolute magnitude of the phase function gradient near each respecitve interface and unity, the theoretical value
observed for a pure distance function.

<float>
    Value of the tolerance, the allowable deviation.

The range of this parameter is any positive real number, however, it is rare to use values
smaller than 0.1 or larger than 5.0. The value of the tolerance defaults to 0.5 if this card
is not specified. Note that a global parameter value is applied to all phase function
fields in the problem. Currently, there is no provision for each phase function field to
have a unique value for this parameter.

This parameter is exactly analogous to the similarly named parameter used in standard
level set interface tracking.

------------
**Examples**
------------

This is a sample renormalization card:
::

	Phase Function Renormalization Tolerance = 0.25

-------------------------
**Technical Discussion**
-------------------------

The reader is referred to the Technical Discussion associated with Level Set
Renormalization Tolerance card as it is virtually identical to the operation of it in the
current context. The only thing to note is that each phase function is evaluted
separately against this tolerance and each function is renormalized independently if the
tolerance is exceeded. That is, exceeding the tolerance by one phase function field only
triggers renormalization for that field. The other phase function fields are left
unaltered.


--------
**FAQs**
--------

What is a proper value for this parameter? Values on the order of unity should work
well. Renormalization based on gradient can be disabled completely by choosing a
very large value for this parameter. Conversely, a very small value will always result in
a renormalization step.

Is it possible to renormalize too often? Yes. Renormalization is an extraphysical
procedure designed solely to improve the numerical performance of the interface
tracker. As such, it can add or subtract volume to or from the phases represented by the
interface contour. Renormalizing too often, therefore, can result in errors being introduced. The renormalization procedure, Huygens_Constrained, attempts to
mitigate this effect.


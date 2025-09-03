Eigen Wave Numbers
==================

.. code-block:: none

    Eigen Wave Numbers = <float1> ... <floatN>

Description/Usage
-----------------

If the type of stability requested in the Linear Stability card indicates that 3D stability 
of a 2D base flow is to be performed, then this card is required. Otherwise it is optional 
(and ignored). This card is used to specify the wave numbers in the normal (z) direction 
for 3D stability of a 2D flow. The valid inputs are:

w1, w2, ..., wN
    Any number (subject to input line length restrictions) of real numbers that are used as normal wave numbers.

There is no default value.

Examples
--------

Here is a sample card:

::

    Eigen Wave Numbers = 0.0 0.01 0.1 1.0 2.0 3.0

Technical Discussion
--------------------

Each wi is used in turn as the wave number in the 3rd direction (the normal direction) 
for which linear stability analysis will be performed. They are equivalent to a 
wavelength of 2π/wi. A wavenumber of wi=0 corresponds to the regular 2D stability of 
a 2D base flow.

This card is applicable to both eggroll and ARPACK eigensolvers. However, 3D of 2D 
stability can only be done in parallel when ARPACK is used.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

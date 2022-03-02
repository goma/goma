*******************
**Element Mapping**
*******************

::

	Element Mapping = <char_string>

-----------------------
**Description / Usage**
-----------------------

This card allows the user to set the order of the finite element shape mapping between the canonical element and each physical element. Valid options for {char_string} are:

================= =========================================================
**isoparametric** This choice sets the element order mapping to the 
                  highest order present in the problem. *However, if a mesh
                  displacement field is present, the element mapping 
                  order is the interpolation order of the mesh 
                  displacement field*.
**Q1**            This choice sets the element mapping order to bilinear.
**Q2**            This choice sets the element mapping order to biquadratic.
**SP**            This choice sets the element mapping to order to
                  subparametric.
================= =========================================================

------------
**Examples**
------------

Some text like this:
::

   Element Mapping = isoparametric

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
*********************
**Polymer Weighting**
*********************

::

   Polymer Weighting = <float> [t/L]

-----------------------
**Description / Usage**
-----------------------

This card is only used if the value of the *Polymer Weight Function* card is **SUPG**. The
single input parameter is defined as

<float> - scale factor for the upwind term in the Petrov-Galerkin
formulation. If this is set to zero, a Galerkin weight
function is used. The correct scaling for this term is
the inverse of the average inflow velocity.

------------
**Examples**
------------

The following is a companion pair of sample input cards that includes setting the
polymer weighting to 0.1:

::

   Polymer Weight Function = SUPG

::

   Polymer Weighting = 0.1

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
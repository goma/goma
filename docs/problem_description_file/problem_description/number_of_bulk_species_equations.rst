************************************
**Number of Bulk Species Equations**
************************************

::

	Number of bulk species equations = <integer>

-----------------------
**Description / Usage**
-----------------------

This card is optional but strongly recommended for each material section in the
*Problem Description File*. It is used to specify the number of species equations in a
phase. The word, bulk, here, refers to its being distributed throughout the domain, not
just at a surface. The single input parameter is defined as:

========= ====================================================================
<integer> The number of species equations; if the value of <integer> is 0, 
          then no conservation equations for species are solved for.
========= ====================================================================

When the number of species is equal to the number of species equations, there is an
implied additional species, i.e., the solute, which is not part of species loops, which fills out the specification of the phase.

------------
**Examples**
------------

Following is a sample card:
::

   Number of bulk species equations = 1

-------------------------
**Technical Discussion**
-------------------------

See the discussion for the “*Number of bulk species*” card.



--------------
**References**
--------------

No References.
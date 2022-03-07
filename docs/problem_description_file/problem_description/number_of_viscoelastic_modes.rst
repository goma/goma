********************************
**Number of Viscoelastic Modes**
********************************

::

	Number of viscoelastic modes = <integer>

-----------------------
**Description / Usage**
-----------------------

This card is required only if you are performing a viscoelastic simulation and have
included stress equations in the equation section and chosen a viscoelastic constitutive
equation in the material file. The integer value denotes how many viscoelastic tensor
stress equations are to be used. The number of modes can vary from a minimum of 1 to
a maximum of 8. The input parameter is defined as

========= ====================================================================
<integer> The number of viscoelastic modes, which must be greater
          than zero, but less than nine.
========= ====================================================================

------------
**Examples**
------------

The following is a sample card, indicating that a calculation with two viscoelastic stress
modes is being undertaken:
::

   Number of viscoelastic modes = 2

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

Please see the viscoelastic tutorial memo for a discussion of multimode viscoelastic
equations:

	GT-014.1: Tutorial for Running Viscoelastic Flow Problems with GOMA, June
	21, 2000, R. R. Rao
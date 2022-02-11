*****************
**Enormsq Field**
*****************

::

   Enormsq Field = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This norm is based on the *ENORM* field variable (which, in turn, is derived from the
*VOLTAGE* field variable).

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the norm.
**no**   Do not calculate the norm.
======== ===============================================

The field is stored in **GENS0**, **GENS1**, and **GENS2** (if 3D) in the output EXODUS II
file.

------------
**Examples**
------------

The following is a sample input card to calculate the field:
::

   Enormsq Field = yes

-------------------------
**Technical Discussion**
-------------------------

This post-processing variable is equal to :math:`\underline{\Delta}` :math:`enorm^2` . This, in turn, should approximate :math:`\underline{\Delta}` :math:`\mid` ( :math:`\underline{\Delta}V`
:math:`\mid^2` ).

See also the *Enormsq Field Norm* post processing option.



--------------
**References**
--------------

No References.
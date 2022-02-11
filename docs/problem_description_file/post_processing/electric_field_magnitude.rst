****************************
**Electric Field Magnitude**
****************************

::

   Electric Field Magnitude = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The magnitude of the electric field is written to the output EXODUS II file. The
electric field is calculated as the negative gradient of the *VOLTAGE* field variable.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the electric field magnitude.
**no**   Do not calculate the electric field magnitude.
======== ===============================================

The electric field magnitude is called EE in the output EXODUS II file.

------------
**Examples**
------------

The following is a sample input card to calculate the Electric Field Magnitude:
::

   Electric Field Magnitude = yes

-------------------------
**Technical Discussion**
-------------------------

See also the *Electric Field* post processing option.



--------------
**References**
--------------

No References.
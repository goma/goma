******************
**Electric Field**
******************

::

   Electric Field = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The electric field vector components are written to the output EXODUS II file. The
electric field is calculated as the negative gradient of the VOLTAGE field variable.

The permissible values for this postprocessing option are

======== ===============================================
**yes**  Calculate the electric field vectors.
**no**   Do not calculate the electric field vectors.
======== ===============================================

The vector components are called **EX**, **EY**, and (for three dimensional problems) **EZ** in the output EXODUS II file.

------------
**Examples**
------------

The following is a sample input card to calculate the Electric Field vector components:
::

   Electric Field = yes

-------------------------
**Technical Discussion**
-------------------------

See also the *Electric Field Magnitude* post processing option.



--------------
**References**
--------------

No References.
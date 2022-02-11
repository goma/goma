********************
**Energy Fluxlines**
********************

::

   Energy Fluxlines= {yes | no}

-----------------------
**Description / Usage**
-----------------------

This post-processing option triggers the energy fluxlines to be calculated. The energy
flux function is analogous to the stream function, its contours representing paths of
energy flow through the domain. This variable is called **TFLUX** in the output
EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate and write the energy fluxlines to the
         output EXODUS II file.
**no**   Do not calculate the energy fluxlines.
======== ===============================================

------------
**Examples**
------------

The energy fluxlines are calculated in this sample input card:
::

   Energy Fluxlines = yes

-------------------------
**Technical Discussion**
-------------------------

The *Energy Conduction Vectors* must also be activated for this post processing option
to work.



--------------
**References**
--------------

No References.
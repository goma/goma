********************
**Time Derivatives**
********************

::

   Time Derivatives = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This option enables writing the time derivative of all the field variables as nodal
variables to the output EXODUS II file. These variables are labeled **XDOT0** (mesh
velocity in x direction), **XDOT1** (mesh velocity in y direction), **XDOT2** (mesh
velocity in z direction), **VDOT0** (fluid acceleration in x direction), **VDOT1** (fluid
acceleration in y direction), **VDOT2** (fluid acceleration in z direction), **TDOT** (rate of
temperature change), **Y0DOT** (rate of 1st species concentration change), **Y1D0T** (rate
of second species concentration change), and so on. The quantities can then be
contoured or displayed by some other means with a visualization or graphics package.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the time derivatives and write them 
         as nodal variables in the output EXODUS II file.
**no**   Do not calculate the time derivatives.
======== ===============================================

------------
**Examples**
------------

The following sample card requests that time derivatives be written to the EXODUS II
file:
::

   Time Derivatives = yes

-------------------------
**Technical Discussion**
-------------------------

Currently, this routine uses the values in the global vector *xdot* to report this data.
During the first time step, all the *xdot* values are zero; by the second time step, these
data should be realistic.



--------------
**References**
--------------

No References.
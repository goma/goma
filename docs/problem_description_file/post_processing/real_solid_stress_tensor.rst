****************************
**Real Solid Stress Tensor**
****************************

::

   Real Solid Stress Tensor = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The real solid stress tensor is associated with the equations of elasticity. If the mesh
motion is of *LAGRANGIAN* type, then these quantities are not available; if is of
*TOTAL_ALE* type, they are available. The stress tensor has six entries (in three
dimensions because it is symmetric) called **T11_RS, T22_RS, T33_RS, T12_RS, T13_RS**, and **T23_RS** in the output EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the real solid stress tensor and write
         to the output EXODUS II file.
**no**   Do not calculate the real solid stress tensor.
======== ===============================================

------------
**Examples**
------------

No stress tensor is written for the following sample input card:
::

   Real Solid Stress Tensor = no

-------------------------
**Technical Discussion**
-------------------------

This option is applicable only to *TOTAL_ALE* mesh motion types (see *Mesh Motion*
card). Compare this with *Mesh Stress Tensor* post processing option for other types of
mesh motion.



--------------
**References**
--------------

No References.
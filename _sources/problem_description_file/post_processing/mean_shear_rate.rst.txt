*******************
**Mean Shear Rate**
*******************

::

	Mean shear rate = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The mean shear rate is defined as

.. figure:: /figures/321_goma_physics.png
	:align: center
	:width: 90%

where :math:`II_D` is the second invariant of **D**, the strain-rate tensor,

.. figure:: /figures/322_goma_physics.png
	:align: center
	:width: 90%

associated with the Navier-Stokes equations. This variable is called **SHEAR** in the
output EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the mean shear rate.
**no**   Do not calculate the mean shear rate.
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   Mean shear rate = yes

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




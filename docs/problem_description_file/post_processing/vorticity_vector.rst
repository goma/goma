********************
**Vorticity Vector**
********************

::

   Vorticity Vector = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This option allows the user to output the vorticity vector to the output EXODUS II file.
It applies to problems with the fluid momentum equations (see *EQ = momentum**
cards). The output nodal variables are named **VORTX, VORTY, VORTZ**.

The permissible values for this postprocessing option are:

============= ================================================================
**yes**       Calculate the vorticity vectors and store in the output
              EXODUS II file.
**no**        Do not calculate the vorticity vectors.
============= ================================================================

------------
**Examples**
------------

This example card requests that vorticity vectors be written to the EXODUS II file:
::

   Vorticity Vector = yes

-------------------------
**Technical Discussion**
-------------------------

The vorticity vector function, :math:`\underline\omega` , is defined in terms of the velocity :math:`\underline\upsilon` as:

.. figure:: /figures/336_goma_physics.png
	:align: center
	:width: 90%




.. 
	TODO - Line 40 is a photo that needs to be replaced with the correct equation. 
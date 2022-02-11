***********************
**Velocity Divergence**
***********************

::

   Velocity Divergence = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The divergence of velocity is associated with local mass conservation or how well the
solenoidal character of the velocity field in *ARBITRARY* mesh motion regions is being
maintained. (Fluid momentum equations are only applied for this *Mesh Motion* option.)
Here we calculate the :math:`L_2` norm of the divergence of velocity so that it is always zero or positive. This variable is called **DIVV** in the output EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the velocity divergence.
**no**   Do not calculate the velocity divergence.
======== ===============================================

------------
**Examples**
------------

A sample input specification for this card is:
::

   Velocity Divergence = no

-------------------------
**Technical Discussion**
-------------------------

The divergence of the fluid velocity field is defined as the scalar :math:`\Delta` • :math:`\underline{v}`.



--------------
**References**
--------------

No References.
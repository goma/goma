****************
Map Conf Stress
****************

::

   Map Conf Stress = {yes | no}

-----------------------
**Description / Usage**
-----------------------


The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Output stress values calculated from conformation tensor
**no**   Disable
======== ===============================================

------------
**Examples**
------------

This is an example of the input to request stress mappings be written to the EXODUS II file.
::

   Map Conf Stress = yes

-------------------------
**Technical Discussion**
-------------------------

Outputs the stress values from the conformation tensor mapping. 

.. math::

   \tau_i = {\mu^i_p}{\lambda_i} (C - I)

Where :math:`C` is the conformation tensor and `i` indicates the mode

Output exodus variables are MSab_i where ab are the components and i is the mode.

--------------
**References**
--------------

No References.
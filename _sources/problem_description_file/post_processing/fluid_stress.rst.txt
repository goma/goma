************
Fluid Stress
************

::

	Fluid Stress = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The Fluid Stress tensor without the pressure contribution. Outputs `FS11..FS33`

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the fluid stress
**no**   Do not calculate the fluid stress
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   Fluid Stress = yes

-------------------------
**Technical Discussion**
-------------------------

For Newtonian and Generalized Newtonian fluids the fluid stress tensor is given by:

.. math::
    \sigma = \mu \left( \nabla u + \nabla u^T \right)

where :math:`\mu` is the dynamic viscosity of the fluid.

For a viscoelastic fluid the fluid stress tensor is given by:

.. math::
    \sigma = \mu \left( \nabla u + \nabla u^T \right) + \tau 

where :math:`\mu` is the dynamic viscosity of the fluid and :math:`\tau` is the viscoelastic stress tensor.


--------------
**References**
--------------

No References.
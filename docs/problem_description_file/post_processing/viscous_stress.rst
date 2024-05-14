**************
Viscous Stress
**************

::

	Viscous Stress = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The Viscous Stress tensor. Output `VS11..VS33`

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the viscous stress
**no**   Do not calculate the viscous stress
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   Viscous Stress = yes

-------------------------
**Technical Discussion**
-------------------------

Calculates the viscous stress tensor:

.. math::
    \sigma = \mu \left( \nabla u + \nabla u^T \right)

where :math:`\mu` is the dynamic viscosity of the fluid.


--------------
**References**
--------------

No References.
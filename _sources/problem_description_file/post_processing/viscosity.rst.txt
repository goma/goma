*************
**Viscosity**
*************

::

   Viscosity = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This option allows you to plot the viscosity, which is written to the *Output EXODUS II*
file as the variable **MU**. This is a useful feature for non-Newtonian fluids such as
Phillip’s model for suspensions, Bingham plastic models, polymerizing solutions and
other materials for which the viscosity may change orders of magnitude, greatly
affecting the velocity and pressure fields. Contouring this variable **MU** over the domain
can be useful in explaining some physical phenomena.

he permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the viscosity and output as a nodal 
         variable in the Output EXODUS II file.
**no**   Do not calculate the viscosity.
======== ===============================================

------------
**Examples**
------------

The following sample card requests **MU** be written to the EXODUS II file:
::

   Viscosity = yes

-------------------------
**Technical Discussion**
-------------------------

See the material file *Viscosity* card for an explanation of the models for which the
viscosity is variable and dependent on the flow field and other variables.




**********************
**Mesh Stress Tensor**
**********************

::

   Mesh Stress Tensor = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The mesh stress tensor is associated with the equations of elasticity. The stress tensor
has six entries (in three dimensions, because it is symmetric) called **T11, T22, T33, T12, T13,** and **T23** in the output EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the mesh stress tensor and write to 
         output EXODUS II file.
**no**   Do not calculate the mesh stress tensor.
======== ===============================================

------------
**Examples**
------------

The following sample card turns on the writing of the stress tensor to the EXODUS II
file:
::

   Mesh Stress Tensor = yes

-------------------------
**Technical Discussion**
-------------------------

The defining constitutive equations for these stresses can be found in the description
for the *Solid Constitutive Equation* card. This option applies to all solid-material types
(see *Mesh Motion* card), viz. *TOTAL_ALE, LAGRANGIAN, ARBITRARY,
DYNAMIC_LAGRANGIAN*. In the *TOTAL_ALE* and *ARBITRARY* mesh motion types,
the mesh stress is exactly that and not the true stress of the material. For *TOTAL_ALE*
mesh motion types, use Real Solid Stress Tensor option to get the true solid material
stresses.



--------------
**References**
--------------

No References.
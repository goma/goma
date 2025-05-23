*************************
**FSI Deformation Model**
*************************

::

   FSI Deformation Model = {model_name}

-----------------------
**Description / Usage**
-----------------------

This card specifies the type of interaction the lubrication shell elements will have with
any surrounding continuum element friends. When not coupling the lubrication
equations to a continuum element, this card should be set to the default value,
FSI_SHELL_ONLY. All models are described below:

**FSI_MESH_BOTH** This model should be used when both the shell
and neighboring continuum elements use deformable meshes and the user wishes to
full couple these behaviors. This model is not currently implemented and should not be
used.

**FSI_MESH_CONTINUUM** In this model, the neighboring continuum
elements use mesh equations, but the lubrication shell does not. This model features a
two-way coupling, where the lubrication pressure can deform the neighboring solid
(through the appropriate boundary condition) and deformations to the mesh in turn
affect the height of the lubrication gap. This is equivalent to the old “toggle = 1”.

**FSI_MESH_SHELL** This model accounts for mesh equations present
in the lubrication shell, but not in the adjoining continuum elements. This model is not
currently implemented and should not be used.

**FSI_SHELL_ONLY** This model can be thought of as the default
behavior, where there is no coupling between the lubrication shell elements and any
neighboring continuum elements. This should also be used if only shells are present.

**FSI_MESH_UNDEF** This model is similar to
FSI_MESH_CONTINUUM, but the normal vectors in the shell are calculated using
the original undeformed configuration, rather than the current deformed state.
Implementation of this model is currently in progress and needs to be fully tested.

**FSI_MESH_ONEWAY** This model is similar to
FSI_MESH_CONTINUUM, but only utilizes a one way coupling. Deformations in the
neighboring continuum element do not affect the lubrication height, but do affect the
calculated normal vectors. This is equivalent to the old “toggle = 0”.

------------
**Examples**
------------

-------------------------
**Technical Discussion**
-------------------------



--------------
**References**
--------------

No References.
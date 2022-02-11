********************************
**Viscoplastic Def_Grad Tensor**
********************************

::

   Viscoplastic Def_Grad Tensor = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The components of this tensor are associated with the elasto-viscoplasticity model,
(described in detail in Schunk, et. al., 2001). If the mesh motion is of *LAGRANGIAN*
type, then this card activates the components of this tensor to be available in the
postprocessing EXODUS II file (see *Mesh Motion* card). The components are called
**FVP11, FVP12, FVP21, FVP22,** and **FVP33**. This tensor is the identity tensor in
regions that have not yielded, and so the diagonal components are unity; in regions that
have yielded, these components deviate from the identity. Contouring them can reveal
regions of plastic flow.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the viscoplastic Def_Grad tensor and 
         write components in the output EXODUS II file.
**no**   Do not calculate the viscoplastic Def_Grad tensor.
======== ===============================================

------------
**Examples**
------------

This sample input card does not activate Def_Grad output to the EXODUS II file:
::

   Viscoplastic Def_Grad Tensor = no

-------------------------
**Technical Discussion**
-------------------------

Please see complete discussion in Schunk, et. el. (2001).



--------------
**References**
--------------

GT-019.1: Elastoviscoplastic (EVP) Constitutive Model in GOMA: Theory, Testing,
and Tutorial, P. R. Schunk, A. Sun, S. Y. Tam (Imation Corp.) and K. S. Chen, January
11, 2001
*************************
**Moving Mesh Residuals**
*************************

::

   Moving Mesh Residuals = {yes | no}

-----------------------
**Description / Usage**
-----------------------

These nodal variables are constructed from the corresponding weighted residual
functions of the solid momentum equations (activated with the *EQ = mesh** cards).
The weighted residuals are formed using a Galerkin finite-element formulation. In the
output EXODUS II file they appear as nodal variables **RDX**, **RDY**, and **RDZ**,
corresponding to each of the independent components of the solid momentum balance
(both pseudo and real).

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Include the moving mesh residuals as nodal 
         variables in the ouput EXODUS II file.
**no**   Do not include moving mesh residuals.
======== ===============================================

------------
**Examples**
------------

Following is a sample card which does not activate writing of mesh residuals:
::

   Moving Mesh Residuals = no

-------------------------
**Technical Discussion**
-------------------------

This option can be used to help understand convergence behavior of a particular
problem, as it allows the user to visualize the pattern of residuals over the
computational domain during a Newton iteration process. The intermediate solutions
of a Newton iteration process can be activated with the *Write Intermediate Results* card.
Contouring these residuals can indicate where the convergence of a problem is being
delayed, and give the user/developer some clues as to the boundary condition or local
region of the mesh which is responsible.




***************************
**Navier Stokes Residuals**
***************************

::

   Navier Stokes Residuals = {yes | no}

-----------------------
**Description / Usage**
-----------------------

These post-processing nodal variables are constructed from the corresponding
weighted residual function of the fluid (e.g. Navier-Stokes) momentum equations,
using a Galerkin finite-element formulation. When activated with this card, variables
named **RMX**, **RMY**, and **RMZ** appear in the output EXODUS II file, corresponding to
each of the independent components of the fluid momentum balance.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the Navier-Stokes residuals and store
         as nodal variables in the output EXODUS II file.
**no**   Do not calculate Navier-Stokes residuals.
======== ===============================================

------------
**Examples**
------------

Following is a sample input card:
::

   Navier Stokes Residuals = no

-------------------------
**Technical Discussion**
-------------------------

This option can be used to help understand convergence behavior of a particular
problem, as it allows the user to visualize the pattern of residuals over the
computational domain during a Newton iteration process. The intermediate solutions of
a Newton iteration process can be activated with the *Write Intermediate Results* card.




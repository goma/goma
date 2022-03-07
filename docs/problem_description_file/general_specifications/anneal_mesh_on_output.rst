*************************
**Anneal Mesh on Output**
*************************

::

	Anneal Mesh on Output = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This optional card enables the user to specify that the mesh displacements should be set
to zero for the next continuation step. Valid options for this input card are

yes                      
    Set the mesh displacements to zero for the next continuation step.
    
no
    Do not set the mesh displacements to zero for the next continuation step.
    This is the default.

There are two important restrictions: 1) annealing the mesh will not work for the
TOTAL_ALE mesh motion types, and 2) only the last time step will be annealed for
transient problems.

------------
**Examples**
------------

Following is a sample card:
::

	Anneal Mesh on Output = yes

-------------------------
**Technical Discussion**
-------------------------

Annealing a mesh is accomplished by adding the displacements (from the
solution) to the base positions (from the FEM file) and writing the resulting
nodal positions to a new EXODUS II file, currently anneal.exoII. During the
annealing process, the displacement field is also set to zero. This file would
be used to restart a subsequent analysis where the anneal.exoII is copied to,
or becomes, the file used in a **read_exoII** option for an *Initial Guess*.



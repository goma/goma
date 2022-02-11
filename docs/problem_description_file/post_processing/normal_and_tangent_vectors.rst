******************************
**Normal and Tangent Vectors**
******************************

::

   Normal and Tangent Vectors = {yes | no}

-----------------------
**Description / Usage**
-----------------------

This option allows one to write the values of the normal and tangent vectors used in
rotating the mesh and momentum equations as nodal variables to the output EXODUS
II file. In two-dimensional problems, the normal and tangent vectors are saved as **N1, N2, N3** and **T1, T2, T3** in the output EXODUS II file; in two dimensions these vectors
are calculated at all the nodes. In three-dimensional problems, the normal and tangent
vectors are saved as **N1, N2, N3, TA1, TA2, TA3,** and **TB1, TB2, TB3**; in three
dimensions, these vectors only exist at nodes with rotation specifications, and the
vectors correspond to the rotation vectors chosen by the *ROT Specifications* for the
given node (see description for *ROT* cards). Thus in three-dimensional problems,
vectors are not necessarily saved for every node, nor do the vectors necessarily
correspond to the normal, first tangent, and second tangent, respectively.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the vectors and store as nodal 
         variables in the output EXODUS II file.
**no**   Do not calculate the vectors.
======== ===============================================

------------
**Examples**
------------

The following sample card produces no output to the EXODUS II file:
::

   Normal and Tangent vectors = no

-------------------------
**Technical Discussion**
-------------------------

This option is mostly used to debug three-dimensional meshes for full threedimensional
ALE mesh motion. The tangent fields in 3D should be smooth across the
surfaces, and *Goma* takes many steps to make them so. The surface normal crossed into any vector that is different will produce one tangent vector. Then the normal crossed (viz. cross product of two vectors) with the first tangent will produce a second tangent vector. Because the surface tangent basis fields are not unique, they must be uniform over a surface when the rotated
Galerkin weighted residuals are formed (see description for *ROT* cards). Imperfections
or defects in the mesh can lead to nonsmooth fields.



--------------
**References**
--------------

GT-018.1: ROT card tutorial, January 22, 2001, T. A. Baer
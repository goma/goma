*********************
**Pressure Contours**
*********************

::

	Pressure contours = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The hydrodynamic pressure is normally a field variable within *Goma*; however, it is
often interpolated in finite element space with discontinuous basis functions (in order
to satisfy the well-known LBB stability criterion, cf. Schunk, et al. 2002). This option
enables interpolating and smoothing the hydrodynamic pressure to nodal values that
most post-processors can deal with (e.g. BLOT, Mustafa). This variable is called
**PRESSURE** in the output EXODUS II file.

The permissible values for this postprocessing option are:

======== ===============================================
**yes**  Calculate the pressure contours.
**no**   Do not calculate the pressure contours.
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   Pressure contours = no

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

SAND2001-3512J: Iterative Solvers and Preconditioners for Fully-coupled Finite
Element Formulations of Incompressible Fluid Mechanics and Related Transport
Problems, P. R. Schunk, M. A. Heroux, R. R. Rao, T. A. Baer, S. R. Subia and A. C.
Sun. (March 2002)
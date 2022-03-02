******************
**APR_PLANE_TRAN**
******************

::

	BC = APR_PLANE_TRAN SS <bc_id> <float1>

-----------------------
**Description / Usage**
-----------------------

**(WIC/SCALAR ACOUS_PREAL)**

This boundary condition card applies the plane wave transmission conditions to the acoustic wave equations. This card concerns the real part while API_PLANE_TRAN concerns the imaginary component. This condition is used to set reflection/ transmission conditions for a surrounded material that is not being meshed. Definitions of the input parameters are as follows:

APR_PLANE_TRAN
    Name of the boundary condition (<bc_name>).
SS             
    Type of boundary condition (<bc_type>), where 
    **SS** denotes side set.

<bc_id>
    The boundary flag identifier, or a side set number which is an integer that
    identifies the boundary location (side set in EXODUS II) in the problem
    domain.

<float1>
    :math:`R_2`, the acoustic impedance (i.e. product of density and wave
    speed) in the surrounded material.

------------
**Examples**
------------

The following is a sample input card:
::

   BC = APR_PLANE_TRAN SS 10 0.1

-------------------------
**Technical Discussion**
-------------------------

This condition should be used to account for transmission/reflection conditions for the external boundaries when the acoustic wave equation is used. It reflects characteristics for an acoustic wave encountering a planar interface between two materials;

.. figure:: /figures/254_goma_physics.png
	:align: center
	:width: 90%

where k is the acoustic wavenumber and R is the acoustic impedance. The subscript 1 refers to the material inside the external boundary and is the material which is meshed. Subscript 2 refers to the material outside of the external boundary. If :math:`R_2` is set equal to
:math:`R_1`, then this condition mimics an infinite boundary condition, i.e. no reflection at the external boundary.


.. TODO - Line 43 has an image that needs to be replaced with the equation. 

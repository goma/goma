****************************
**Streamwise Normal Stress**
****************************

::

	Streamwise normal stress = {yes | no}

-----------------------
**Description / Usage**
-----------------------

The stream-wise normal stress, :math:`T_{tt}`, is defined as tt: :math:`\tau`, where t is the unit tangent vector
to the streamlines computed as v ⁄ :math:`\mid` v :math:`\mid` and :math:`\tau` is the deviatoric part of the dissipative stress
tensor,

.. figure:: /figures/320_goma_physics.png
	:align: center
	:width: 90%

associated with the Navier-Stokes equations. This variable is called SNS in the output
EXODUS II file.

The permissible values for this postprocessing option are

======== ===============================================
**yes**  Calculate the stream-wise normal stress.
**no**   Do not calculate the stream-wise normal stress.
======== ===============================================

------------
**Examples**
------------

Following is a sample card:
::

   Streamwise normal stress = yes

-------------------------
**Technical Discussion**
-------------------------

As of 2/9/02 this function is computed with the based viscosity, and not the strain-rate
dependent viscosity as might be the case for viscosity models other than *NEWTONIAN*
(see *Viscosity* card).




..
	TODO - Line 17 contains a photo that needs to be written as an equation.
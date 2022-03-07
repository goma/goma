*******************************
**Diffusion Coefficient Model**
*******************************

::

   Diffusion Coefficient Model = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card takes the specification of the diffusion coefficient model for the conservation
of particles inside film-flow capability, i.e. equation describing shell_partc.
Currently two models for {model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model specifies a constant diffusion coefficient. This option requires one      |
|                          |floating point values                                                                |
|                          |                                                                                     |
|                          | * <float1> is the diffusion coefficient                                             |
+--------------------------+-------------------------------------------------------------------------------------+
|**STOKES_EINSTEIN**       |This model specifies diffusion coefficient that depends on particles radius and the  |
|                          |film viscosity. The functional form is:                                              |
|                          |                                                                                     |
|                          | * <float1> is the Boltzmann constant kb where the magnitude depends on the units    |
|                          |   chosen by the user.                                                               |
|                          | * <float2> is temperature T in unit of Kelvin.                                      |
|                          | * <float3> is the particles radii R.                                                |
+--------------------------+-------------------------------------------------------------------------------------+

.. figure:: /figures/489_goma_physics.png
	:align: center
	:width: 90%

------------
**Examples**
------------

Following is a sample card:

::

   Diffusion Coefficient Model = STOKES_EINSTEIN 1.3807e-16 298 1.0e-6

This results in diffusion coefficient calculated with Stokes Einstein model with
Bolztmann constant of 1.3807e-16 in CGS units, 298 K temperature, and 1.0e-6 cm
radius particles.

-------------------------
**Technical Discussion**
-------------------------

Viscosity dependence of diffusion coefficient can be exploited to relate particles
concentration (or volume fraction in this case) to diffusion coefficient by employing
SUSPENSION viscosity model in the material file. See SUSPENSION viscosity model for
further detail




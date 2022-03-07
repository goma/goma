**************************
**Shear Rate Diffusivity**
**************************

::

   Shear Rate Diffusivity = {model_name} <species> <float>

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the coefficient for the shear-rate gradient term when
**HYDRO** is specified in the *Diffusivity* card. Definitions of the input parameters follow
for the {model_name} options **CONSTANT** and **LINEAR** based on the model:

.. figure:: /figures/447_goma_physics.png
	:align: center
	:width: 90%

+----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**          |Name of the model for constant shear rate diffusivity.                               |
|                      |                                                                                     |
|                      | * <species> - an integer designating the species equation.                          |
|                      | * <float> - Dc when there is no concentration dependency.                           |
+----------------------+-------------------------------------------------------------------------------------+
|**LINEAR**            |Name of the model in which shear rate diffusivity is a linear function of            |
|                      |concentration.                                                                       |
|                      |                                                                                     |
|                      | * <species> - an integer designating the species equation.                          |
|                      | * <float> - kc when the diffusivity is a linear function of concentration.          |
+----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Shear Rate Diffusivity = CONSTANT 0 0.

-------------------------
**Technical Discussion**
-------------------------

Please refer to the technical discussion given under **HYDRO** section of the *Diffusivity* card.



--------------
**References**
--------------

No References.
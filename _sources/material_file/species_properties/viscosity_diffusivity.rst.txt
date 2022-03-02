*************************
**Viscosity Diffusivity**
*************************

::

   Viscosity Diffusivity = {model_name} <species> <float>

-----------------------
**Description / Usage**
-----------------------

This card is used to specify Dr when the model in the *Diffusivity* card is **HYDRO**.
Definitions of the input parameters follow for the {model_name} options **CONSTANT**
and **LINEAR** based on the model:

.. figure:: /figures/449_goma_physics.png
	:align: center
	:width: 90%

+----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**          |Name of the model for a constant curvature diffusivity.                              |
|                      |                                                                                     |
|                      | * <species> - An integer designating the species equation.                          |
|                      | * <float> - Dμ when there is no concentration dependency.                           |
+----------------------+-------------------------------------------------------------------------------------+
|**LINEAR**            |Name of the model in which the diffusivity is a linear function of concentration.    |
|                      |                                                                                     |
|                      | * <species> - an integer designating the species equation.                          |
|                      | * <float> - kμ when the diffusivity is a linear function of concentration.          |
+----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Viscosity Diffusivity = CONSTANT 0 0.

-------------------------
**Technical Discussion**
-------------------------

Please refer to the technical discussion given under **HYDRO** section of the Diffusivity
card.



--------------
**References**
--------------

No References.
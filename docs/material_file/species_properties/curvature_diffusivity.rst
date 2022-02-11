*************************
**Curvature Diffusivity**
*************************

::

   Curvature Diffusivity = {model_name} <species> <float>

-----------------------
**Description / Usage**
-----------------------

This card is used to specify Dr when the model in the *Diffusivity* card is **HYDRO**.
Definitions of the input parameters follow for the {model_name} options **CONSTANT**
and **LINEAR** based on the model:

.. figure:: /figures/448_goma_physics.png
	:align: center
	:width: 90%

+----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**          |Name of the model for a constant curvature diffusivity.                              |
|                      |                                                                                     |
|                      | * <species> - An integer designating the species equation.                          |
|                      | * <float> - Dr when there is no concentration dependency.                           |
+----------------------+-------------------------------------------------------------------------------------+
|**LINEAR**            |Name of the model in which the diffusivity is a linear function of concentration.    |
|                      |                                                                                     |
|                      | * <species> - an integer designating the species equation.                          |
|                      | * <float> - kr when the diffusivity is a linear function of concentration.          |
+----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Curvature Diffusivity = CONSTANT 0 0.

-------------------------
**Technical Discussion**
-------------------------

It was proposed that adding a curvature contribution of the diffusive flux for
suspension particles would correct suspension migration behavior in parallel-plate and
cone-and-plate. However, this correction term is not frame-invariant; hence, it cannot
be used in generalized flow geometry. It is therefore not recommended.



--------------
**References**
--------------

No References.
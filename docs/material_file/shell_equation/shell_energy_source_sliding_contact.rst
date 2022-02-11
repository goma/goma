***************************************
**Shell Energy Source Sliding Contact**
***************************************

::

   Shell Energy Source Sliding_Contact = {model_name} <float_list>

-----------------------
**Description / Usage**
-----------------------

This card activates a heat source (or sink, as it were) in the shell_energy
equation. The functional form of this source/sink is a sliding contact model derived in
the frame-of-reference of the slider on a stationary surface, so that the surface is
moving in the simulation. In this case, the conditions for the flux vary from the leading
edge to the trailing edge of the slider as a thermal boundary layer builds up. Think of
this as a hot slider moving over a cold stationary wall, so that the flux at the leading
edge of the slider into the cold wall will be larger, due to a steeper thermal boundary
layer. Clearly the contact time will play a role. Currently two models {model_name}
are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**LOCAL_CONTACT**         |This model invokes the following functional form:                                    |
|                          |                                                                                     |
|                          |Commensurately there are seven floats required:                                      |
|                          |                                                                                     |
|                          | * <float1> - length l of slider.                                                    |
|                          | * <float2> - Sink temperature of substrate.                                         |
|                          | * <float3> - Thermal conductivity of substrate.                                     |
|                          | * <float4> - Density of substrate                                                   |
|                          | * <float5> - Heat capacity of substrate.                                            |
|                          | * <float6> - Delta L, or L1 – L2. This parameter sets the segment size (less than   |
|                          |   the total slider length) over which the heat flux is resolved.                    |
|                          | * <float7> - Leading edge coordinate of slider.                                     |
+--------------------------+-------------------------------------------------------------------------------------+

.. figure:: /figures/492_goma_physics.png
	:align: center
	:width: 90%

+--------------------------+-------------------------------------------------------------------------------------+
|**AVERAGE_CONTACT**       |This model invokes the following functional form:                                    |
|                          |                                                                                     |
|                          |Commensurately there are seven floats required:                                      |
|                          |                                                                                     |
|                          | * <float1> - length l of slider.                                                    |
|                          | * <float2> - Sink temperature of substrate.                                         |
|                          | * <float3> - Thermal conductivity of substrate.                                     |
|                          | * <float4> - Density of substrate                                                   |
|                          | * <float5> - Heat capacity of substrate.                                            |
+--------------------------+-------------------------------------------------------------------------------------+

.. figure:: /figures/493_goma_physics.png
	:align: center
	:width: 90%

------------
**Examples**
------------

Following is a sample card:

::

   Shell Energy Source Sliding_Contact = LOCAL_CONTACT {L= 2.5} {t_r=20} {t_cond_cu_cgs} {density_cu_cgs} {heat_capacity_cu_cgs} {delta_L = 0.1} {leading_edge_coordx = 2.5}

-------------------------
**Technical Discussion**
-------------------------

Technical Discussion
This boundary condition was derived using the analytical solution for heat conduction
into an infinite slab, as derived by Carslaw and Jaeger. The modification here is that
the temperature source accommodates a motion relative to the substrate, which is what
leads to the need to segment the slider into bins over which a local heat flux solution is
derived.

NOTE: If this card is used and there is no upper-wall or lower wall sliding motion, and
error is thrown.




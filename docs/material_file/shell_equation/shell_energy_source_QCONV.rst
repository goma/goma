*****************************
**Shell Energy Source QCONV**
*****************************

::

   Shell Energy Source QCONV = {model_name} <float_list>

-----------------------
**Description / Usage**
-----------------------

This card activates a heat source (or sink, as it were) in the shell_energy
equation. The functional form of this source/sink is a lumped heat-transfer coefficient
model, hence the QCONV in its name (see BC = QCONV card in main user manual).
Currently two models {model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model invokes a simple constant heat-transfer coefficient and reference         |
|                          |temperature, viz.                                                                    |
|                          |                                                                                     |
|                          |q = H(T − Tref ). Commensurately there are two floats required:                      |
|                          |                                                                                     |
|                          |<float1> - Heat transfer coefficient in units of Energy/time/L2/deg T E.g. W/m2-K in |
|                          |MKS units.                                                                           |
|                          |                                                                                     |
|                          |<float2> - Reference temperature.                                                    |
+--------------------------+-------------------------------------------------------------------------------------+
|**MELT_TURB**             |This model also invokes a lumped parameter model, but the heat-transfer coefficient  |
|                          |depends on the flow strength (Reynolds number), viz.                                 |
|                          |                                                                                     |
|                          |q = H(T − Tref ). Three floats are required:                                         |
|                          |                                                                                     |
|                          |<float1> - Thermal conductivity in units of Energy time/L/deg (e.g. W/m/k).          |
|                          |                                                                                     |
|                          |<float2> - Reference temperature.                                                    |
|                          |                                                                                     |
|                          |<float3> - Latent heat of melting (Energy/M, e.g. J/Kg). This quantity is required   |
|                          |due to the cross use of this in the shell_deltah equation (viz. EQ =  shell_deltah). |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Shell Energy Source QCONV = MELT_TURB {thermal_conductivity} {Tref} {latent_heat}

-------------------------
**Technical Discussion**
-------------------------

The MELT_TURB model warrants further discussion. The functional form of the heat
transfer coefficient H is

.. figure:: /figures/491_goma_physics.png
	:align: center
	:width: 90%

Here cf is the coefficient of friction, which for now is taken as 8./Re.




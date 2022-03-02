******************************
**Turbulent Lubrication Mode**
******************************

::

   Turbulent Lubrication Model = {model_name}

-----------------------
**Description / Usage**
-----------------------

This card activates a turbulent model for viscosity in the lub_p equation. Currently
one model {model_name} is permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**PRANDTL_MIXING**        |This model is used to determine the pre-multiplier on the molecular viscosity in the |
|                          |Reynolds lubrication equation. For confined, laminar flow, this multiplier is 12. For|
|                          |turbulent flow it is taken as K(Re), where Re is the local Reynolds number.          |
|                          |Specifically, invoking a analytical approximation for K from Hirs (1973), we set k0  |
|                          |according to the:                                                                    |
|                          |                                                                                     |
|                          |Reynolds number Re= :ρh U μ                                                          |
|                          |                                                                                     |
|                          |For 0 < Re < 2000 K0=12 (Laminar case),                                              |
|                          |                                                                                     |
|                          |Else 2000 < Re < 100000 K0 = 0.05Re 3/4                                              |
|                          |                                                                                     |
|                          |Here the wall velocity is used to compute The Reynolds number, as this turbulence    |
|                          |model is specific to turbulent Couette flow.                                         |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Turbulent Lubrication Model = PRANDTL_MIXING

-------------------------
**Technical Discussion**
-------------------------

Several other models can be implemented in this instance. We chose this simple
model which derives from Prandtl mixing length theory.



--------------
**References**
--------------

G.G. Hirs, “Bulk flow theory for turbulence in lubricant films”, Trans. ASME, ser. F,
95, pp 137-146, 1973.


********************************
**Porous Shell Gas Diffusivity**
********************************

::

   Porous Shell Gas Diffusivity = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card is used to set the gas diffusivity for the trapped gas in the
porous_sat_closed equation. Basically, the gas trapped in closed pores during
the imbibition process is allowed to diffuse into the liquid, and this property is a part of
that model gas inventory equation R_SHELL_SAT_GASN. Only one model is
available for this property:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model applies a constant cross-region permeability. It requires a single        |
|                          |floating point input:                                                                |
|                          |                                                                                     |
|                          | * <float1> is the cross region permeability                                         |
+--------------------------+-------------------------------------------------------------------------------------+
|**EXTERNAL_FIELD**        |This model applies a constant gas diffusivity. It requires a single floating point   |
|                          |input:                                                                               |
|                          |                                                                                     |
|                          | * <float1> is the gas diffusivity (L2/t)                                            |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Porous Shell Gas Diffusivity = CONSTANT 1.e-5




--------------
**References**
--------------

S. A. Roberts and P. R. Schunk 2012. in preparation.
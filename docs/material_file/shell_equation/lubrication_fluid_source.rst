****************************
**Lubrication Fluid Source**
****************************

::

   Lubrication Fluid Source = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card sets a fluid mass source term in the lub_p equation. Can be used to specify
inflow mass fluxes over the entire portion of the lubrication gap in which the lub_p
equation is active (over the shell material). This flux might be the result of an injection
of fluid, or even melting. Currently two models {model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model is used to set a constant fluid source in units of velocity. Only one     |
|                          |floating point value is required.                                                    |
|                          |                                                                                     |
|                          | * <float1> is the velocity of the fluid source.                                     |
+--------------------------+-------------------------------------------------------------------------------------+
|**MELT**                  |This model is used to set fluid source in units of Velocity which results from an    |
|                          |analytical model of lubricated melt bearing flow due to Stiffler (1959). Three       |
|                          |floating point values are required.                                                  |
|                          |                                                                                     |
|                          | * <float1> is load on the slider in units of pressure                               |
|                          | * <float2> is the Stiffler delta factor. Unitless but depends on the aspect ratio.  |
|                          | * <float3> is the length of the slider in the direction of the motion.              |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Lubrication Fluid Source = CONSTANT 180.






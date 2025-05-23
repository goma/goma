************************************
**Porous Shell Henrys Law Constant**
************************************

::

   Porous Shell Henrys Law Constant = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card is used to set the partitionconstant in Henry’s law for the trapped gas in the
porous_sat_closed equation. Basically, the gas trapped in closed pores during
the imbibition process is allowed to diffuse into the liquid, and this property is a part of
that model for the dissolution constant of gas in pressurized liquid. It is only needed if
the equation R_SHELL_SAT_GASN is used. Please consult the references for a
detailed explanation. Only one model is available for this property:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model sets the Henry’s law constant. It requires a single floating point        |
|                          |input:                                                                               |
|                          |                                                                                     |
|                          | * <float1> is the Henry’s law constant                                              |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Porous Shell Gas Temperature Constant= CONSTANT 1.e10




--------------
**References**
--------------

S. A. Roberts and P. R. Schunk 2012. in preparation.
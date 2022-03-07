*****************
**Heat Capacity**
*****************

::

   Heat Capacity = {model_name} {float_list} [E/MT]

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for the heat capacity. Definitions of the
input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|{model_name}     |Name of the model for the heat capacity. This parameter can have one of the following values: **ONSTANT**,  |
|                 |**USER**, or **ENTHALPY**.                                                                                  |
+-----------------+------------------------------------------------------------------------------------------------------------+
|{float_list}     |One or more floating point numbers (<float1> through <floatn>) whose values are determined by the selection |
|                 |for {model_name}. These are identified in the discussion of each model below.                               |
+-----------------+------------------------------------------------------------------------------------------------------------+

Thus,

+---------------------+--------------------------------------------------------------------------------------------------------+
|**CONSTANT** <float> |This option specifies a constant heat capacity. The {float_list} has a single value:                    |
|                     |                                                                                                        |
|                     | * <float1> - Heat capacity                                                                             |
+---------------------+--------------------------------------------------------------------------------------------------------+
|**USER** <float1>... |the heat capacity will be a user-defined model. This user-defined model must be incorporated by         |
|                     |modifying the routine “usr_heat_capacity” in the file user_mp.c. The model parameters are entered in the|
|                     |{float_list} as <float1> through <floatn> and passed to the routine as an array.                        |
+---------------------+--------------------------------------------------------------------------------------------------------+
|**ENTHALPY** <float1>|a model of heat capacity that uses the latent heat of fusion parameter. The model goes as follows:      |
|                     |Here the {float_list} requires two values, where:                                                       |
|                     |                                                                                                        |
|                     | * <float1> - Base heat capacity in the solid state, cp                                                 |
|                     | * <float2> - Latent heat of fusion Hf                                                                  |
|                     |                                                                                                        |
|                     |The liquidus temperature Tl and the solidus temperature Ts are taken from the material file. This model |
|                     |is currently available for single species only, and is used for rapid melting problems in alloys.       |
+---------------------+--------------------------------------------------------------------------------------------------------+
|**TABLE** <integer1> |Please see discussion at the beginning of the material properties chapter 5 for input description and   |
|<character_string1>  |options. Most often character_string1 will be **TEMPERATURE** or maybe **MASS_FRACTION**.               |
|{LINEAR | BILINEAR}  |                                                                                                        |
|[integer2]           |                                                                                                        |
|[FILE = filenm]      |                                                                                                        |
+---------------------+--------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Heat Capacity = CONSTANT 1.

-------------------------
**Technical Discussion**
-------------------------

When the **ENTHALPY** option is used, the liquidus (Tl) and solidus (Ts) temperatures
must be added through the *Liquidus Temperature and Solidus Temperature* cards.



--------------
**References**
--------------

No References.
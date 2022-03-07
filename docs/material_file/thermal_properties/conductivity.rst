****************
**Conductivity**
****************

::

   Conductivity = {model_name} {float_list} [E/LtT]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for thermal conductivity. Definitions of the input
parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|{model_name}     |Name of the model for thermal conductivity; this parameter can have the value **CONSTANT** or **USER**.     |
+-----------------+------------------------------------------------------------------------------------------------------------+
|{float_list}     |One or more floating point numbers (<float1> through <floatn>) whose values are determined by the selection |
|                 |for {model_name}. These are identified in the discussion of each model below.                               |
+-----------------+------------------------------------------------------------------------------------------------------------+

Thus,

+---------------------+--------------------------------------------------------------------------------------------------------+
|**CONSTANT** <float> |a constant thermal conductivity model, {float_list} is a single value:                                  |
|                     |                                                                                                        |
|                     | * <float1> - Standard value of *k*                                                                     |
+---------------------+--------------------------------------------------------------------------------------------------------+
|**USER** <float1>... |a user-defined model. With the USER option the appropriate modifications to the routine                 |
|<floatn>             |“usr_thermal_conductivity” in the user_mp.c file must be undertaken. The {float_list} can be of         |
|                     |arbitrary length and is used to parameterize the model. These parameters are made available in the      |
|                     |subroutine via <float1> through <floatn>.                                                               |
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

   Conductivity = USER 1. 1. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
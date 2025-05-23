***************************
**Shell Bending Stiffness**
***************************

::

   Shell bending stiffness = {model_name} {float_list} [M/Lt2]

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for the Shell bending stiffness property
D which is defined as D=Et3/12(1-ν2), where E is the elastic modulus, ν Poisson’s
ratio, and t the shell thickness. The units are M-L2/t2 (or F-L). The elastic modulus is
set through the *Lame MU* and *Lame Lambda* cards. This property is needed for the
inextensible cylindrical shell equations (see *EQ = Shell Tension*).

Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|{model_name}     |Name of the Shell bending stiffness coefficient model. This parameter can have one of the following values: |
|                 |**CONSTANT**.                                                                                               |
|                 |                                                                                                            |
|                 | * {float1} - The value of the shell bending stiffness.                                                     |
+-----------------+------------------------------------------------------------------------------------------------------------+

The details of each model option are given below:

+---------------------+--------------------------------------------------------------------------------------------------------+
|**CONSTANT** <float1>|For the **CONSTANT** model, {float_list} is a single value:                                             |
|                     |                                                                                                        |
|                     | * <float1> - Standard value of the coefficient D.                                                      |
+---------------------+--------------------------------------------------------------------------------------------------------+



--------------
**References**
--------------

GT-027.1: GOMA’s Shell Structure Capability: User Tutorial (GT-027.0). P. R.
Schunk and E. D. Wilkes.

GT-033.0: Structural shell application example: tensioned-web slot coater (GT-033.0).
P. R. Schunk and E. D. Wilkes.
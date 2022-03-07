************************************
**Porous Shell Closed Gas Pressure**
************************************

::

   Porous Shell Closed Gas Pressure = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card specifies the gas pressure used in porous_shell_closed equations.
Currently one model for {model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model applies a constant gas pressure for the entire model. It requires a single|
|                          |floating point value.                                                                |
|                          |                                                                                     |
|                          | * <float1> is the gas pressure                                                      |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Porous Shell Closed Gas Pressure = CONSTANT 0.5




--------------
**References**
--------------

No References.
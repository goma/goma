***********************
**Porous Shell Height**
***********************

::

   Porous Shell Height = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card specifies height of the pores used in porous_shell_closed and
porous_shell_open equations. Currently two models for {model_name} are
permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model applies a constant pore height for the entire model. It requires a single |
|                          |floating point value.                                                                |
|                          |                                                                                     |
|                          | * <float1> is the pore height. L.                                                   |
+--------------------------+-------------------------------------------------------------------------------------+
|**EXTERNAL_FIELD**        |This model reads in an array of values for the height from an initial exodus file.   |
|                          |This allows for spatial variations in the parameter value.                           |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Porous Shell Height = CONSTANT 0.00001




--------------
**References**
--------------

No References.
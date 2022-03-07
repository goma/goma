********************************
**Porous Shell Closed Porosity**
********************************

::

   Porous Shell Closed Porosity = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card specifies the porosity used in porous_shell_closed equations.
Currently two models for {model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model applies a constant porosity for the entire model. It requires a single    |
|                          |floating point value.                                                                |
|                          |                                                                                     |
|                          | * <float1> is the porosity                                                          |
+--------------------------+-------------------------------------------------------------------------------------+
|**EXTERNAL_FIELD**        |FIELDThis model reads in an array of values for the porosity from an initial exodus  |
|                          |file. This allows for spatial variations in the parameter value.                     |
|                          |                                                                                     |
|                          | * <float1> scale factor for scaling field value                                     |
|                          |                                                                                     |
|                          |The ExodusII field variable name should be “SH_SAT_CL_POROSITY”, viz.                |
+--------------------------+-------------------------------------------------------------------------------------+

::

   External Field = SH_SAT_CLOSED_POROSITY Q1 name.exoII (see this card)

------------
**Examples**
------------

Following is a sample card:

::

   Porous Shell Closed Porosity= CONSTANT 0.5




--------------
**References**
--------------

No References.
***************************
Extinction Index
***************************

::

   Extinction Index = {model_name} {float} []

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for electrical permittivity when using REFRACTIVE_INDEX model.
Definitions of the input parameters are as follows:

CONSTANT
   Name of the model for constant extinction index
      * <float> -the value 
TABLE
    See TABLE instructions, WAVELENGTH is probably what you want to interpolate with

------------
**Examples**
------------

Following are sample cards:

::

   Extinction Index = CONSTANT 1.
   Extinction Index = TABLE 2 WAVELENGTH 0 LINEAR FILE=k.txt

-------------------------
**Technical Discussion**
-------------------------

This card is utilized to set the extinction index, used for acoustic and time-harmonic Maxwell.


--------------
**References**
--------------

No References.

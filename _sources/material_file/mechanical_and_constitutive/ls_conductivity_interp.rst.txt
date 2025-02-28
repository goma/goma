*****************************************************
Level Set Interface Conductivity Interpolation Method
*****************************************************

::

   Level Set Interface Conductivity Interpolation Method = {char_string}

-----------------------
**Description / Usage**
-----------------------

This specifies the interpolation method of the level set properties for Conductivity

Possible choiches for char_string:

LINEAR 
   Default linear interpolation :math:`p = (1-H(\phi)) p_- + H(\phi) p_+`

LOG
   Log interpolation :math:`p = p_-^{(1-H(\phi))} p_+^{H(\phi)}`


------------
**Examples**
------------

The following is a usage example for this card:

::

   Level Set Interface Conductivity Interpolation Method = LOG


-------------------------
**Technical Discussion**
-------------------------

The log interpolation has been found to behave better in some cases where
material property differences are large.




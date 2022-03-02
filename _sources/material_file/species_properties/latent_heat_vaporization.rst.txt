****************************
**Latent Heat Vaporization**
****************************

::

   Latent Heat Vaporization = CONSTANT <species> <float> [E/M]

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for the latent heat of vaporization for
each species. Definitions of the input parameters are as follows:

+-----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**           |Name of the constant latent heat of vaporization model.                              |
|                       |                                                                                     |
|                       | * <species> - an integer designating the species equation.                          |
|                       | * <float> - the value of the latent heat of vaporization.                           |
+-----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Latent Heat Vaporization = CONSTANT 0 0.0

-------------------------
**Technical Discussion**
-------------------------

See the discussion for the *Latent Heat Fusion* model.



--------------
**References**
--------------

No References.
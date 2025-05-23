**********************
**Latent Heat Fusion**
**********************

::

   Latent Heat Fusion = CONSTANT <species> <float> [E/M]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the latent heat of fusion for each species.
Thus an input deck may include several of these cards. Definitions of the input
parameters are as follows:

+-----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**           |Name of the latent heat of fusion model, the only one available.                     |
|                       |                                                                                     |
|                       | * <species> - an integer designating the species equation.                          |
|                       | * <float> - the value of the latent heat of fusion.                                 |
+-----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Latent Heat Fusion = CONSTANT 0 0.0

-------------------------
**Technical Discussion**
-------------------------

This card is used on a species-basis and is unrelated to the latent heat of fusion
specification for the **ENTHALPY** model of heat capacity. It is used to calculate the
standard state heat of formation for the species. A related important card is the *Latent
Heat Vaporization*.



--------------
**References**
--------------

No References.
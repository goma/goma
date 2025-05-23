****************************
**Species Volume Expansion**
****************************

::

   Species Volume Expansion = CONSTANT <species> <float> [1/T]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the coefficient of volume expansion
associated with the concentration of a particular species. This property is optional for
the **BOUSS** and **BOUSSINESQ** option on the *Navier-Stokes Source* card and, if
nonzero, will result in a buoyancy term to be added to the Navier-Stokes equation that
is apportioned to the species volume expansion coefficient, defined as the logarithmic
sensitivity of density to concentration, or (dlnr) ⁄ (dC).

+-----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**           |Name of the constant volume expansion coefficient model.                             |
|                       |                                                                                     |
|                       | * <species> - an integer designating the species equation.                          |
|                       | * <float> - the value of the constant expansion coefficient.                        |
+-----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Species Volume Expansion = CONSTANT 0 0.

-------------------------
**Technical Discussion**
-------------------------

WARNING: Please be aware that if the thermal volume expansion coefficient is also
nonzero, the buoyancy force will be augmented.



--------------
**References**
--------------

No References.
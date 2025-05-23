***************************
**Reference Concentration**
***************************

::

   Reference Concentration = CONSTANT <species> <float> []

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for the reference concentration, which is
required by the BOUSS option on the *Navier-Stokes Source* card. Definitions of the
input parameters are as follows:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |Model for a constant reference concentration.                                        |
|                          |                                                                                     |
|                          | * <species> - the species equation to which this specification applies.             |
|                          | * <float> - the value of the reference concentration.                               |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Reference Concentration = CONSTANT 0 0.

-------------------------
**Technical Discussion**
-------------------------

The Boussinesq model subtracts out the pressure head in its final equations. Thus, to
zeroth order, hydrodynamic pressure field doesn’t include a static variation in the
gravity direction due to the pressure head. But, the source term in the momentum
equations then becomes –g(ρ – ρo) instead of simply –gρ. The reference
concentration values entered via this card are used to evaluate ρo for use in calculating
the natural convective force due to concentration differences.

The card is also used in various places where a value for a species concentration is
needed. However, the species unknown variable is not included in the solution vector.



--------------
**References**
--------------

No References.
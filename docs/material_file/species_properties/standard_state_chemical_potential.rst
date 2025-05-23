*************************************
**Standard State Chemical Potential**
*************************************

::

   Standard State Chemical Potential = CONSTANT <integer> <float>

-----------------------
**Description / Usage**
-----------------------

This card sets the standard state chemical potential of a species, <integer>, in the
current material to a specified value, <float>. Currently, only the generic CONSTANT
model is implemented. However, extensions to polynomial expressions in the
temperature are easily implemented and forthcoming.

+-----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**           |Model name for the standard chemical state chemical potential model.                 |
|                       |                                                                                     |
|                       | * <species> - an integer designating the species equation.                          |
|                       | * <float> - the value of the chemical potential                                     |
+-----------------------+-------------------------------------------------------------------------------------+

The standard state chemical potential, μk, o(T) , which is defined to be only a function
of the temperature, is used in the evaluation of the definition of the pure species
chemical potential of species k, μk (T, P) , which in turn is used in the evaluation of
the mixture chemical potential of species k, μk(T, P, Xi) .

------------
**Examples**
------------

The following is a sample input card:

::

   Standard State Chemical Potential = CONSTANT 0 1.0

-------------------------
**Technical Discussion**
-------------------------

The values in this card are currently only applicable to the *IS_EQUIL_PSEUDORXN*
boundary condition.



--------------
**References**
--------------

No References.
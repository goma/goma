***********************************
**Pure Species Chemical Potential**
***********************************

::

   Pure Species Chemical Potential = {model_name} <integer>

-----------------------
**Description / Usage**
-----------------------

This card takes the specification of the standard state chemical potential, which is
defined as a function of temperature only, and completes the definition of the pure
species chemical potential by possibly adding in a pressure dependence. Two model
values are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**PRESSURE_INDEPENDENCE** |No pressure dependence to the pure species state when this value of {model_name} is  |
|                          |specified. The standard state chemical potential is equal to the pure species        |
|                          |chemical potential. The <integer> argument specifies the species subindex, k.        |
+--------------------------+-------------------------------------------------------------------------------------+
|**PRESSURE_IDEALGAS**     |The following expression holds for the pressure dependence:                          |
|                          |                                                                                     |
|                          | μk*(T, P) μ= k, o(T) + RTln(P ⁄ 1 atm )                                             |
|                          |                                                                                     |
|                          | The <integer> argument specifies the species subindex, k.                           |
|                          |                                                                                     |
|                          |The standard state chemical potential, μk, o(T), which is defined to be only a       |
|                          |function of the temperature, is used in the evaluation of the definition of the pure |
|                          |species chemical potential of species k, μk*(T, P) , which in turn is used in the    |
|                          |evaluation of the mixture chemical potential of species k, μk (T P Xi).              |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Pure Species Chemical Potential = PRESSURE INDEPENDENT 0

-------------------------
**Technical Discussion**
-------------------------

The values in this card are only applicable to the *IS_EQUIL_PSEUDORXN* boundary
condition currently.



--------------
**References**
--------------

No References.
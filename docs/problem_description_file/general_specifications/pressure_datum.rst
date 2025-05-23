******************
**Pressure Datum**
******************

::

	Pressure Datum = <float> { atm | torr | cgs }

-----------------------
**Description / Usage**
-----------------------

This card is used to set a thermodynamic pressure datum on fluid or solid mechanics
problems that calculate equations of state requiring a true value for the total pressure.
The total pressure is then defined as the sum of a constant base thermodynamic
pressure, specified by this card, and a variable hydrodynamic pressure calculated via
the pressure unknown. Definitions of the input parameters are as follows:

<float>                       
    Value of the thermodynamic pressure datum.

{ atm | torr | cgs }          
    Units of the float specified above.

------------
**Examples**
------------

Following is a sample card:
::

	Pressure Datum = 1.0 atm

-------------------------
**Technical Discussion**
-------------------------

The value of this variable is stored in the unified problem description structure in cgs
units. It is then used in consistency checks and as input into some equation of state
routines, such as the ideal gas equation of state routine.


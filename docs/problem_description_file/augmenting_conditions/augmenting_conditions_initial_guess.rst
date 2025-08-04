***********************************
Augmenting Conditions Initial Guess
***********************************

::

    Augmenting Conditions Initial Guess = {read | other}

-----------------------
Description / Usage
-----------------------

This card is used to direct the specification of the initial guess for values of the 
augmenting condition parameters.

**read**
    Initial guesses are read from the ASCII file identified on the 
    GUESS file card.

**other**
    Any other string that is not read, or if this card is missing 
    altogether, the initial guesses for the augmenting condition 
    values will be obtained from the input deck itself or the 
    material file.

Note that the computed values of the augmenting conditions are automatically written 
at the end of the SOLN file, so it is usually necessary when employing this option to 
copy the SOLN file to the GUESS file before restarting.

------------
Examples
------------

In the following example, the initial guesses for the augmenting condition values will 
be obtained from the input deck itself or the material file.

::

    Augmenting Conditions Initial Guess = none

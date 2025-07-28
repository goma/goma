***********************************
Number of Augmenting Conditions
***********************************

::

    Number of augmenting conditions = <integer>

-----------------------
Description / Usage
-----------------------

This card allows the user to specify the number of augmenting condition cards to be
read. This card must be present in order to run a problem employing augmenting
conditions.

<integer> 
    N, an integer parameter that is either the number of AC
    cards to be read or -1.

Goma will try to read N cards that start with: AC = between this card and a line having
only the string END OF AC. If Goma finds fewer than N cards before encountering this
string, it will stop with an error. More than N cards is fine; the excess are ignored. If N
is -1, Goma will read all the augmenting condition cards between this card and the
END OF AC card.

---------
Examples
---------

In the following example, the initial guesses for the augmenting condition values will 
be obtained from the input deck itself or the material file.

::

    Number of augmenting conditions = -1
    AC = FC 1 0 0 VOLUME_FLUX 4 {-PI}
    AC = FC  1 1 0 FORCE_TANGENT1 3 {3*5.0}
    END OF AC

In this example, Goma will read and execute two different augmenting conditions.

--------------------
Technical Discussion
--------------------

When doing arc length continuation with augmenting conditions, one additional AC
record will be automatically created so the number of reported augmenting conditions
will be incremented to N+1. In this case, N is still the number of AC cards to be read,
which would not include this additional AC (for which there is no input card).
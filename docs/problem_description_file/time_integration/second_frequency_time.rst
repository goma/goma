*************************
Second Frequency Time
*************************

::

	Second frequency time = <float1> <float2>

-----------------------
Description / Usage
-----------------------

This card allows the time between successive writings of the solution to change after a
specified time and is only used if the <integer> in the *Printing Frequency* card is set to
0. Definitions of input parameters are as follows:

<float1>
    Any number indicating the time at which the printing frequency should shift
    from that specified in the *Printing Frequency* card to <float2>.

<float2>
    Printing frequency in time units (same units as specified in the *delta_t*
    card) for printing the solution at times greater than <float1>.

------------
Examples
------------

The following is a sample card that will change the printing frequency to print every 3
time units after 15 time units:
::

	Second frequency time = 15. 3.


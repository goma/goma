**********************************
Turbulence Calculate Wall Distance
**********************************

::

    Turbulence Calculate Wall Distance = {yes, no}

-----------------------
Description / Usage
-----------------------

Use built in goma function to calculate wall distance for Turbulence models
otherwise an external field is required.

yes|on
    Use built in goma function to calculate wall distance for Turbulence models

no|off
    Default. Use external field to calculate wall distance for Turbulence models

------------
Examples
------------

This is a sample card:
::

    Turbulence Calculate Wall Distance = yes
    Turbulence Wall Node Sets = 3
    Turbulence Wall Side Sets = 5

-------------------------
Technical Discussion
-------------------------

--------------
References
--------------

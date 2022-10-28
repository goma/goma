***********************
Polymer Shock Capturing
***********************

::

   Polymer Shock Capturing = {NONE | DCDD <float1>}

-----------------------
**Description / Usage**
-----------------------

This optional card is used to specify the shock function for the polymer stress
equation. Valid options are

NONE
   No polymer shock capturing is used, default
DCDD
   DCDD shock capturing term is used (currently only in LOG_CONF/LOG_CONF_GRADV)
   * <float1> the scaling value of the shock capturing term
YZBETA
   YZBETA shock capturing term is used (currently only in SQRT_CONF)
   * <float1> the scaling value of the shock capturing term

Currently only available for log-conformation / sqrt-conformation formulations.

------------
Examples
------------

The following is a sample card that set the polymer shock function to **SUPG** and
demonstrates the required cards.

::

   Polymer Shock Capturing = DCDD 1.0

The following is a sample card that set the polymer shock function to **NONE**.
::

   Polymer Shock Capturing = NONE

-------------------------
**Technical Discussion**
-------------------------

No Discussion.





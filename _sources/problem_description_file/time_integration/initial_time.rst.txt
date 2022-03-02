****************
Initial Time
****************

::

	Initial Time = <float>

-----------------------
Description / Usage
-----------------------

This card sets the time at which the calculation starts. The input parameter is defined as

<float>
    Any number indicating the initial solution time (in the same units as
    specified in the *delta_t* card). An additional feature can be triggered if
    this float is specified to be negative, which triggers GOMA to look for the
    nearest restart time in the restart ExodusII database to use as the start
    time. Note that this option can only be used with Initial Guess options of
    read_exoII_file or read_exoII.

Normally, the value of <float> will be set to zero unless the problem is a continuation
of a previous transient problem.

------------
Examples
------------

The following is a sample card that shows a restart at 45 time units:
::

	Initial Time = 45.0

The following is a sample card that triggers Goma to look for a restart time of 10 time
units, or the closest time value to 10 time units, to start from:
::

	Initial Time = -10.0


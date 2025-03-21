***************************
Level Set Contact Tolerance
***************************

::

	Level Set Contact Tolerance = <float>

-----------------------
Description / Usage
-----------------------

This card specifies a tolerance for when VELO_NORMAL like BCs will let the gas leave the domain.

The gas is determined by the less dense fluid.

This is enabled when the level set in the gas phase is closer than `tolerance * length_scale * 2` 

<float>
    Specified tolerance, defaults to 0.


------------
Examples
------------

This is a sample card:

::

    # Undocumented behavior in Goma 7.9.0 and earlier
	Level Set Contact Tolerance = 0.5

-------------------------
Technical Discussion
-------------------------

--------------
References
--------------

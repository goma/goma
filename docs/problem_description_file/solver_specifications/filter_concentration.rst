************************
Filter Concentration
************************

::

	Filter Concentration = <integer> <float1> <float2>

-----------------------
Description / Usage
-----------------------

This optional card allows the user to enforce strict bounds on the concentration of a
specific species. The input parameters are defined as:

<integer>
    **i,** this integer indicates which species ( i ≥ 0 ) receives this special
    restriction.
<float1>
    **min,** a real number indicating the minimum concentration.
<float2>
    **max,** a real number indicating the maximum concentration.

There are no default values; concentrations take on whatever values are naturally
dictated by the Newton iterations.

------------
Examples
------------

The following is a sample card:
::

	Filter Concentration = 0 0.0 1.0

-------------------------
Technical Discussion
-------------------------

Although a correct solution should not have concentrations less than 0 or greater than
1.0, such values may arise in the solution vector due to various sources. Intermediate
solutions during the Newton iteration may cause non-physical values to arise.
Numerical error due to inexact linear solves, rounding, etc., may cause the values to be
inexact. This card allows the user to force the concentration of species **i** to be corrected
to fall within a strict concentration range **[min,max]** after the Newton iterations have
terminated.


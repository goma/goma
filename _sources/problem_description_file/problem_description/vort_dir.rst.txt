************
**vort_dir**
************

::

	EQ = vort_dir{1|2|3} {Galerkin_wt} {VD1|VD2|VD3} {Interpol_fnc}

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for one component of
the vorticity equation. Definitions of the input parameters are defined below; there is
no <float> input for this equation. The Galerkin weight and the interpolation function
must be the same for the code to work properly.

+-------------------------------------+------------------------------------------------+
|**vort_dir1 | vort_dir2 | vort_dir3**|Name of the equation to be solved, where the 1, | 
|                                     |2 and 3 components correspond to one of the     |
|                                     |principal coordinate directions, e.g. X, Y and  |
|                                     |Z for Cartesian geometry.                       |
+-------------------------------------+------------------------------------------------+
|{Galerkin_wt}                        |Two--character value that defines the type of   |
|                                     |weighting function for this equation, where:    |
|                                     |                                                |
|                                     | * **Q1**-Linear                                |
|                                     | * **Q2**-Quadratic                             |
+-------------------------------------+------------------------------------------------+
|**VD1 | VD2 | VD3**                  |Name of the variable associated with the 1, 2   |
|                                     |or 3 principal coordinate direction for this    |
|                                     |component equation.                             |
+-------------------------------------+------------------------------------------------+
|{Interpol_fnc}                       |Two-character value that defines the            |
|                                     |interpolation function used to represent the    |
|                                     |variable **VD1, VD2** or **VD3**, where:        |
|                                     |                                                |
|                                     | * **Q1**-Linear                                |
|                                     | * **Q2**-Quadratic                             |
+-------------------------------------+------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that uses a linear continuous interpolation and weight
function:
::

   EQ = vort_dir1 Q1 VD1 Q1

-------------------------
**Technical Discussion**
-------------------------

This equation type is used for a research capability involving the flows of suspensions
in curvilinear coordinates and is not currently being used for production computations.




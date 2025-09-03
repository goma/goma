*********************************
AC (Periodic Boundary Condition)
*********************************

::

    AC = PBC <var_name> <sideset_id1> <sideset_id2>

-----------------------
Description / Usage
-----------------------

This augmenting condition card allows the user to specify that a field variable has the 
same value at both ends of a given domain. The ends are specified as side sets, which 
must have the same number of nodes with the same spacing between them.

Definitions of the input parameters are as follows:

PBC
    Mandatory string identifying the augmenting condition 
    card as being of type Periodic Boundary Condition.

<var_name>
    A string indicating the variable whose value will be 
    matched at the specified periodic boundaries of the 
    problem domain. This variable must be active 
    throughout the region bordered by these two boundaries.

<sideset_id1>
    An integer parameter that identifies the sideset at one of 
    the two periodic boundaries for the specified variable. 
    This sideset must be coextensive with this domain 
    boundary.

<sideset_id2>
    An integer parameter that identifies the sideset at the 
    other periodic boundary for the specified variable. This 
    sideset must be coextensive with this domain boundary.

------------
Examples
------------

The following is an AC (Periodic Boundary Condition) card which specifies periodic 
boundaries for voltage at sidesets 20 and 30:

::

    AC = PBC VOLTAGE 20 30

-------------------------
Technical Discussion
-------------------------

When this augmenting condition is invoked, the degrees of freedom (DOF's) of the 
specified variable at corresponding points on the two periodic boundaries are matched. 
During equation assembly, a Lagrange multiplier constraint is applied to each pair of 
DOF's to enforce the constraint, which in residual form is:

.. math::

    \lambda(x_2 - x_1) = 0

where :math:`x_1` and :math:`x_2` are the unknown values at corresponding points along the two 
specified periodic boundaries, and λ is the added unknown associated with this 
constraint. This residual and all relevant sensitivities to problem unknowns are 
assembled into the appropriate augmenting condition arrays and solved along with any 
other active augmenting conditions via the internal bordering algorithm within Goma's 
nonlinear solver.


--------------
References
--------------

SAND2000-2465: Gates, I.D., Labreche, D. A., Hopkins, M. M. and Wilkes, E. D., 
2001. "Advanced Capabilities in GOMA 3.0 - Augmenting Conditions, Automatic 
Continuation, and Linear Stability Analysis," Sandia Technical Report.
**********************
AC (Volume Constraint)
**********************

::

    AC = VC <mat_id> <vc_type> <bc_id> <data_float_index> <species_id> <volume_value>

-----------------------
Description / Usage
-----------------------

This augmenting condition card allows the user to specify a boundary condition 
parameter as an unknown parameter in order to satisfy an integrated volumetric 
constraint. Three types of constraint are allowed: mesh volume, mass in mesh, and 
mass of species in mesh. All three integrated constraints are computed automatically 
and so no user-defined function is needed. See below for a mathematical description of 
each type of integrated constraint.

This card attaches a boundary condition parameter to these integrated constraints as an 
additional unknown degree of freedom. This parameter is specified in exactly the same 
way as detailed in the AC (Boundary Condition) section.

Definitions of the input parameters are as follows:

VC
    Mandatory string identifying the augmenting condition 
    card as being of type Volume Constraint.

<mat_id>
    An integer parameter giving the material identification
    number (element block id in the mesh) over which the 
    integrated volume constraint is to be applied. That is to 
    say, if more than one material occupies a problem 
    domain, a volume constraint can only be applied to the 
    space occupied by one of the materials for any given 
    augmenting condition.

<vc_type>
    An integer parameter that sets the type of volume 
    constraint as follows:

    vc_type = 1 (mesh volume):

    .. math::

        V_T = \int_V dV

    vc_type = 2 (mass in mesh):

    .. math::

        V_T = \int_V \rho dV

    where ρ is the fluid density.

    vc_type = 3 (mass of species in mesh):

    .. math::

        V_T = \int_V \rho y_i dV

    where :math:`y_i` is the mass fraction of the ith species.

<bc_id>
    An integer parameter giving the index of the boundary 
    condition card whose parameter is being used as the 
    new degree of freedom. Numbering begins with zero, 
    starting with the first boundary condition in the input 
    file and proceeding sequentially upward with each read 
    boundary condition card.

<data_float_index>
    An integer parameter that identifies the boundary 
    condition parameter that will be varied. It is an index 
    that starts at zero with leftmost float value on the 
    <bc_id> boundary condition card and increments 
    upward from right to left.

<species_id>
    An integer parameter that identifies the species number 
    to be used when evaluating a vc_type = 3 constraint 
    equation. For other types of constraints, it is still 
    necessary for the syntax, but is unused. Good procedure 
    sets it to zero.

<volume_value>
    A float parameter that fixes the value of V\ :sub:`T` used in the 
    expressions above. That is, it fixes the value that the 
    integrated quantity will have at the end of the 
    computation.

------------
Examples
------------

The following is an AC (Volume Constraint) card with its accompanying boundary 
condition card.

::

    AC = VC 1 1 2 1 0 3.14156

which has the following set of boundary conditions

::

    Number of BC = -1
    BC = U NS 1 0.0
    BC = V NS 1 1.0
    BC = CAPILLARY SS 10 1.0 10.0 0.0
    ...

This augmenting condition then applies to material 1 and will fix the mesh volume of 
this material to 3.14156. To do this it will vary the external pressure applied via the 
CAPILLARY card (the second float parameter). Note that the value given for this 
pressure in the input deck serves as a starting guess. At the end of the calculation, a new 
value for pressure will be in the memory location assigned to this float parameter. If an 
ASCII solution file is requested (via SOLN file =), the parameter value will be written 
following the output of the nodal unknown vector.

-------------------------
Technical Discussion
-------------------------

See the technical discussion appearing in the documentation for the AC (Boundary 
Condition) card.

--------------
References
--------------

SAND2000-2465: Gates, I.D., Labreche, D. A., Hopkins, M. M. and Wilkes, E. D., 
2001. "Advanced Capabilities in GOMA 3.0 - Augmenting Conditions, Automatic 
Continuation, and Linear Stability Analysis," Sandia Technical Report.
*************
GD_LINEAR
*************

::

	BC = GD_LINEAR SS <bc_id> <equation_name> <integer1> <variable_name> <integer2> <float1> <float2>

-----------------------
Description / Usage
-----------------------

**(PCC/VARIED)**

This boundary condition of type *Category 1* (see discussion) is used to impose a linear
function for any nodal variable, using the residual function form

.. math::

   C_1 + C_2x = 0

where :math:`C_1` and :math:`C_2` being the constant values and :math:`x` representing any variable
(<variable_name>). This boundary condition card can be used in combination with any
of the other *GD_** conditions as a building block to construct more complicated
conditions. Moreover, the resulting boundary condition can be applied as a strong
residual replacement to any differential equation type. Please see the examples on all of
these cards for details and instructive uses. Definitions of the input parameters are as
follows:

GD_LINEAR
    Name of the boundary condition (<bc_name>).
SS
    Type of boundary condition (<bc_type>), where **SS** denotes side set in
    the EXODUS II database.
<bc_id>
    The boundary flag identifier, an integer associated with <bc_type> that
    identifies the boundary location (side set in EXODUS II) in the problem
    domain.
<equation_name>
    A character string indicating the equation to which this boundary condition
    is applied (see the list of permissible values in the discussion above for
    *Category 1*).
<integer1>
    Species number of the mass transport equation. The value should be 0 unless
    the <equation_name> is of type R_MASS.
<variable_name>
    A character string indicating the variable that should be fixed (see the
    list of permissible values in the discussion above for *Category 1*).
<integer2>
    Species number of the concentration variable.The value should be 0 unless
    the <variable_name> is of type MASS_FRACTION.
<float1>
    Intercept :math:`C_1`
<float2>
    Slope :math:`C_2`

------------
Examples
------------

Following is a sample card:
::

	BC = GD_LINEAR SS 1 R_MESH1 0 MESH_POSITION1 0 -1. 2.

This boundary condition results in the equation :math:`2.0*x - 1.0 = 0` to be applied as a
boundary condition to the x-component of the mesh motion equation. x is the xcomponent
of the mesh position (N.B. not displacement, as *MESH_POSITION1* would
be replaced by *MESH_DISPLACEMENT1* in the above). The equation is actually
applied as a replacement to the x-component of the mesh motion equation and in this
case would lead to the mesh surface, defined by side set 1, to move or position itself
according to this linear relationship.

-------------------------
Technical Discussion
-------------------------

Note that this collocated boundary condition may be applied as a rotated, vector or
scalar condition depending on the equation to which this condition applies. Please
consult the example discussions on the other *GD_** options for more detailed
examples, as this boundary condition card can be used in an additive way with those.



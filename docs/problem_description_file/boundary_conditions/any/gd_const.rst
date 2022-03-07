************
GD_CONST
************

::

	BC = GD_CONST SS <bc_id> <equation_name> <integer1> <variable_name> <integer2> <float>

-----------------------
Description / Usage
-----------------------

**(PCC/VARIED)**

This boundary condition of type *Category 1* (see discussion) is used to impose a
constant value for any nodal variable, using the residual function form

.. math::

   x - C_1 = 0

:math:`C_1` being the constant value (<float>) and :math:`x` being the <variable_name>. This boundary
condition card can be used in combination with any of the other *GD_** conditions as a
building block to construct more complicated conditions. Please see the examples on
all of these cards for details and instructive uses. Definitions of the input parameters are
as follows:

GD_CONST
    Name of the boundary condition (<bc_name>).
SS
    Type of boundary condition (<bc_type>), where SS denotes side set in the
    EXODUS II database.
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
<float>
    Value of variable :math:`C_1`.

------------
Examples
------------

Following is a sample card:
::

	BC = GD_CONST SS 2 R_MESH_NORMAL 0 MASS_FRACTION 0 0.2

This boundary condition results in the equation :math:`C_1 - 0.2 = 0` being applied as a
boundary condition to the mesh-motion equation and being rotated into a normal-tangential
basis. :math:`C_1` is the concentration of the zeroth species. The equation is actually
applied as a replacement to the normal component of the mesh motion equation and in
this case would cause the mesh surface, defined by side set 2, to move as the
isoconcentration surface of :math:`C_1 = 0.2`.

-------------------------
Technical Discussion
-------------------------

Note that this collocated boundary condition may be applied as a rotated, vector or
scalar condition depending on the equation to which this condition applies. The
example above is a powerful demonstration of this boundary condition as a
distinguishing condition. Please consult the example discussions on the other *GD_**
options for more detailed examples, as this boundary condition card can be used in an
additive way with other *GD_** cards.

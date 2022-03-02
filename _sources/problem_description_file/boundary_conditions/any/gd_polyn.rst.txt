************
GD_POLYN
************

::

	BC = GD_POLYN SS <bc_id> <equation_name> <integer1> <variable_name> <integer2> <float_list>

-----------------------
Description / Usage
-----------------------

**(PCC/VARIED)**

This boundary condition of type *Category 1* (see discussion) is used to impose a
polynomial function for any nodal variable, using the residual function form of a 6\ :sup:`th` order
polynomial dependence on a variable

.. math::

   C_1 + C_2x + C_3x^2 + C_4x^3 + C_5x^4 + C_6x^5 + C_7x^6 = 0

There are three required and four optional parameters in the <float_list>; definitions of
the input parameters are as follows:

GD_POLYN
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
    Species number of the concentration variable. The value should be 0 unless
    the <variable_name> is of type MASS_FRACTION.
<float1>
    Intercept :math:`C_1`.
<float2>
    Slope :math:`C_2`.
<float3>
    Acceleration :math:`C_3`.
<float4>
    Coefficient for 3rd-order term :math:`C_4`.
<float5>
    Coefficient for 4th-order term :math:`C_5`.
<float6>
    Coefficient for 5th-order term :math:`C_6`.
<float7>
    Coefficient for 6th-order term :math:`C_7``.

------------
Examples
------------

Following is a set of sample cards:
::

	BC = GD_POLYN SS 2 R_ENERGY 0 MESH_POSITION1 0 {c1} {c2} {c3} {c4} {c5} {c6} {c7}

::

	BC = GD_LINEAR SS 2 R_ENERGY 0 TEMPERATURE 0   0. -1.

This boundary condition results in the equation

.. math::

   C_1 + C_2x + C_3x^2 + C_4x^3 + C_5x^4 + C_6x^5 + C_7x^6 = T

to be applied as a boundary condition on the energy equation, i.e., made a boundary
condition on temperature with second card, which brings in a dependence on
temperature. Here the coefficients are set by APREPRO, :math:`x` is the x-component of the
mesh position (N.B. not displacement, as *MESH_POSITION2* would be replaced by
*MESH_DISPLACEMENT2* in the above).

-------------------------
Technical Discussion
-------------------------

This condition is not used as often as *GD_LINEAR* and *GD_PARAB*, and in fact
supersedes those conditions. Please consult the example discussions on the other *GD_**
options and the example section after *GD_TABLE* for more descriptive examples.



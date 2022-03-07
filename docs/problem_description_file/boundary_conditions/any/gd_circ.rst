***********
GD_CIRC
***********

::

	BC = GD_CIRC SS <bc_id> <equation_name> <integer1> <variable_name> <integer2> <float1> <float2> <float3>

-----------------------
Description / Usage
-----------------------

**(PCC/VARIED)**

This boundary condition of type *Category 1* (see discussion) is used to impose a
quadratic function for any nodal variable using the residual function form

.. math::

   -{C_1}^2 + C_3 \left( x - C_2\right)^2 = 0.

where :math:`C_1`, :math:`C_2`, and :math:`C_3` are the constant values (<*float*>) and :math:`x` represents any variable
(<variable_name>). This boundary condition card can be used in combination with any
of the other *GD_** conditions as a building block to construct more complicated
conditions. *GD_CIRC* happens to be a convenient building block for circles or
elliptical functions (see examples below). Moreover, the resulting boundary condition
can be applied as a strong residual replacement to any differential equation type. Please
see the examples on all of these cards for details and instructive uses. Definitions of the
input parameters are as follows: (convenient for circles):

GD_CIRC
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
    is applied. See the list of permissible values in the introduction to the
    *Category 1 BCs* following the *Number of BC* card.
<integer1>
    Species number of the mass transport equation. The value should be 0 unless
    the <equation_name> is of type R_MASS.
<variable_name>
    A character string indicating the variable that should be fixed. See the
    list of permissible values in the introduction to the *Category 1 BCs*
    following the *Number of BC&* card.
<integer2>
    Species number of the concentration variable. The value should be 0 unless
    the <variable_name> is of type *MASS_FRACTION*.
<float1>
    Radius :math:`C_1`. This *should appear in only one* GD_CIRC *condition on
    each boundary*.
<float2>
    Origin :math:`C_2`.
<float3>
    Ellipticity :math:`C_3`.

------------
Examples
------------

Following is a sample set of cards:
::

	BC = GD_CIRC SS 1 R_MESH_NORMAL 0 MESH_POSITION1 0 1. 1. 1.
	BC = GD_CIRC SS 1 R_MESH_NORMAL 0 MESH_POSITION2 0 0. 1. 1.

This set of cards can be used to prescribe a mesh distinguishing condition for a mesh
surface with a quadratic dependence on x and y, a circle center at [1., 1.], and a radius
of 1.0 (note the radius only appears on one card).

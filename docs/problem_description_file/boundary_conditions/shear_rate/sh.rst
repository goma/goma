******
**SH**
******

::

	BC = SH NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/SHEAR_RATE)**

This boundary condition is used to set a Dirichlet condition for the scalar shear rate
unknown field.

Description of the input parameters is as follows:

========== ======================================================================
**SH**     Name of the boundary condition (<bc_name>).
**NS**     Type of boundary condition (<bc_type>), where **NS** denotes
           node set in the EXODUS II database.
<bc_id>    The boundary flag identifier, an integer associated with
           <bc_type> that identifies the boundary location (node set in
           EXODUS II) in the problem domain.
<float1>   Value at which the scalar shear rate unknown will be fixed
           on node set <bc_id>.
[float2]   An optional parameter (that serves as a flag to the code for a
           Dirichlet boundary condition). If a value is present, and is
           not -1.0, the condition is applied as a residual equation.
           Otherwise, it is a “hard set” condition and is eliminated
           from the matrix. *The residual method must be used when
           this Dirichlet boundary condition is used as a parameter in
           automatic continuation sequences*.
========== ======================================================================

------------
**Examples**
------------

An example of its used:
::

   BC = SH NS 10 0.5

This boundary condition sets the scalar shear rate unknown to 0.5 on nodeset 10.

-------------------------
**Technical Discussion**
-------------------------

The scalar shear rate unknown field is otherwise known as the second invariant of the
rate of deformation tensor.


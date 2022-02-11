*****
**F**
*****

::

	BC = F NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/FILL)**

This Dirichlet boundary condition specifies a value of the fill or level set unknown field
on a node set.

A description of the input parameters is as follows:

=========== ============================================================================
**F**       Name of the boundary condition (<bc_name>).
**NS**      Type of boundary condition (<bc_type>), where **NS** denotes
            node set in the EXODUS II database.
<bc_id>     The boundary flag identifier, an integer associated with
            <bc_type> that identifies the boundary location (node set in
            EXODUS II) in the problem domain.
<float1>    Value at which the fill or level set unknown will be fixed on
            this node set.
[float2]    An optional parameter (that serves as a flag to the code for a
            Dirichlet boundary condition). If a value is present, and is
            not -1.0, the condition is applied as a residual equation.
            Otherwise, it is a “hard set” condition and is eliminated
            from the matrix. *The residual method must be used when
            this Dirichlet boundary condition is used as a parameter in
            automatic continuation sequences*.
=========== ============================================================================

------------
**Examples**
------------

An example:
::

   BC = F NS 100 1.0

-------------------------
**Technical Discussion**
-------------------------

This boundary condition finds most of its use in the VOF/FILL interface tracking
algorithm where it is used to fix the value of the color function at an inlet or outlet
boundary. In the level set formulation, it is used less but is still useful in defining the
absolute fixed location of an interface by setting the value assigned to 0 on a node set.




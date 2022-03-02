***********
**SH_TENS**
***********

::

	BC = SH_TENS NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/shell_tension)**

This Dirichlet boundary condition specification is used to set a tension (in 
stress per unit length) to the inextensible shell equations (see EQ = 
shell_tension and EQ = shell_curvature) at an endpoint. This boundary
condition can be applied in two dimensions only, and only to the endpoint of a
bar-type element. In put is as follows :

=========== ===================================================================
**SH_TENS** Boundary condition name (<bc_name>) that defines the shell tension 
            (compressive or expansion depending on the sign). 
**NS**      Type of boundary condition (<bc_type>), where **NS** denotes
            node set in the EXODUS II database. Note that this must be
            a single-node node set representing and endpoint to a bar
            element type.
<bc_id>     The boundary flag identifier, an integer associated with
            <bc_type> that identifies the boundary location (node set in
            EXODUS II) in the problem domain.
<float1>    Value of tension.
[float2]    An optional parameter (that serves as a flag to the code for a
            Dirichlet boundary condition). If a value is present, and is
            not -1.0, the condition is applied as a residual equation.
            Otherwise, it is a "hard set" condition and is eliminated
            from the matrix. *The residual method must be used when
            this Dirichlet boundary condition is used as a parameter in
            automatic continuation sequences*.
=========== ===================================================================

------------
**Examples**
------------

The following is a sample card for applying a Dirichlet condition on the tension for a shell equation:
::

   BC = SH_TENS NS 100 10.

This condition sets a tension of 10.0 at Nodeset 100.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

GT-027.1: GOMA’s Shell Structure Capability: User Tutorial (GT-027.1)
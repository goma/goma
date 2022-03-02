***************
**SHELL_PARTC**
***************

::

	BC = SHELL_PARTC NS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(DC/R_SHELL_PARTC)**

This boundary condition card applies a particle volume fraction to the boundary of a shell-element sheet. The corresponding equation is EQ=shell_filmh. The boundary condition is applied to a node set.

=========== ============================================================
SHELL_PARTC Name of boundary condition.
NS          Type of boundary condition (<bc_type>), where NS
            denotes node set in the EXODUS II database.
<bc_id>     The boundary flag identifier, an integer associated with
            <bc_type> that identifies the boundary location (node
            set in EXODUS II) in the problem domain.
<float1>    shell_partc, the value of film thickness at the boundary.
[float2]    Optional floating point number set between 0.0 and 1.0
            which serves as a flag to the code for a Dirichlet
            boundary condition. If this value is present, and is not
            1.0, the condition is applied as a residual equation.
            Otherwise, it is “hard-set” condition and is eliminate
            from the matrix. The residual method must beused 
            when this Dirichlet boundary condition is used as a
            parameter in automatic continuation sequences.
=========== ============================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = SHELL_PARTC NS   100 0.

This condition applies a particles volume fractioin of 0.0 at nodeset 100.

-------------------------
**Technical Discussion**
-------------------------

The equation applied at the specified nodeset in place of the particles conservation equation.



--------------
**References**
--------------

No References.
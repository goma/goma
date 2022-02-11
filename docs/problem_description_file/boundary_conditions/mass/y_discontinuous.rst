*******************
**Y_DISCONTINUOUS** 
*******************

::

	BC = Y_DISCONTINUOUS NS <bc_id> <integer1> <float1> [float2 integer2]

-----------------------
**Description / Usage**
-----------------------

**(DC/MASS)**

This card is used to set a constant valued Dirichlet boundary condition for the species
unknown. The condition only applies to interphase mass, heat, and momentum transfer
problems applied to discontinuous (or multivalued) species unknown variables at an
interface, and it must be invoked on fields that employ the **Q1_D** or **Q2_D**
interpolation functions to “tie” together or constrain the extra degrees of freedom at 
the interface in question.

Definitions of the input parameters are as follows:

=================== ========================================================
**Y_DISCONTINUOUS** Name of the boundary condition (<bc_name>).
**NS**              Type of boundary condition (<bc_type>), where **NS**
                    denotes node set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (node
                    set in EXODUS II) in the problem domain.
<integer1>          Species subvariable number.
<float1>            Value of the species unknown on the boundary. Note,
                    the units depend on the specification of the type of the
                    species unknown.
[float2]            An optional parameter (that serves as a flag to the code
                    for a Dirichlet boundary condition). If a value is present,
                    and is not -1.0, the condition is applied as a residual
                    equation. Otherwise, it is a “hard set” condition and is
                    eliminated from the matrix. *The residual method must
                    be used when this Dirichlet boundary condition is used
                    as a parameter in automatic continuation sequences*.
[integer2]          Element block ID; only applicable to node sets, this
                    optional parameter specifies the element block on which
                    to impose the boundary condition, if there is a choice, as     
                    occurs at discontinuous variable interfaces where there
                    may be more that one unknown corresponding to
                    species 0 at a single node. This parameter allows the
                    user to specify which unknown to set the boundary
                    condition on, and allows for a jump discontinuity in
                    species value across a discontinuous variables interface.
=================== ========================================================

------------
**Examples**
------------

The following is a sample input card with no Dirichlet flag:
::

   BC = Y_DISCONTINUOUS SS 3   0   0.00126

-------------------------
**Technical Discussion**
-------------------------

Typically, this boundary condition may be used to set the species unknown variable on
one side of a discontinuous variables interface, while the species unknown variable on
the other side of the interface is solved for via a *KINEMATIC_SPECIES* boundary
condition. Note, this boundary condition is not covered by the test suite, and thus, 
may or may not work.




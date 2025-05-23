*****
**Y**
*****

::

	BC = Y NS <bc_id> <integer> <float1> [float2] [integer2]

-----------------------
**Description / Usage**
-----------------------

**(DC/MASS)**

This card is used to set the Dirichlet boundary condition of constant concentration 
for a given species. Definitions of the input parameters are as follows:

=========== ===============================================================
**Y**       Name of the boundary condition (<bc_name>).
**NS**      Type of boundary condition (<bc_type>), where **NS** denotes
            node set in the EXODUS II database.
<bc_id>     The boundary flag identifier, an integer associated with
            <bc_type> that identifies the boundary location (node set in
            EXODUS II) in the problem domain.
<integer1>  Species number of concentration.
<float1>    Value of concentration, in user’s choice of units, e.g. moles/
            :math:`cm^3`.
[float2]    An optional parameter (that serves as a flag to the code for a
            Dirichlet boundary condition). If a value is present, and is
            not -1.0, the condition is applied as a residual equation.
            Otherwise, it is a “hard set” condition and is eliminated
            from the matrix. *The residual method must be used when
            this Dirichlet boundary condition is used as a parameter in
            automatic continuation sequences*.
[integer2]  Element block ID; only applicable to node sets, this optional
            parameter specifies the element block on which to impose
            the boundary condition, if there is a choice, as occurs at
            discontinuous variable interfaces where there may be more
            that one unknown corresponding to species 0 at a single
            node. This parameter allows the user to specify which
            unknown to set the boundary condition on, and allows for a
            jump discontinuity in species value across a discontinuous
            variables interface.
=========== ===============================================================

------------
**Examples**
------------

The following is a sample card with no Dirichlet flag:
::

   BC = Y NS 3   0   0.00126

-------------------------
**Technical Discussion**
-------------------------

No Discussion.





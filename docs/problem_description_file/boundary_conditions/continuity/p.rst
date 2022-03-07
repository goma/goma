*****
**P**
*****

::

	BC = P NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/CONTINUITY)**

This Dirichlet boundary condition specification is used to set a constant pressure on a
node set. It is generally used for specifying a pressure datum on a single-node node set.
The pressure datum is useful for setting the absolute value of the pressure, which, for
many problems, is indeterminate to a constant. Pressure datums are especially
important for closed flow problems, such as the lid driven cavity, where there is no
inflow or outflow. Mass conservation problems can arise if this card is used to specify
the pressure along a group of nodes, since this equation replaces the continuity
equation. To specify pressure for a group of nodes, it is preferable to use the flow
pressure boundary condition, which is applied in a weak sense to the momentum
equation and does not cause mass conservation problems. Definitions of the input
parameters are as follows:

========== ================================================================
**P**      One-character boundary condition name (<bc_name>) that
           defines the pressure.
**NS**     Type of boundary condition (<bc_type>), where **NS** denotes
           side set in the EXODUS II database.
<bc_id>    The boundary flag identifier, an integer associated with
           <bc_type> that identifies the boundary location (node set in
           EXODUS II) in the problem domain.
<float1>   Value of pressure.
[float2]   An optional parameter (that serves as a flag to the code for a
           Dirichlet boundary condition). If a value is present, and is
           not -1.0, the condition is applied as a residual equation.
           Otherwise, it is a “hard set” condition and is eliminated
           from the matrix. *The residual method must be used when
           this Dirichlet boundary condition is used as a parameter in
           automatic continuation sequences*.
========== ================================================================

------------
**Examples**
------------

The following are sample cards for specifying a pressure Dirichlet card:
::

   BC = P NS 7   0.

::

   BC = P NS 7   0. 1.0

where the second form is an example using the “residual” method for applying the
same Dirichlet condition.

-------------------------
**Technical Discussion**
-------------------------

See the technical discussion for the *UVW* velocity for a discussion of the two ways of
applying Dirichlet boundary conditions.








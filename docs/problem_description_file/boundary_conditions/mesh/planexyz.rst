************
**PLANEXYZ**
************

::

	BC = {PLANEX | PLANEY | PLANEZ} SS <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(PCC/ MESH)**

This boundary condition card is used to specify a planar surface (solid) boundary
description as a replacement on the X, Y, or Z-component (*PLANEX, PLANEY,
PLANEZ*, respectively) of the mesh equations (see *EQ* cards *mesh1, mesh2*, or *mesh3*).
The form of this equation is given by

.. math::

   f(x, y, z) = ax + by + cz + d = 0

This mathematical form and its usage is exactly like the BC = PLANE boundary
condition card (see PLANE for description), but is applied to the mesh motion
equations without rotation. Definitions of the input parameters are given below; note
that <floatlist> has four parameters corresponding to the four constants in the equation:

============================== ==========================================================
**{PLANEX | PLANEY | PLANEZ}** Boundary condition name (<bc_name>) where:
                               
                               * **PLANEX** - normal predominantly in X direction
                               * **PLANEY** - normal predominantly in Y direction
                               * **PLANEZ** - normal predominantly in Z direction
**SS**                         Type of boundary condition (<bc_type>), where **SS**
                               denotes side set in the EXODUS II database
<bc_id>                        The boundary flag identifier, an integer associated with
                               <bc_type> that identifies the boundary location (side set
                               in EXODUS II) in the problem domain.
<float1>                       :math:`a` in function :math:`f(x, y, z)`
<float2>                       :math:`b` in function :math:`f(x, y, z)`
<float3>                       :math:`c` in function :math:`f(x, y, z)`
<float4>                       :math:`d` in function :math:`f(x, y, z)`
============================== ==========================================================

------------
**Examples**
------------

Following is a sample input card for a predominantly X-directed surface (viz, as planar
surface whose normal has a dominant component in the positive or negative X
direction):
::

     BC = PLANEX SS 101 1.0 1.0 -2.0   100.0

This boundary condition leads to the application of the equation :math:`1.0x + 1.0y – 2.0z = –100.0`
to the *mesh1* equation on EXODUS II side set number 101.

-------------------------
**Technical Discussion**
-------------------------

These conditions are sometimes used instead of the more general *PLANE* boundary
condition in situations where *ROTATION* (see *ROT* command section) leads to poor
convergence of the matrix solvers or is not desirable for some other reason. In general,
the *PLANE* condition should be used instead of these, but in special cases these can be
used to force the application of the planar geometry to a specific component of the
mesh stress equation residuals. Full understanding of the boundary rotation concept is
necessary to understand these reasons (see *Rotation Specifications*).



--------------
**References**
--------------

GT-001.4: GOMA and SEAMS tutorial for new users, February 18, 2002, P. R. Schunk
and D. A. Labreche

GT-007.2: Tutorial on droplet on incline problem, July 30, 1999, T. A. Baer

GT-013.2: Computations for slot coater edge section, October 10, 2002, T.A. Baer

GT-018.1: ROT card tutorial, January 22, 2001, T. A. Baer

.. 
	TODO - In line 20 the picture needs to be changed into the equation.

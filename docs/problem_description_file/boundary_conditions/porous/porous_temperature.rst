**********************
**POROUS_TEMPERATURE**
**********************

::

	BC = POROUS_TEMPERATURE NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/POR_TEMP)**

This Dirichlet boundary condition is used to set the temperature for a nonisothermal
porous media problem at a node set. It can be applied to a node set on a boundary of a
*POROUS_SATURATED, POROUS_UNSATURATED* or *POROUS_TWO_PHASE*
medium type (see *Media Type* card).

====================== ===============================================================
**POROUS_TEMPERATURE** Boundary condition name (bc_name).
**NS**                 Type of boundary condition (<bc_type>), where **NS**
                       denotes node set in the EXODUS II database.
<bc_id>                The boundary flag identifier, an integer associated with
                       <bc_type> that identifies the boundary location (node
                       set in EXODUS II) in the problem domain.
<float1>               Value of temperature at the NS in the porous medium.
[float2]               An optional parameter (that serves as a flag to the code
                       for a Dirichlet boundary condition). If a value is present,
                       and is not -1.0, the condition is applied as a residual
                       equation. Otherwise, it is a "hard set" condition and is
                       eliminated from the matrix. *The residual method must
                       be used when this Dirichlet boundary condition is used
                       as a parameter in automatic continuation sequences*.
====================== ===============================================================

------------
**Examples**
------------

An example input card for this boundary condition follows:
::

   BC = POROUS_TEMPERATURE NS 101   1.0   1.0

This card sets the temperature(*p_temp* in the output EXODUS II file) in element block
1 at the nodes defined by nodeset 101. Also, the second 1.0 float is to instruct goma to
apply this condition in a residual form.

-------------------------
**Technical Discussion**
-------------------------

This condition is used to set a temperature boundary condition for nonisothermal
porous media problems, viz. problems that use the R_POR_ENERGY equation (called
*EQ = porous_energy*). This energy equation is written in multiphase enthalpy form
and hence requires a different equatioin that for continuous media.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

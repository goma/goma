*******
**API**
*******

::

	BC = API NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/ACOUS_PIMAG)**

This Dirichlet boundary condition card is used to set constant amplitude of the
imaginary part of the acoustic pressure. Definitions of the input parameters are as follows:

**API**     
    Name of the boundary condition (<bc_name>).

**NS**      
    Type of boundary condition (<bc_type>), where **NS** denotes node set in
    the EXODUS II database.

<bc_id>
    The boundary flag identifier, an integer associated with <bc_type> that
    identifies the boundary location (node set in EXODUS II) in the problem
    domain.

<float1>
    Value of the imaginary part of the acoustic pressure amplitude.

[float2]
    An optional parameter (that serves as a flag to the code for a Dirichlet
    boundary condition). If a value is present, and is not -1.0, the condition
    is applied as a residual equation. Otherwise, it is a “hard set” condition
    and is eliminated from the matrix. *The residual method must be used when
    this Dirichlet boundary condition is used as a parameter in automatic
    continuation sequences*.

------------
**Examples**
------------

The following is a sample input card:
::

   BC = API NS 100 1000.0


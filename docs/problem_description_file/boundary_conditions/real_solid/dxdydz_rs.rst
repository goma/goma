*************
**DXDYDZ_RS**
*************

::

	BC = {DX_RS | DY_RS | DZ_RS} NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/REALSOLID)**

This boundary condition format is used to set a constant X, Y, or Z real-solid
displacement on the real-solid mesh motion equations (see *TOTAL_ALE* option on the
*Mesh Motion* card). Each such specification is made on a separate input card. These
boundary conditions are of the Dirichlet type and must be applied on EXODUS II node
sets. Definitions of the input parameters are as follows:

+---------------------------+----------------------------------------------------------------+
|{**DX_RS | DY_RS | DZ_RS**}| Boundary condition name (<bc_name>) that defines the           |
|                           | displacement, where:                                           |
|                           |                                                                |
|                           |   * **DX_RS** - real solid X displacement                      |
|                           |   * **DY_RS** - real solid Y displacement                      |
|                           |   * **DZ_RS** - real solid Z displacement                      |
+---------------------------+----------------------------------------------------------------+
|**NS**                     | Type of boundary condition (<bc_type>), where **NS** denotes   |
|                           | node set in the EXODUS II database.                            |
+---------------------------+----------------------------------------------------------------+
|<bc_id>                    | The boundary flag identifier, an integer associated with       |
|                           | <bc_type> that identifies the boundary location (node set in   |
|                           | EXODUS II) in the problem domain.                              |
+---------------------------+----------------------------------------------------------------+
|<float1>                   | Value of the real_solid displacement (X, Y, or Z) defined      |
|                           | above.                                                         |
+---------------------------+----------------------------------------------------------------+
|[float2]                   | An optional parameter (that serves as a flag to the code for a |
|                           | Dirichlet boundary condition). If a value is present, and is   |
|                           | not -1.0, the condition is applied as a residual equation.     |
|                           | Otherwise, it is a “hard set” condition and is eliminated      |
|                           | from the matrix. *The residual method must be used when        |
|                           | this Dirichlet boundary condition is used as a parameter in    |
|                           | automatic continuation sequences*.                             |
+---------------------------+----------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card which applies in an X-displacement boundary condition to
the real-solid to the nodes in node set 100, specifically an X- real-solid Displacement of
0.1. These displacements are applied immediately to the unknowns, and hence result in
immediate mesh displacement from the initial state.

::

     BC = DX_RS NS 100 1.0

This sample card applies the same condition above, but as a residual equation that is
iterated upon with Newton’s method.

::

     BC = DX_RS NS 100 1.0   1.0

The second float 1.0 forces this application. This approach is advisable in most
situations, as the nodes are gradually moved as a part of the mesh deformation process.
Sudden movements, as in the first example, can lead to folds in the mesh.

-------------------------
**Technical Discussion**
-------------------------

This condition performs the same function as *DX|DY|DZ* boundary conditions, except
that it is applied to the real-solid of a *TOTAL_ALE* solid mesh motion model (see *Mesh
Motion* card). More than likely, these conditions are applied together with geometry
conditions on the mesh equations, e.g. *PLANE, DX, DY, GEOM*, etc., on the same
boundary. *TOTAL_ALE* mesh motion involves two sets of elasticity equations: mesh
motion equations (*mesh1* and *mesh2*), and real-solid elasticity equations (*mom_solid1*
and *mom_solid2*).



--------------
**References**
--------------

GT-005.3: THE NEW TOTAL-ARBITRARY-LAGRANGIAN-EULERIAN (TALE)
CAPABILITY and its applicability to coating with/on deformable media, August 6,
1999, P. R. Schunk
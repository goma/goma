**********
**DXDYDZ**
**********

::

	BC = {DX | DY | DZ} NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/MESH)**

This boundary condition format is used to set a constant X, Y, or Z displacement. Each
such specification is made on a separate input card. These boundary conditions must be
applied to node sets. Definitions of the input parameters are as follows:

=================== =============================================================================
**{DX | DY | DZ}**  Two-character boundary condition name (<bc_name>) that
                    defines the displacement, where:

                    	* **DX** - X displacement
                    	* **DY** - Y displacement
                    	* **DZ** - Z displacement
**NS**              Type of boundary condition (<bc_type>), where **NS** denotes
                    node set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (node set in
                    EXODUS II) in the problem domain.
<float1>            Value of the displacement (X, Y, or Z) defined above.
[float2]            An optional parameter (that serves as a flag to the code for a
                    Dirichlet boundary condition). If a value is present, and is
                    not -1.0, the condition is applied as a residual equation. 
                    Otherwise, it is a “hard set” condition and is eliminated
                    from the matrix. *The residual method must be used when
                    this Dirichlet boundary condition is used as a parameter in
                    automatic continuation sequences.*
=================== =============================================================================

------------
**Examples**
------------

Following is a sample card which applies an X-displacement boundary condition to the
nodes in node set 100, specifically an X-Displacement of 1.0. These displacements are
applied immediately to the unknowns and hence result in immediate mesh
displacements from the initial state.

::

     BC = DX NS 100   1.0

This sample card applies the same condition as above, except as a residual equation
that is iterated upon with Newton’s method.

::

     BC = DX NS 100 1.0   1.0

The second float 1.0 forces this application. This approach is advisable in most
situations, as the nodes are gradually moved as a part of the mesh deformation process;
sudden movements, as in the first example, can lead to folds in the mesh.

-------------------------
**Technical Discussion**
-------------------------

Application of boundary conditions of the Dirichlet type on mesh motion requires
different considerations than those on non-mesh degrees of freedom. Sudden
displacements at a point, without any motion in the mesh surrounding that point, can
lead to poorly shaped elements. It is advisable to apply these sorts of boundary
conditions as residual equations, as discussed above. Examples of how these conditions
are used to move solid structures relative to a fluid, as in a roll-coating flow, are
contained in the references below.



--------------
**References**
--------------

GT-003.1: Roll coating templates and tutorial for GOMA and SEAMS, February 29,
2000, P. R. Schunk and M. S. Stay

GT-005.3: THE NEW TOTAL-ARBITRARY-LAGRANGIAN-EULERIAN (TALE)
CAPABILITY and its applicability to coating with/on deformable media, August 6,
1999, P. R. Schunk.


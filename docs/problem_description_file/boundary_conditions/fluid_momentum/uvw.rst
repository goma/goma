*******
**UVW**
*******

::

	BC = {U | V | W} NS <bc_id> <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

**(DC/MOMENTUM)**

This Dirichlet boundary condition specification is used to set a constant velocity in the
X-, Y-, or Z-direction. Each such specification is made on a separate input card.
Definitions of the input parameters are as follows:

+----------------+----------------------------------------------------------------+
|**{U | V | W}** | One-character boundary condition name (<bc_name>) that         |
|                | defines the velocity direction, where:                         |
|                |                                                                |
|                |   * **U** - Indicates X velocity component                     |
|                |   * **V** - Indicates Y velocity component                     |
|                |   * **W** - Indicates Z velocity component                     |
+----------------+----------------------------------------------------------------+
|**NS**          | Type of boundary condition (<bc_type>), where **NS**           |
|                | denotes node set in the EXODUS II database.                    |
+----------------+----------------------------------------------------------------+
|<bc_id>         | The boundary flag identifier, an integer associated with       |
|                | <bc_type> that identifies the boundary location                |
|                | (node set in EXODUS II) in the problem domain.                 |
+----------------+----------------------------------------------------------------+
|<float1>        | Value of velocity component.                                   |
+----------------+----------------------------------------------------------------+
|[float2]        | An optional parameter (that serves as a flag to the code for a |
|                | Dirichlet boundary condition). If a value is present, and is   |
|                | not -1.0, the condition is applied as a residual equation.     |
|                | Otherwise, it is a “hard set” condition and is eliminated      |
|                | from the matrix. *The residual method must be used when        |
|                | this Dirichlet boundary condition is used as a parameter in    |
|                | automatic continuation sequences*.                             |
+----------------+----------------------------------------------------------------+

------------
**Examples**
------------

The following are sample input cards for the X velocity component Dirichlet card:
::

     BC = U NS 7   1.50

::

     BC = U NS 7   1.50   1.0

where the second example uses the “residual” method for applying the same Dirichlet
condition.

-------------------------
**Technical Discussion**
-------------------------

This class of card is used to set Dirichlet conditions on the velocity components. When
the second optional float parameter is not present, the matrix rows corresponding to the
appropriate velocity component for nodes on this node set are filled with zeros, the
diagonal element is set to one, the corresponding residual entry is also set to zero, and
in the solution vector the appropriate degree of freedom is set to the value specified by
<float1>. This is the so-called “hard set” method for specifying Dirichlet conditions.

An alternate method for specifying Dirichlet conditions is applied when the second
float parameter is present (the actual value is not important except that it be different
from -1.0). In this case, the Dirichlet constraint is applied as a residual equation. That
is, the momentum equation for the appropriate component at each node in the nodeset
is replaced by the residual equation,

.. figure:: /figures/074_goma_physics.png
	:align: center
	:width: 90%

This residual equation is included in the Newton’s method iteration scheme like any
other residual equation. Note that in this case, nothing is set in the solution vector since
that will occur automatically as part of the iteration method.



--------------
**References**
--------------

No References.
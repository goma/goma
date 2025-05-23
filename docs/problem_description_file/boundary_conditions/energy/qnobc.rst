*******************
QNOBC
*******************

::

	BC = QNOBC SS <bc_id>

-----------------------
Description / Usage
-----------------------

**(WIC/ENERGY)**

This boundary condition card applies the open boundary condition equivalent
on the energy equation (see momentum open BC FLOW_STRESSNOBC)

+-------------------+------------------------------------------------------------+
|**QNOBC**          | Name of the boundary condition.                            |
+-------------------+------------------------------------------------------------+
|**SS**             | Type of boundary condition (<bc_type>), where **SS**       | 
|                   | denotes side set in the EXODUS II database.                |
+-------------------+------------------------------------------------------------+
|<bc_id>            | The boundary flag identifier, an integer associated with   |
|                   | <bc_type> that identifies the boundary location (side set  |
|                   | in EXODUS II) in the problem domain.                       |
+-------------------+------------------------------------------------------------+

------------
Examples
------------

Following is a sample card:
::

     BC = QNOBC SS 3

Here the BC is applied to SS 3

-------------------------
Technical Discussion
-------------------------

The finite element formulation of the energy equation generates
boundary integrals of the form:

.. math::
   
   \int_\Gamma \mathbf{n} \cdot \mathbf{q}\ dS

where :math:`\mathbf{q}` is the heat flux.

This boundary condition adds this term back in at the boundary avoiding the
natural boundary condition.

--------------
**References**
--------------

Griffiths, D.F., “The ‘no boundary condition’ outflow boundary condition,” IJNMF, 24,
393-411, (1997)

Sani, R.L., and P.M. Gresho, “Resume and remarks on the open boundary condition
minisymposium,” IJNMF, 18, 983-1008, (1994).

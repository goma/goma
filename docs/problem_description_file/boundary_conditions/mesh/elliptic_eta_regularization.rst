***************************
ELLIPTIC_ETA_REGULARIZATION
***************************

::

	BC = ELLIPTIC_ETA_REGULARIZATION SS <bc_id> <float1>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MESH2)**

This Boundary condition applies a regularization term on the boundary for spacing in the :math:`\eta` direction

ELLIPTIC_ETA_REGULARIZATION
   Name of the boundary condition (<bc_name>).
SS 
   Type of boundary condition (<bc_type>), where **SS** denotes
   side set in the EXODUS II database.
<bc_id>
   The boundary flag identifier, an integer associated with
   <bc_type> that identifies the boundary location (side set in
   EXODUS II) in the problem domain.
<float1> 
   Penalty on the regularization term (float).

------------
**Examples**
------------

The following sample card
::

     BC = ELLIPTIC_ETA_REGULARIZATION SS 7 100.0


-------------------------
**Technical Discussion**
-------------------------

.. math::

    \mathbf{n} \cdot \nabla \eta = -M * log(x_\eta^2 + y_\eta^2)

Here :math:`M` is the penalty from the input card.


--------
**FAQs**
--------

--------------
**References**
--------------

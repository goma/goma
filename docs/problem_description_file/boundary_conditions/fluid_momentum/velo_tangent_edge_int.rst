*************************
**VELO_TANGENT_EDGE_INT**
*************************

::

	BC = VELO_TANGENT_EDGE_INT SS <bc_id1> <bc_id2> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC-EDGE/ROTATED MOMENTUM)**

This boundary condition card is used to make the velocity component tangent to the
contact line in the plane of the web equal to the component of web velocity
(:math:`W_x`, :math:`W_y`, :math:`W_z`) along the contact line. It imposes the identical constraint as the
*VELO_TANGENT_EDGE* card, but applies it as a strongly integrated condition rather
than a point collocated condition.

Definitions of the input parameters are as follows:

========================= =============================================================
**VELO_TANGENT_EDGE_INT** Name of the boundary condition.
**SS**                    Type of boundary condition (<bc_type>), where **SS**
                          denotes side set in the EXODUS II database.
<bc_id1>                  The boundary flag identifier, an integer associated with
                          <bc_type> that identifies the boundary location (side set
                          in EXODUS II) of the primary side set defining the edge
                          geometry in the problem domain. When applied to
                          dynamic contact lines, this side set should correspond to
                          the free surface.
<bc_id2>                  The boundary flag identifier, an integer associated with
                          <bc_type> that identifies the boundary location (side set
                          in EXODUS II) of the secondary side set defining the
                          edge geometry in the problem domain. The boundary
                          condition is applied to the curve defined as the
                          intersection of this side set with the primary side set
                          When applied to dynamic contact lines, this side set
                          should correspond to the substrate.
<float1>                  :math:`W_x`, x-component of the substrate (or web) velocity.
<float2>                  :math:`W_y`, y-component of the substrate (or web) velocity.
<float3>                  :math:`W_z`, z-component of the substrate (or web) velocity.
========================= =============================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = VELO_TANGENT_EDGE_INT SS 5 4   -1.0 0.0 0.0

This card imposes a tangent velocity component along the curve formed by the
intersection of sidesets 5 and 4. The value of the component is the projection of the
substrate velocity (-1.0, 0. ,0.) into the tangent direction. The tangent direction is along
the curve itself.

-------------------------
**Technical Discussion**
-------------------------

* This equation imposes the following constraint as a point collocated condition at
  the integration points of the elements along the curve:

.. figure:: /figures/087_goma_physics.png
	:align: center
	:width: 90%

where :math:`t_{cl}` is a vector tangent to the curve, *v* is the fluid velocity, *W* is the (constant)
velocity of the moving substrate, :math:`\phi_i` is the shape function each node along the curve
C. This integral condition is imposed strongly at each node. The reader is referred
to the sketch appearing with the *VELO_NORMAL_EDGE* card for a depiction of
these vectors.

* The reader is referred to the *VELO_TANGENT_EDGE* discussion for information
  about the context in which this condition is applied. Because it is applied in a
  different fashion than the former condition, it sometimes is the case that it will
  allow more flexibility in situations involving many boundary conditions applied in
  close proximity. There may also be situations where an integrated constraint
  results in better matrix conditioning that a collocated constraint.




.. TODO - In line 67, the photo needs to be repalced by the proper equation.
************************
**VELO_NORMAL_EDGE_INT**
************************

::

	BC = VELO_NORMAL_EDGE_INT SS <bc_id1> <bc_id2> <float>

-----------------------
**Description / Usage**
-----------------------

**(SIC-EDGE/ROTATED MOMENTUM)**

This boundary condition card is used to specify the normal velocity component on a
dynamic contact line in three-dimensions. The velocity component is normal to the
contact line in the plane of the web and is equal to :math:`V_n`. The free-surface side set should
always be <bc_id1>, the primary side set, and the web side set should be <bc_id2>, the
secondary side set. This boundary condition is identical in function to
*VELO_NORMAL_EDGE*. It differs only in that is applied as a strongly integrated
condition along the curve defined by <bc_id1> and <bc_id2>.

Definitions of the input parameters are as follows:

========================= ==========================================================
**VELO_NORMAL_EDGE_INT**  Name of the boundary condition.
**SS**                    Type of boundary condition (<bc_type>), where **SS**
                          denotes side set in the EXODUS II database.
<bc_id1>                  The boundary flag identifier, an integer associated with
                          <bc_type> that identifies the boundary location (side set
                          in EXODUS II) for the primary side set in the problem
                          domain. This side set should also be the side set
                          associated with the capillary free surface if used in the
                          context of a dynamic contact line.
<bc_id2>                  The boundary flag identifier, an integer associated with
                          <bc_type> that identifies the boundary location (side set
                          in EXODUS II) for the secondary side set defining the
                          edge in the problem domain. Together with <bc_id1>,
                          this defines the curve on which the boundary condition
                          applies as the intersection of the two side sets. In
                          problems involving dynamic contact lines, this side set
                          should correspond to the moving substrate.
<float>                   :math:`V_n`, a parameter supplying the imposed normal 
                          velocity
                          component value. This component is taken normal to
                          the edge curve parallel to <bc_id2>. See below for a
                          more detailed description.
========================= ==========================================================

------------
**Examples**
------------

The following is a sample card:
::

     BC = VELO_NORMAL_EDGE_INT SS 5 4   0.0

This card sets the normal-to-contact line component of the velocity to zero along the
curve defined by the intersections of side set 5 and 4.

-------------------------
**Technical Discussion**
-------------------------

* This boundary condition imposes a strongly
  integrated constraint of the form:

.. figure:: /figures/081_goma_physics.png
	:align: center
	:width: 90%

where :math:`\phi_i` is the velocity trial function, *v* is the
fluid velocity, :math:`v_m` is the mesh velocity and :math:`n_cl`
is the normal to the contact line in the plane of
the moving substrate <bc_id2>. The sketch at right depicts the orientation of this
latter vector.

.. figure:: /figures/082_goma_physics.png
	:align: center
	:width: 90%

* As noted above, this boundary condition functions nearly identically to the
  *VELO_NORMAL_EDGE* condition (except for its manner of application within
  *Goma*) and all comments appearing for the latter apply equally well for this
  boundary condition.



--------------
**References**
--------------

No References.
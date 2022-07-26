********************
**VELO_NORM_COLLOC**
********************

::

	BC = VELO_NORM_COLLOC SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MOMENTUM)**

This boundary condition allows the user to set the outward velocity component normal
to a surface. It is identical in function to the *VELO_NORMAL* boundary condition, but
differs in that it is applied as a point collocated condition.

Definitions of the input parameters are as follows:

===================== ==========================================================
**VELO_NORM_COLLOC**  Boundary condition designation.
**SS**                Type of boundary condition (<bc_type>), where **SS**
                      denotes side set in the EXODUS II database.
<bc_id>               The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location 
                      (side set in EXODUS II) in the problem domain.
<float>               :math:`v_n`, value of normal velocity component. Note that this
                      velocity component is relative to the motion of the
                      underlying mesh.
===================== ==========================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = VELO_NORM_COLLOC SS 20 0.0

This boundary condition will enforce an impenetrability constraint over side set 20 
as it
excludes normal velocity of the fluid relative to the mesh. This is by far the most
common context for this boundary condition.

-------------------------
**Technical Discussion**
-------------------------

* The actual equation that is applied to a node, *j*, on the surface in question is 
  as follows:

.. math::

  n \cdot \left(v_j - v_s\right) = v_n

  

where :math:`v_j` is the fluid velocity at the node, *n* the outward-pointing normal to the
surface, :math:`v_s` the velocity of the underlying mesh at the node, and :math:`v_n` is the normal
velocity set by <float> above.

* This constraint is a rotated collocated equation so that it will replace one of 
  the
  rotated components of the fluid momentum equation. This component should
  generally always be the normal rotated component. In two dimensions, this
  replacement is automatic. In three dimensions, this replacement must be specified
  by a *ROT* condition.

* As noted above this boundary condition applies exactly the same constraint as the
  *VELO_NORMAL* condition but via a point collocated method instead of as a
  strongly integrated condition. This might be advantageous at times when it is
  desirable to enforce a normal velocity component unambiguously at a point in the
  mesh.






******************
**VELO_NORMAL_LS**
******************

::

	BC = VELO_NORMAL_LS SS <bc_id> 0.0 <blk_id> <float1> <float2>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MOMENTUM)**

This boundary condition relaxes the VELO_NORMAL condition in the light phase of a
level-set simulation, thereby allowing gas to escape from a confined space.

Definitions of the input parameters are as follows:

=================== ===========================================================
**VELO_NORMAL_LS**  Boundary condition designation.
**SS**              Type of boundary condition (<bc_type>), where **SS**
                    denotes side set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (side set
                    in EXODUS II) in the problem domain.
<*blk_id*>          *blk_id*, an optional parameter that is the element block
                    number in conjugate problems that identifies the
                    material region where the *VELO_NORMAL_LS*
                    condition will be applied (usually the liquid element
                    block in solid/liquid conjugate problems). For external
                    boundaries, this optional parameter can be set to unity to
                    force the condition to be kept at a corner between two
                    side sets (2D only). This is handy for corner conditions.
                    Please see GTM-004.0 for details.1
<float>             L=interface half-width over which the
                    VELO_NORMAL bc changes.
<float2>            alpha=shift in the VELO_NORMAL change relative to
                    the LS interface. With alpha=0, VELO_NORMAL
                    begins to be enforced when the LS interface reaches a
                    distance L from a wall. With alpha=1,
                    VELO_NORMAL begins to be enforced when the LS
                    inteface reaches the wall.
=================== ===========================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = VELO_NORMAL_LS SS 10   0.0 {blk_id=1} 0.05 0.4.

-------------------------
**Technical Discussion**
-------------------------

The technical discussion under VELO_NORMAL largely applies here as well.



--------------
**References**
--------------

GT-001.4: GOMA and SEAMS tutorial for new users, February 18, 2002, P. R. Schunk
and D. A. Labreche

GTM-004.1: Corners and Outflow Boundary Conditions in Goma, April 24, 2001, P. R.
Schunk
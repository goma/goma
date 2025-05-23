********************
**LUBP_SH_FP_MATCH**
********************

::

	C = LUBP_SH_FP_MATCH SS <bc_id> <int1> <int2>

-----------------------
**Description / Usage**
-----------------------

**(STRONG_INT_SURF/R_LUBP)**

This boundary condition card matches the pressure in one region of confined flow (lubp) to the pressure from a second region of film flow (shell_filmp). The pressure matching is handled as a sideset between two shell regions. In this way both equations can be coupled for exit or entrance flows. The boundary condition is applied in collocated form, and replaces the R_LUBP equation.

=============== ==================================================
LUBP_SH_FP_FLUX Name of boundary condition.
SS              Type of boundary condition (<bc_type>), where SS
                denotes node set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set
                in EXODUS II) in the problem domain.
<int1>          Block id of mesh material which invokes the lubp
                equation.
<int2>          Block id of mesh material which invokes the
                shell_filmp equation.
=============== ==================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = LUBP_SH_FP_MATCH SS   100 2 1

This condition applies the matching tie condition at a side set boundary between block 2 (which invokes the EQ = lubp equation) and block 1 (which invokes the EQ=lubp equation).

-------------------------
**Technical Discussion**
-------------------------

The best example of the use of this equation is the exit of a metered coating flow. It must be used together with a flux-matching condition LUBP_SH_FP_FLUX.


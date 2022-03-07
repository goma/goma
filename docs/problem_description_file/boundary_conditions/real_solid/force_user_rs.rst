*****************
**FORCE_USER_RS**
*****************

::

	BC = FORCE_USER_RS SS <bc_id> <float1>...<floatn>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR REALSOLID)**

This boundary condition card applies a user-defined force per unit area (traction) on a
*TOTAL_ALE* real solid region (see *Mesh Motion* card). It differs from its counterpart
*FORCE_USER* only in the type of material to which the force is applied, as described
on the *Mesh Motion* card. The functional form of the force is programmed in the
function force_user_surf in bc_user.c, and can be made a function of any of
the independent or dependent variables of the problem, including position (see example
below). The force per unit area is applied to boundary delineated by the side set ID.
Definitions of the input parameters are as follows:

==================== =============================================================
**FORCE_USER_RS**    Name of the boundary condition (<bc_name>)
**SS**               Type of boundary condition (<bc_type>), where **SS**
                     denotes side set in the EXODUS II database.
<bc_id>              The boundary flag identifier, an integer associated with
                     <bc_type> that identifies the boundary location (side set
                     in EXODUS II) in the problem domain.
<float1>...<floatn>  Parameters list (length arbitrary) for parameterizing the
                     user defined force. These parameters are accessed
                     through the p[]array in force_user_surf.
==================== =============================================================

------------
**Examples**
------------

The input card
::

     BC = FORCE_USER_RS SS 3   {delta_t} 0.   1000.0 0.

used in conjuction with the following snippet of code in force_user_surf:

::

     /* Comment this out FIRST!!!!! */
     /* EH(-1,"No FORCE_USER model implemented"); */
     /**************************** EXECUTION BEGINS
     *******************************/
        if (time <= p[0])
           {
              func[0] = p[1]*time/p[0];
                func[1] = p[2]*time/p[0];
                func[2] = p[3]*time/p[0];
           }
       else
           {
                func[0] = p[1];
                func[1] = p[2];
                func[2] = p[3];
           }

applies a time-dependent force ramped from zero to 1000.0 in the +y direction over the
time period {*delta_t*}. Note how p[0] is the time period, viz. {*delta_t*}, over which the
force is ramped up.

-------------------------
**Technical Discussion**
-------------------------

Used commonly to apply a force per unit area to an external surface of a solid region
(*TOTAL_ALE* type, cf. *FORCE_USER*), that is nonconstant, viz. time varying or
spatially varying. The *FORCE_RS* and *NORM_FORCE_RS* boundary conditions can
be used for constant forces. This condition is applied as a weak integrated condition in
Goma, and hence will be additive with others of its kind.



--------------
**References**
--------------

No References.
***************
**VAR_CA_USER**
***************

::

	BC = VAR_CA_USER SS <bc_id1> <bc_id2> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC-EDGE/ROTATED MESH)**

This card is used to set a variable contact angle on a dynamic three-dimensional contact
line. It is identical in function to the *VAR_CA_USER* except that it allows the user to
provide a contact angle model to relate local contact angle to local capillary number.

Definitions of the input parameters are as follows:

================= =========================================================
**VAR_CA_USER**   Name of the boundary condition.
**SS**            Type of boundary condition (<bc_type>), where **SS**
                  denotes side set in the EXODUS II database.
<bc_id1>          The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set
                  in EXODUS II) in the problem domain. This identifies
                  the *primary side set*; it should be a free surface.
<bc_id2>          The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set
                  in EXODUS II) in the problem domain. This identifies
                  the *secondary side set*, which should be a “fixed”
                  geometric entity, e.g. PLANE or SPLINE. Taken
                  together, the primary and secondary sidesets define an
                  edge over which this boundary is applicable.
<float1>          W\ :sub:`x`, x-component of the substrate velocity.
<float2>          W\ :sub:`y`, y-component of the substrate velocity.
<float3>          W\ :sub:`z`, z-component of the substrate velocity.
[float4-floatn]   An optional list of floats which will be passed to the
                  user-supplied function for use with the user model.
================= =========================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = VAR_CA_USER SS 60 20   -1. 0. 0. 1.e-3 135.0

This card sets a variable contact angle condition on the edge between side sets 60 and
20. The solid substrate is moving at the fixed velocity (-1., 0., 0.). The var_CA_user
function is passed the constants 1.e-3 and 135.0 in variable locations p[0] and p[1],
respectively.

-------------------------
**Technical Discussion**
-------------------------

* *VAR_CA_USER* function is identical to *VAR_CA_EDGE*. It is applied to threedimensional
  dynamic contact lines in order to set a variable contact angle. The user
  must supply internal coding for the function var_CA_user in the file *user*_bc.c.
  This function receives as parameters the local capillary number as described under
  *VAR_CA_EDGE* and a double array containing the optional list of float
  parameters. It should return the cosine of the desired contact angle.

* What follows is an example that implements the linear contact angle model
  described in *VAR_CA_EDGE*.

::

     double
     var_CA_user(double Ca_local,
                int num,
                const double *a,
                double *d_cos_CA_Ca_local)
    {
     double cos_CA;
     double static_CA;
     double cT;
     static_CA = a[0]*M_PIE/180.0;
     cT = a[1];
     cos_CA = cos(static_CA) - cT * Ca_local;
     *d_cos_CA_Ca_local = cT;
     return ( cos_CA );
    }



--------------
**References**
--------------

No References.
*****************
**PRESSURE_USER**
*****************

::

	BC = PRESSURE_USER SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition allows the user to specify an arbitrary functional form for the
pressure field on a boundary via a user-defined subroutine. The boundary condition is
identical in form to the *FLOW_PRESSURE* and *FLOW_HYDROSTATIC* conditions,
but whereas the latter conditions have constant and linear spatial dependencies for the
pressure, this boundary condition allows for any dependency, including dependencies
on other degrees of freedom and time.

Definitions of the input parameters are as follows:

================== ===========================================================
**PRESSURE_USER**  Name of the boundary condition.
**SS**             Type of boundary condition (<bc_type>), where **SS** 
                   denotes side set in the EXODUS II database.
<bc_id>            The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (side set
                   in EXODUS II) in the problem domain.
<float_list>       A list of float values separated by spaces which will be
                   passed to the user-defined subroutine so the user can
                   vary the parameters of the boundary condition. This list
                   of float values is passed as a one-dimensional double
                   array to the appropriate C function.
================== ===========================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = PRESSURE_USER SS 10 0.0 3.0 0.5

-------------------------
**Technical Discussion**
-------------------------

* Frequently, it is desired to be able to set a pressure on a boundary that is more
  complicated than constant or linear; this boundary condition is used for this
  purpose. By modifying a function in user_bc.c (fn_dot_T_user), any
  functional dependence of pressure can be installed. This dependence may entail a
  more complicated spatial dependence, variability in time, and/or dependence on
  other degrees of freedom.

* An example is supplied in fn_dot_T_user that illustrates how this boundary
  condition can be used to set a sinusoidal-type of spatial dependence. A similar
  function could be used to set a temporal sinusoidal variation. The only caveat is
  that when inserting a function, it is very important that the sensitivities of the
  function with respect to position (and other degrees of freedom if they exist) be
  added to the array d_func. This does not apply to the time variable however.

* Like *FLOW_PRESSURE* and *FLOW_HYDROSTATIC*, this boundary condition is
  a weakly integrated condition. Therefore, it is additive with other weak conditions,
  but is superseded by strong conditions or Dirichlet conditions.







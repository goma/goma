**************
**FILL_INLET**
**************

::

	BC = FILL_INLET SS <bc_id> <float1>

-----------------------
**Description / Usage**
-----------------------

**(SPECIAL/FILL)**

This boundary condition allows the user to specify a value on a inlet boundary from
VOF problems employing discontinuous interpolation of the color function, *F*.

Description of the input parameters is as follows:

============== ===================================================================
**FILL_INLET** Name of the boundary condition.
**SS**         Type of boundary condition (<bc_type>), where **SS** denotes
               side set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set in
               EXODUS II) in the problem domain.
<float1>       The value of the fill function, *F*, as it flows across <bc_id>
               into the domain.
============== ===================================================================

------------
**Examples**
------------

An example:
::

   BC = FILL_INLET SS 10 1.0

-------------------------
**Technical Discussion**
-------------------------

* This boundary condition is useful only in problems involving VOF interface
  tracking in which the fill function is interpolated discontinuously. In this
  formulation, communication of the fill function value can only be made by finding
  the value of the fill function in the element upstream of the current position. While
  this is a stable formulation for the advective VOF method, it does introduce the
  complexity of determining which element is actually upstream.

* When there is no element upstream, as in the case of an inlet boundary, this
  boundary condition must be present to establish the value of the fill function that is
  flowing across the inlet boundary into the domain. Consequently, this boundary
  condition should be present on all inlet boundaries of the problem. It sometimes is
  also useful to have it on outflow boundaries as well, just in case a backflow
  situation arises.


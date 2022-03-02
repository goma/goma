************
**FLOWRATE**
************

::

	BC = FLOWRATE SS <bc_id> <float> <float | char_string>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition allows the user to specify a single value for the volumetric
flowrate of material across an inflow (or outflow) boundary. The pressure and velocity
fields on this boundary are then computed implicitly by *Goma*.

Definitions of the input parameters are as follows:

+-------------------+-----------------------------------------------------------------+
|**FLOWRATE**       | Boundary condition name.                                        |
+-------------------+-----------------------------------------------------------------+
|**SS**             | Type of boundary condition (<bc_type>), where **SS**            |
|                   | denotes side set in the EXODUS II database.                     |
+-------------------+-----------------------------------------------------------------+
|<bc_id>            | The boundary flag identifier, an integer associated with        |
|                   | <bc_type> that identifies the boundary location (side set       |
|                   | in EXODUS II) in the problem domain.                            |
+-------------------+-----------------------------------------------------------------+
|<float>            | *flowrate*, a parameter fixing the value of volumetric          |
|                   | flowrate across this boundary. For two-dimensional              |
|                   | *CARTESIAN* simulations, this value should be per unit          |
|                   | distance in the out-of-plane coordinate. For                    |
|                   | *CYLINDRICAL* and *SWIRLING* coordinate systems, this           |
|                   | flowrate value should include integration in the                |
|                   | azimuthal direction.                                            |
+-------------------+-----------------------------------------------------------------+
|<float|char_string>| This parameter can either be a <float> or a <char_string>.      |
+-------------------+------------+----------------------------------------------------+
|                   | float      |:math:`P_{guess}`, an initial guess for the pressure|
|                   |            |on the inlet                                        |
+-------------------+------------+----------------------------------------------------+
|                   | char_string|*read*, indicating that the initial guess           |
|                   |            |for the pressure file should be read                |
|                   |            |from the ASCII file identified on the               |
|                   |            |*GUESS file* card.                                  |
+-------------------+------------+----------------------------------------------------+

------------
**Examples**
------------

Specifying the average velocity on the inlet to a tube of radius 1.0:
::

     BC = FLOWRATE SS 10 3.1415 10.0

Since the radius is 1.0, the area of the surface is 3.1415 so the volumetric flowrate must
be specified as shown. An initial pressure guess of 10.0 is also supplied. Note this does
not specify the pressure on the boundary as the final value will generally be different
than specified here.

Continuing in the flowrate by reading the last pressure value from *GUESS file*:

::

     BC = FLOWRATE SS 10 3.2 read

-------------------------
**Technical Discussion**
-------------------------

* The requirement that is imposed by this boundary condition is the following
  integral:

.. figure:: /figures/099_goma_physics.png
	:align: center
	:width: 90%

where *U* is the flowrate value supplied on the card. It is imposed by the addition of
a Lagrange multiplier unknown on the boundary in question which will be
determined as a part of the solution process. For Newtonian and generalized
Newtonian models, the value of the multiplier is the inverse of the pressure value
on the boundary. Thus, a boundary condition nearly identical to a
*FLOW_PRESSURE* condition is applied to the sideset, but it takes as its pressure
the value of the inverse of the Lagrange multiplier unknown as it is computed.

The augmenting condition capability in *Goma* is used to impose the above integral.
When the boundary condition is invoked, an augmenting condition of the
appropriate type is automatically created. Its associated degree of freedom is the
Lagrange multiplier. During the iteration sequence, the user will see updates and
residuals for this augmenting condition.

* Originally, the initial guessed value for the pressure over the side set is read from
  the float value specified on this card, or from the *GUESS file* (if the parameter *read*
  is specified on this card). However, it can also be read from an EXODUS II
  database file. This is the same file the rest of the solution vector is read from if the
  problem is being restarted from a previous computation. If a value for the
  augmenting condition is present in this EXODUS II file, it will be read in. This
  value will override the float value specified on this card. The initial guess may still
  be read from the ASCII *GUESS file* by specifying *read* on the *Initial Guess* card
  and on the *Augmenting Conditions Initial Guess* card.



--------------
**References**
--------------

No References.
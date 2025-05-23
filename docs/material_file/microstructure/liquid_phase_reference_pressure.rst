***********************************
**Liquid phase reference pressure**
***********************************

::

   Liquid phase reference pressure = CONSTANT <float> [M/L-t2] or [N/L2]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model and model parameters for the liquid-phase
compressibility reference pressure. See *Liquid phase compressibility* card for
discussion and theory.

+-----------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**                 |model for the reference pressure, currently the only available option. It requires a |
|                             |single floating point value:                                                         |
|                             |                                                                                     |
|                             | * <float> - The reference pressure, in units of pressure.                           |
+-----------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The cards

::

   Liquid phase compressibility = CONSTANT {beta_liquid}

::

   Liquid phase reference pressure = CONSTANT {p_not}

leads to the application of a linearized compressibility model for the density of liquid in
the time-derivative capacitance term. This is useful for rigid porous media when the
conditions are such that the saturation front is sharp.

-------------------------
**Technical Discussion**
-------------------------

See discussion on *Liquid phase compressibility* card.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s capabilities for partially saturated flow in porous media,
September 1, 2002, P. R. Schunk

****************
**porous_unsat**
****************

::

	EQ = porous_unsat {Galerkin_wt} P_LIQ {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This equation cannot be invoked in an element block in which the media type is set to
POROUS_TWO_PHASE (cf. *Microstructure Properties, Media Type* card).
Otherwise, it is used exactly as the *porous_liq* equation card; please consult that section
for a detailed discussion.

See *porous_liq* card for description of input requirements.

------------
**Examples**
------------

See *porous_liq* card.

-------------------------
**Technical Discussion**
-------------------------

This card is used for single phase (viz. constant gas pressure) simulations of partially
saturated flow, as described by the *Media Type* material property card. The equation it
invokes is one of Darcy flow in a partially saturated medium in which the gas phase
pressure is taken as constant. The dependent variable here is the liquid phase pressure.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk
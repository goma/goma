***********************************
**Porous Shell Cross Permeability**
***********************************

::

   Porous Shell Cross Permeability = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card is used to set the permeability in the thin direction of a shell porous region.
The property is used for the porous_sat_open equation. The in-shell (in-plane for
a flat shell) permeabilities are set on the Permeability card. Please consult the
references for the equation form. The property can take on one of two models:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model applies a constant cross-region permeability. It requires a single        |
|                          |floating point input:                                                                |
|                          |                                                                                     |
|                          | * <float1> is the cross region permeability                                         |
+--------------------------+-------------------------------------------------------------------------------------+
|**EXTERNAL_FIELD**        |This model is used to read a finite element mesh field representing the cross-term   |
|                          |permeability. Please consult tutorials listed below for proper usage. This model     |
|                          |requires one float:                                                                  |
|                          |                                                                                     |
|                          | * <float1> scale factor for incoming exodusII field and desired level.              |
|                          |                                                                                     |
|                          |The ExodusII field variable name should be “CROSS_PERM”, viz.                        |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

No Examples.




--------------
**References**
--------------

S. A. Roberts and P. R. Schunk 2012. in preparation.

Randy Schunk 2011. GT-038 “Pixel-to-Mesh-Map Tool Tutorial for GOMA”. Memo to
distribution.
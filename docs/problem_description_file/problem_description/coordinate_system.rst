*********************
**Coordinate System**
*********************

::

	Coordinate System = {char_string}

-----------------------
**Description / Usage**
-----------------------

This card is required for each material section in the *Problem Description File*. It is used to specify formulation of the equations to be solved. Valid options for {char_string} are as follows:

======================= ==================================================
**CARTESIAN**           For a two (x-y) or three (x-y-z) dimensional Cartesian
                        formulation.
**CYLINDRICAL**         For an axisymmetric (z-r) or three-dimensional
                        cylindrical (z-r-:math:`\theta`) formulation; the three-dimensional option has not been tested.
**SPHERICAL**           For a spherical (r-:math:`\theta`-:math:`\phi`) 
                        formulation.
**SWIRLING**            For a two-dimensional formulation (z-r-:math:`\theta`) 
                        with a swirling velocity component that is independent of azimuthal coordinate.
**PROJECTED_CARTESIAN** For use in the analysis of the three-dimensional 
                        stability of a two-dimensional flow field. The formulation (x-yz) has a z-velocity component that is independent of the z-direction.
======================= ==================================================

------------
**Examples**
------------

The following is a sample card that sets the coordinate system to Cartesian:
::

   Coordinate System = CARTESIAN

-------------------------
**Technical Discussion**
-------------------------

Note the coordinate ordering for the **CYLINDRICAL** and **SWIRLING** options where the z-direction is first followed by the r-component (which in lay terms means the modeled region/part will appear to be ”lying down.”) If the **SWIRLING** option is activated, *Goma* expects a third momentum equation for the :math:`\theta`-direction, i.e. *EQ = momentum3*, as explained in the equation section. The third component is basically the azimuthal :math:`\theta`-velocity component, and the appropriate boundary conditions must be applied, e.g., on the w-component as described in the *Category 4* boundary conditions for *Fluid Momentum Equations*.



--------------
**References**
--------------

No References.
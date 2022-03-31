*****************************
Electromagnetic Incident Wave
*****************************

::

   Electromagnetic Incident Wave = {model_name} {float} []

-----------------------
**Description / Usage**
-----------------------

This is to specify the incident wave for boundary conditions and post processing.
This should probably be the same for all materials.


PLANE_Z_WAVE
      * <float1> magnitude :math:`E_0` Plane Z wave in the x direction:
        :math:`(E_inc)_z = E_0 exp(i\omega x)`


------------
**Examples**
------------

Following are sample cards:

::

   Electromagnetic Incident Wave = PLANE_Z_WAVE 1.

-------------------------
**Technical Discussion**
-------------------------

This is used to set an incident wave for time-harmonic Maxwell problems, currently done through the absorbing boundary condition, and used in post processing the scattered and incident wave.


--------------
**References**
--------------

No References.

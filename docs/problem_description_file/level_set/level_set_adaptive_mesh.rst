***************************
Level Set Adaptive Mesh
***************************

::

	Level Set Adaptive Mesh = {yes|no}

-----------------------
Description / Usage
-----------------------

This card specifies whether to enable level set adaptive mesh through :code:`Omega_h` library

yes|on
    Adaptivity is enabled

no|off
    Default, no adaptive mesh

------------
Examples
------------

This is a sample card:
::

    Level Set Adaptive Mesh = on

-------------------------
Technical Discussion
-------------------------

This requires first order trianglular or tetrahedral meshes

--------------
References
--------------

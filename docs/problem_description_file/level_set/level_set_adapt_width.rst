***************************
Level Set Adapt Width
***************************

::

	Level Set Adapt Width = <float>

-----------------------
Description / Usage
-----------------------

This card specifies the width to apply an "Inner" adaptive size using :code:`Omega_h` library

<float>
    Width around the level set interface to apply the inner size


------------
Examples
------------

This is a sample card:
::

    Level Set Adapt Width = 0.5
    Level Set Adapt Inner Size = 0.1
    Level Set Adapt Outer Size = 0.8

-------------------------
Technical Discussion
-------------------------

This requires first order trianglular or tetrahedral meshes

--------------
References
--------------

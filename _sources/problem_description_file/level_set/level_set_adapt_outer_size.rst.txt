***************************
Level Set Adapt Outer Size
***************************

::

	Level Set Adapt Outer Size = <float>

-----------------------
Description / Usage
-----------------------

This card specifies the size outside of the "Adapt Width" range near the interface, using :code:`Omega_h` library

<float>
    ISO element size outside of the Adapt Width region


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

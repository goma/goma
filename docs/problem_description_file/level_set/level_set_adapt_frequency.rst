***************************
Level Set Adapt Frequency
***************************

::

	Level Set Adapt Frequency = <integer>

-----------------------
Description / Usage
-----------------------

This card specifies how often to adapt the mesh

<integer>
    Will adapt the mesh after <integer> time steps


------------
Examples
------------

This is a sample card:
::

    Level Set Adapt Width = 0.5
    Level Set Adapt Inner Size = 0.1
    Level Set Adapt Outer Size = 0.8
    Level Set Adapt Frequency = 10

-------------------------
Technical Discussion
-------------------------

This requires first order trianglular or tetrahedral meshes

--------------
References
--------------

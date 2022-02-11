************************
Pseudo-Solid Lame LAMBDA
************************

::

   Pseudo-Solid Lame LAMBDA = CONSTANT <float> [M/Lt2]

-----------------------
**Description / Usage**
-----------------------

This card is required only for *TOTAL_ALE* mesh motion types (see *Mesh Motion* card)
and is used to specify the model for the Lame coefficient :math:`\lambda` for the mesh motion
elasticity (see Sackinger et al., 1995).

This material parameter currently has only one possible model type (CONSTANT)
with only a single required input value, as follows:

+-----------------+------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the *Lame LAMBDA* coefficient model.                                              |
|                 | * <float1> - Standard value of μ (or the shear modulus G for the mesh). See              |
|                 |   *Pseudo-Solid Constitutive Equation* card.                                             |
+-----------------+------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Pseudo-Solid Lame LAMBDA = CONSTANT 1.

This card specifies that the current material have a constant shear modulus of 0.5 for
the mesh elasticity. Note that the real-solid mesh Lame MU is set with the *Lame MU*
card.

-------------------------
**Technical Discussion**
-------------------------

See discussion on *Lame LAMBDA* card and *Solid Constitutive Equation* card for more
details. The main difference here is that this modulus is applied only to the moving
mesh, and not the real solid as in an ALE solid mechanics simulation.



--------------
**References**
--------------

GT-005.3: THE NEW TOTAL-ARBITRARY-LAGRANGIAN-EULERIAN (TALE)
CAPABILITY and its applicability to coating with/on deformable media, August 6,
1999, P. R. Schunk

Sackinger, P. A., Schunk, P. R. and Rao, R. R. 1995. "A Newton-Raphson Pseudo-Solid
Domain Mapping Technique for Free and Moving Boundary Problems: A Finite
Element Implementation", J. Comp. Phys., 125 (1996) 83-103.

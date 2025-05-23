***************************
Usage Notes on the GD Cards
***************************

Following are several examples of uses of the Generalized Dirichlet conditions:

* For a circular boundary ( with radius 1, center at (0,0), :math:`x^2 + y^2 = 1` ):

::

	BC = GD_PARAB SS   1   R_MESH2 0 MESH_POSITION2     0   -1.   0.   1.
	BC = GD_PARAB SS   1   R_MESH2 0 MESH_POSITION1     0   -0.   0.   1.

* For a planar boundary ( :math:`2x + y = 1` )

::

	BC = GD_LINEAR SS   1   R_MESH1 0 MESH_POSITION1     0   -1.   2.
	BC = GD_LINEAR SS   1   R_MESH1 0 MESH_POSITION2     0   0.   1.

* For a parabolic inflow velocity profile ( :math:`u = 1 – 2y – 3y^2` ):

::

	BC = GD_LINEAR SS   4   R_MOMENTUM1 0 VELOCITY1     0   0.   -1.
	BC = GD_PARAB SS   4   R_MOMENTUM1 0 MESH_POSITION2     0   1.   -2.   -3.

* For a distinguishing condition where the mesh is an iso-concentration surface (:math:`C = 0.2` with 
  mesh equations rotated):

::

	BC = GD_CONST SS   2   R_MESH_NORMAL 0 MASS_FRACTION     0   0.  2

* For a temperature boundary condition with APREPRO constants :math:`C_i` of the form

.. math::

   T = C_1 + C_2x + C_3x^2 + C_4x^3 + C_5x^4 + C_6x^5 + C_7x^6 

::

	BC = GD_LINEAR SS   2   R_ENERGY 0 TEMPERATURE     0   -1
	BC = GD_POLYN SS   2   R_ENERGY 0 MESH_POSITION1     0   {c1 c2 c3 c4 c5 c6 c7}

Note, in the first three examples, two cards are combined to create a single boundary
condition that is a function of two variables. Thus, with a little creativity, the Generalized
Dirichlet conditions can replace many of the other boundary condition types.

To help generalize the Dirichlet conditions even more, GD_TIME can be used to modulate
any combination of spatial GD conditions (the CONST, LINEAR, PARAB,
POLYN, CIRC and TABLE options above) which appears prior to the set. Some examples
here are warranted:

* For a parabolic inflow velocity profile which is ramped from zero to a linearly growing multiplier times ( :math:`u = 1 – 2y – 3y^2` ):

::

	BC = GD_PARAB SS   4   R_MOMENTUM1 0 MESH_POSITION2     0   1.   -2.   -3.
	BC = GD_TIME SS   4   R_MOMENTUM1   0   LINEAR     0   0.   1.
	BC = GD_LINEAR SS   4   R_MOMENTUM1   0   VELOCITY1      0   0.   -1.

(This set of 3 conditions actually applies :math:`f(x, y, z, t, u) = 1t \left(1 – 2y – 3y^2 \right) – u = 0` in place of the x-momentum equation. )

* For a sinusoidally time-varying roller surface with equation :math:`(x – x_0)^2 + (y – y_0)^2 = {R_0}^2` with a frequency of 2. and a phase lag of 10:

::

	BC = GD_PARAB SS   1   R_MESH_NORMAL   0   MESH_POSITION2 0 {x0*x0 + y0*y0 - R0*R0} {-2.*y0} 1
	BC = GD_PARAB SS   1   R_MESH_NORMAL   0   MESH_POSITION1     0 {0.}   {-2.*x0} 1
	BC = GD_TIME SS   1   R_MESH_NORMAL   0   SINUSOIDAL     0   10.   2.

This set of cards applies :math:`f(x, y, z, t) (x – x_0)^2 + (y – y_0)^2 - \sin(2t + 10)
{R_0}^2 = 0` to the normal component of the mesh equations along side set 1.

* For a sinusoidally varying gap on a slot coater, the substrate has been made to oscillate according to :math:`f(x, y, t) = y – 3 \sin(t/4 + 5) = 0` :

::

	BC = GD_LINEAR SS   9   R_MESH2   0   MESH_POSITION1     0   -3.  0   0.
	BC = GD_TIME SS   9   R_MESH2   0   SINUSOIDAL     0   5.   0.25
	BC = GD_LINEAR SS   9   R_MESH2   0   MESH_POSITION2     0   0.   1.0

* Setting the u-velocity on an inlet boundary for a power law fluid:

::

        BC = GD_LINEAR SS   1   R_MOMENTUM1   0   VELOCITY1     0   0.   -1.
        BC = GD_TABLE SS   1   R_MOMENTUM1   0   MESH_POSITION2 0 1.0 LINEAR

        $ r/R         Uz
        0.000000      1.666667
        0.050000      1.666458
        0.100000      1.665000
        0.150000      1.661042
        0.200000      1.653333
        0.250000      1.640625
        0.300000      1.621667
        .             .
        .             .
        .             .
        0.900000      0.451667
        0.950000      0.237708
        1.000000      0.000000
        END TABLE

* Setting the inlet concentration profile for species 0 from data in y0.table

::

	BC = GD_LINEAR SS   1   R_MASS   0   MASS_FRACTION     0   0.0   -1.0
	BC = GD_TABLE SS   1   R_MASS   0   MESH_POSITION2 0 1.0 LINEAR FILE = y0.table

* Setting the inlet concentration profile for species 0 from an implicit relation.

Occasionally, we have analytic representations that are in the wrong form. For example, in particulate
suspension modelling, a relation exists that gives the radial coordinate as a function of the concentration,
i.e. :math:`r = F(C)`, where :math:`F` is a non-linear relation. We would prefer it the other way around.
We can use GD_TABLE to solve this dilemma. First, a file is prepared with the two columns,
eqn.table for example:

==== =======
C_0  F(C_0)
C_1  F(C_1)
.    .
.    . 
C_N  F(C_N)
==== =======

This just requires function evaluation. In the input deck, we then use the following cards

::

	BC = GD_LINEAR SS   1   R_MASS   0   MESH_POSITION2     0   0.0   -1.0
	BC = GD_TABLE SS   1   R_MASS   0   MASS_FRACTION 0 1.0 LINEAR FILE = eqn.table

and the right inlet concentration profile results.


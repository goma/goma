*************
**LS_QLASER**
*************

::

	BC = LS_QLASER LS <integer> <float1> <float2> <float3> <float4>

-----------------------
**Description / Usage**
-----------------------

**(EMB/ENERGY)**

This boundary condition card specifies heat flux model derived from a laser welding
application. This heat flux value is applied as an “embedded” source term on the heat
conservation equation at the zero level set contour (cf. *BC = Q_LASER_WELD* for
ALE surfaces). It can be used both when subgrid or subelement integration is being
used. The <float_list> has twenty-seven parameters; definitions of the input parameters
are as follows:

A description of the input parameters follows:

============= ==================================================================
**LS_QLASER** Name of the boundary condition.
**LS**        This string is used to indicated that this is a “boundary”
              condition is applied at an internal phase boundary defined
              by the zero contour of the level set function.
<integer>     An integer parameter than is permitted to take one of three
              values -1, 0, or 1. Depending upon the choice of this
              parameter the heat flux value is applied to the negative
              phase, both phase, or the positive phase, respectively.
              Details are given below.
<float 1>     Nominal power of laser.
<float 2>     Power of laser at base state (simmer).
<float 3>     Base value of surface absorptivity.
<float 4>     Switch to allow tracking of normal component of liquid
              surface relative to laser beam axis for surface absorption
              (0 = OFF, 1 = ON)
<float 5>     Cutoff time for laser power.
<float 6>     Time at which laser power drops to 1/e.
<float 7>     For pulse weld, the laser power overshoot (%) of peak
              power at time to reach peak laser power.
<float 8>     Radius of laser beam.
<float 9>     For pulse weld, the time for laser pulse to reach peak power.
<float 10>    For pulse weld, the time for laser pulse to reach steady
              state in power.
<float 11>    Switch to either activate laser power distribution from
              beam center based on absolute distance (0) or based on
              radial distance in 2D plane (1).
<float 12>    Location of laser beam center (x-coordinate).
<float 13>    Location of laser beam center (y-coordinate).
<float 14>    Location of laser beam center (z-coordinate).
<float 15>    Laser beam orientation, normal to x-coordinate of body.
<float 16>    Laser beam orientation, normal to y-coordinate of body.
<float 17>    Laser beam orientation, normal to z-coordinate of body.
<float 18>    For pulse weld, spot frequency.
<float 19>    For pulse weld, total number of spots to simulate.
<float 20>    Switch to set type of weld simulation. (0=pulse weld,
              1=linear continuous weld, -1=pseudo pulse weld,
              2=sinusoidal continous weld)
<float 21>    For pulse weld, spacing of spots.
<float 22>    For radial traverse continuous weld, radius of beam travel.
<float 23>    Switch to activate beam shadowing for lap weld
              (0=OFF, 1=ON). Currently only active for ALE simulations.
<float 24>    Not active, should be set to zero.
<float 25>    For continuous weld, laser beam travel speed in xdirection
              (u velocity).
<float 26>    For continuous weld, laser beam travel speed in ydirection
              (v velocity).
<float 27>    For continuous weld, laser beam travel speed in zdirection
              (w velocity).
============= ==================================================================

------------
**Examples**
------------

An example:
::

   BC = LS_QLASER LS -1 4.774648293 0 0.4 1 1 1.01 4.774648293 0.2
   0.01 0.01 1 0.005 0 -0.198 -1 0 0 0.025 1 1 0.2032 -1000 0 0 0 0
   0.0254

-------------------------
**Technical Discussion**
-------------------------

This is the level-set counterpart to *BC = Q_LASER* which is the same boundary
condition applied to a parameterized mesh surface. Please see the discussion of that
input record for the functional form of this boundary condition.



--------------
**References**
--------------

No References. 

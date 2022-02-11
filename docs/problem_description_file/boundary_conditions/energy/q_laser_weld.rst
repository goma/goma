****************
**Q_LASER_WELD**
****************

::

	BC = Q_LASER_WELD SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition card specifies the thermal boundary conditions for laser
welding. The <float_list> requires twenty-seven values be specified; definitions of the
input parameters are as follows:

================ ===================================================================
**Q_LASER_WELD** Name of the boundary condition (<bc_name>).
**SS**           Type of boundary condition (<bc_type>), where **SS**
                 denotes side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain.
<float1>         Nominal power of laser.
<float2>         Power of laser at base state (simmer).
<float3>         Base value of surface absorptivity.
<float4>         Switch to allow tracking of normal component of liquid
                 surface relative to laser beam axis for surface absorption
                 (0 = OFF, 1 = ON)
<float5>         Cutoff time for laser power.
<float6>         Time at which laser power drops to 1/e.
<float7>         For pulse weld, the laser power overshoot (%) of peak
                 power at time to reach peak laser power.
<float8>         Radius of laser beam.
<float9>         For pulse weld, the time for laser pulse to reach peak power.
<float10>        For pulse weld, the time for laser pulse to reach steady
                 state in power.
<float11>        Switch to either activate laser power distribution from
                 beam center based on absolute distance (0) or based on
                 radial distance in 2D plane (1).
<float12>        Location of laser beam center (x-coordinate).
<float13>        Location of laser beam center (y-coordinate).
<float14>        Location of laser beam center (z-coordinate).
<float15>        Laser beam orientation, normal to x-coordinate of body.
<float16>        Laser beam orientation, normal to y-coordinate of body.
<float17>        Laser beam orientation, normal to z-coordinate of body.
<float18>        For pulse weld, spot frequency.
<float19>        For pulse weld, total number of spots to simulate.
<float20>        Switch to set type of weld simulation. (0=pulse weld,
                 1=linear continuous weld, -1=pseudo pulse weld,
                 2=sinusoidal continous weld)
<float21>        For pulse weld, spacing of spots.
<float22>        For radial traverse continuous weld, radius of beam travel.
<float23>        Switch to activate beam shadowing for lap weld
                 (0=OFF, 1=ON). Currently only active for ALE simulations.
<float24>        Not active, should be set to zero.
<float25>        For continuous weld, laser beam travel speed in xdirection
                 (u velocity).
<float26>        For continuous weld, laser beam travel speed in ydirection
                 (v velocity).
<float27>        For continuous weld, laser beam travel speed in zdirection
                 (w velocity).
================ ===================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = Q_LASER_WELD SS 10 4.774648293 0 0.4 1 1 1.01 4.774648293 0.2
   0.01 0.01 1 0.005 0 -0.198 -1 0 0 0.025 1 1 0.2032 -1000 0 0 0 0 0.0254

-------------------------
**Technical Discussion**
-------------------------

Several required pieces of information to use this boundary condition are not in final
form, and the user can expect future changes and improvements. Below is a listing of
some of these parameters:

* This boundary condition requires that node sets 1001 is defined in the EXODUS II
  file. NS 1001 should include the point at the center of the keyhole on the surface
  closest to the beam.

* Currently the laser flux distribution is set as a fixed exponential distribution.
  Plans are to include more options including a user-defined exponential and a TABLE
  option.

* Correlations are used to specify the evaporation energy loss. Currently only iron
  and ice correlations exist; the appropriate correlation is selected based on the 
  value set for the *Solidus Temperature* (in *Thermal Properties* portion of the 
  material file).




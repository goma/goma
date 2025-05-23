****************
**VELO_SLIP_LS**
****************

::

	BC = VELO_SLIP_LS SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is applied only in problems involving embedded interface
tracking, that is, level set or volume of fluid. The boundary condition serves two major
purposes: first to allow for slip in the vicinity of a moving contact line and second to
facilitate impact of a dense fluid on a substrate, displacing a less dense fluid (e.g. water
drop and air displacement). . Elsewhere, this boundary condition enforces a no-slip
condition between fluid and substrate. A more detailed description is given below.

This boundary condition is most often used in conjunction with the FILL_CA,
WETTING_SPEED_LINEAR, or WETTING_SPEED_BLAKE boundary conditions.
These apply forces to contact lines in order to simulate wetting line motion. These
forces are applied in a weak sense to the same regions near the interface so it is
necessary to use VELO_SLIP_LS with a large slipping coefficient so that effectively
no-slip is relaxed completely near the interface.

Definitions of the input parameters are as follows:

================ ======================================================
**VELO_SLIP_LS** Name of the boundary condition.
**SS**           Type of boundary condition (<bc_type>), where **SS**
                 denotes side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain.
<float1>         *alpha*, or *slip_width*, a characteristic length scale around
                 the contact line that will be used to apply Navier Slip
                 Condition with :math:`\beta` 0 coefficient. This length scale
                 is also used to detect the thickness of light-phase (gas) between
                 the substrate denoted by the sideset, and the zero-level-set 
                 contour (or the boundary between liquid and gas). If this 
                 distance is less than 8*slip_width, then perfect slip in the
                 gas phase is allowed to help facilitate contact. See
                 discussion below.
<float2>         :math:`\beta` 0, the slip coefficient near the contact line. 
                 The inverse of :math:`\beta` 0 defines the scaling between stress and
                 slip. The parameter supplied on the input deck is used only within
                 a lengthscale *slip_width* setting around the contact line..
                 Elsewhere, the slip coefficient is uniformly set to :math:`\beta` 1.
                 Hence, this parameter is usually set to a large value to
                 allow for perfect slip.
<float3>         :math:`v_{s,x}`, the x-component of surface velocity vector. This
                 would be the x-component of the fluid velocity if a noslip
                 condition were applied.
<float4>         :math:`v{s,y}`, the y-component of surface velocity vector. This
                 would be the y-component of the fluid velocity if a noslip
                 condition were applied.
<float5>         :math:`v_{s,z}`, the z-component of surface velocity vector. This
                 would be the z-component of the fluid velocity if a noslip
                 condition were applied.
<float6>         :math:`\beta` 1, the slip coefficient away from the contact line. The
                 inverse of :math:`\beta` 1 defines the scaling between stress and slip.
                 Hence, this parameter is usually set to a small value
                 (like 1e-6) to allow for no-slip.
================ ======================================================

------------
**Examples**
------------

Following is a sample card without the optional parameters:
::

   BC = VELO_SLIP_LS SS 10 0.05 100000.0 0.0 0.0 1.e-6

The large value of slip coefficient ensures nearly perfect slip in the region around the
interface, a region that has a half-width of 0.05 centered about the contact line. Away
from the contact line (outside the hat function of width 0.05), the slip coefficient is 1.e-6, which corresponds to significantly less slip. Note also that if the substrate defined by SS 10 is in contact with gas (light phase), and a liquid front (zero-level set) is nearby, or within a distance 8 times 0.05, then the light phase is allowed to slip along the wall with a coefficient of 100000.0. This is to help facilitate contact.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition was originally developed to allow for fluid slip near a
dynamic contact line, a necessary condition for dynamic wetting line motion when the
contact angle is not 180 degrees (viz. rolling motion condition). The slippage
mechanism was deployed through the use of Navier’s slip condition, which basically
goes as

.. figure:: /figures/243_goma_physics.png
	:align: center
	:width: 90%

Here :math:`\beta` is the slip coefficient, which is taken to be variable depending on its proximity
to the contact line (through the “slip_width” parameter). Note that the smaller the :math:`\beta`, the
more no-slip is enforced. The left hand side of this condition is the fluid traction on the
substrate. :math:`\underline{v}_s` is the velocity of the substrate, specified component-wise with {vx} {vy}
{vz}. This base functionality of applying the Navier slip condition still exists in this
condition, but in addition it was furbished to allow for complete slip on the boundary if
a gas film is being displaced by liquid. In this latter case, complete slip is a mechanism
(subgrid event) that allows for the otherwise infinite stress to be relieved so that the
liquid can make contact with the solid. The perfect slip condition at the substrate/gas
surface is activated by just setting the slip coefficient to the large value, as this
condition does anyway in the vicinity of a contact line. The “gas phase” is determined
by determining which phase is the lighter one based on the density specification. The
figure below details more on how this condition works for wetting/dewetting and for
incipient liquid/solid impact.

Some more usage notes as follows:

* The slip coefficient function is computed as 
  :math:`\beta` = :math:`\beta` 0 :math:`\delta` ( :math:`\phi` ) + :math:`\beta_\infty`,
  where the delta
  function is a level-set hat function centered around the zero level set contour where
  it intersects the boundary. It has a length sacale associated with it which is called
  “alpha”, and that basically sets the length over which the :math:`\beta` 0 is applied as the slip
  parameter and it is large, leading to a shear-stress-free or slippery region in both
  the gas and liquid phases. :math:`\beta_{INF}` is taken as real small (typically 1.e-6 or less) and is
  applied away from the contact line, and hence forces a true “no-slip” condition.

* The most recent addition to this condition is the functionality that adds perfect slip
  to a wall in the gas phase as it is displaced during near contact state by a liquid
  phase. This capability is of course applicable only to level-set capillary
  hydrodynamics problems. Level-set methods have been plagued by the fact that it
  is hard to break down the displaced phase (e.g. gas phase) as a liquid phase surface
  flows towards a solid boundary. Theoretically this event requires an infinite stress,
  in the continuum. To relieve this stress and promote a collapse and wetting, we add
  perfect slip in the gas phase at near contact conditions, which reduces the
  lubrication pressure in the gas film and promotes breakdown. This of course
  introduces more length scales. First, the length scale over which slip is applied
  (this is the alpha parameter described above) and seond is the length scale over
  which “nearness” of the liquid phase to the substrate is considered to be “close
  enough” to allow for perfect slip. Right now this “nearness” length scale is
  arbitrarily set to 8*alpha. A third length scale is that which we use to declare
  contact. We currently have that set to 1.e-6*alpha. After contact is declared,
  VELO_SLIP_LS reverts to the form under the first bullet. The figure below
  hopefully clarifies the condition a little better.

.. figure:: /figures/244_goma_physics.png
	:align: center
	:width: 90%




.. TODO -Lines 92 and 141 have pictures that need to be swapped with the correct equations.

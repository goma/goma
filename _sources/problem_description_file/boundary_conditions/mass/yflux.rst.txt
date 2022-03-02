*********
**YFLUX**
*********

::

	BC = YFLUX SS <bc_id> <integer1> <float1> <float2>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary condition card is used to specify the mass flux of a given species normal
to the boundary (or interface) using a mass transfer coefficient. When used in
conjunction with the *KIN_LEAK* card, the *YFLUX* card also enables the determination
of velocity normal to the moving boundary at which the *YFLUX* boundary condition is
applied.

Definitions of the input parameters are as follows:

========== ======================================================================
**YFLUX**  Name of the boundary condition (<bc_name>).
**SS**     Type of boundary condition (<bc_type>), where **SS** denotes
           side set in the EXODUS II database.
<bc_id>    The boundary flag identifier, an integer associated with
           <bc_type> that identifies the boundary location (side set in
           EXODUS II) in the problem domain.
<integer>  *i*, species number of concentration.
<float1>   :math:`k_i`, value of mass transfer coefficient of species *i*.
<float2>   :math:`c_i^\infty` value of reference concentration of species *i*.
========== ======================================================================

------------
**Examples**
------------

Following are two sample cards:
::

   BC = YFLUX SS   3 0   0.12   0.

::

   BC = YFLUX SS   3 1   0.05   0.

-------------------------
**Technical Discussion**
-------------------------

Specifically, the species mass flux is given by

.. figure:: /figures/136_goma_physics.png
	:align: center
	:width: 90%

where *n* is the unit vector normal to the boundary, :math:`J_i` is mass flux of species *i*, *v* is the
fluid velocity, :math:`v_m` is the mesh displacement velocity, :math:`k_i` is mass transfer coefficient of
species *i*, :math:`c_i` is concentration of species *i* at the boundary surface,
:math:`c_i^\infty` and is reference
concentration of species *i*. The units of :math:`J_i`, :math:`k_i`, :math:`c_i` and 
:math:`c_i^\infty` depend on the user’s choice. For
example, if :math:`c_i` and :math:`c_i^\infty` are chosen to have units of moles/ :math:`cm^3`, then :math:`k_i` has the unit of cm/
s, and :math:`J_i` has the units of moles/ :math:`cm^2` /s.

For the *KIN_LEAK* and *VNORM_LEAK* cards, the information from *YFLUX* boundary
conditions corresponding to each species is needed. *Goma* automatically searches for
these boundary conditions and uses an extra variable in the BC data storage to record
the boundary condition number of the next *YFLUX* condition in a linked list; when the
extra storage value is -1, there are no more *YFLUX* conditions at this boundary.


--------
**FAQs**
--------

A question was raised regarding the use of volume flux in Goma; the following portion
of the question and response elucidate this topic and the subject of units. Note the
references in the response are to the Version 2.0 Goma User’s Manual.

**Question**: ... I know what you are calling volume flux is mass flux divided by
density. The point I am trying to make is that the conservation equations in the books I
am familiar with talk about mass, energy, momentum, and heat fluxes. Why do you not
write your conservation equations in their naturally occurring form? If density just so
happens to be common in all of the terms, then it will be obvious to the user that the
problem does not depend on density. You get the same answer no matter whether you
input rho=1.0 or rho=6.9834, provided of course this does not impact iterative
convergence. This way, you write fluxes in terms of gradients with the transport
properties (viscosity, thermal conductivity, diffusion coefficient, etc.) being in familiar
units.

**Answer**: First let me state the only error in the manual that exists with regard to
the convection-diffusion equation is the following:

   :math:`J_t` in the nomenclature table ... should be described as a volume flux with units
   of L/t, i.e., :math:`D ⋅ \Delta_yi` , where D is in :math:`L^2` ⁄ t units.

Now, ... this is actually stated correctly, as it states the :math:`J_i` is a diffusion flux (without
being specific); to be more specific here, we should say it is a "volume flux of species
*i*." So, in this case D is in L ⋅ L ⁄ t units, :math:`y_i` is dimensionless and it is immaterial that (the
mass conservation equation) is multiplied by density or not, *as long as density is
constant*.

Now, in *Goma* we actually code it up EXACTLY as in the ... (mass conservation
equation), i.e., there are no densities anywhere for the *FICKIAN* diffusion model. For
the *HYDRO* diffusion model, we actually compute a :math:`J_i` ⁄ ρ in the code, and handle variable density changes through that p . In that case :math:`J_i` as computed in Goma is a mass
flux vector, not a volume flux vector, but by dividing it by p and sending it back up to
the mass conservation equation it changes back into a volume flux. i. e., everything is the same.

Concerning the units of the mass transfer coefficient on the *YFLUX* boundary
condition, the above discussion now sets those. *Goma* clearly needs the flux in the
following form:

.. figure:: /figures/137_goma_physics.png
	:align: center
	:width: 90%

and dimensionally for the left hand side

.. figure:: /figures/138_goma_physics.png
	:align: center
	:width: 90%

where D is in units :math:`L^2` /t , the gradient operator has units of 1/L so K HAS to be in units
of L/t (period!) because :math:`y_i` is a fraction.

So, if you want a formulation as follows:

.. figure:: /figures/139_goma_physics.png
	:align: center
	:width: 90%

then K’s units will have to accommodate for the relationship between :math:`p_i` and :math:`y_i` in the
liquid, hopefully a linear one as in Raoult’s law, i.e. if :math:`p_i` = :math:`P_V`:math:`y_i` where :math:`P_v` is the vapor
pressure, then

.. figure:: /figures/140_goma_physics.png
	:align: center
	:width: 90%

and so K on the *YFLUX* command has to be :math:`KP_v` ....and so on.

Finally, you will note, since we do not multiply through by density, you will have to
take care of that, i. e., in the Price paper he gives K in units of t/L. So, that must be converted as follows:

.. figure:: /figures/141_goma_physics.png
	:align: center
	:width: 90%

This checks out!

--------------
**References**
--------------

Price, P. E., Jr., S. Wang, I. H. Romdhane, 1997. “Extracting Effective Diffusion
Parameters from Drying Experiments”, AIChE Journal, 43, 8, 1925-1934.

.. TODO - Lines 54, 120, 126, 135, 143, and 152 all contain picture that need to be changed into the correct equations. 
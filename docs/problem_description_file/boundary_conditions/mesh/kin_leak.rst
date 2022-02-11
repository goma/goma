************
**KIN_LEAK**
************

::

	BC = KIN_LEAK SS <bc_id> <float1> <float2>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MESH)**

This boundary condition card is used as a distinguishing condition - kinematic with
mass transfer on mesh equations. The flux quantity is specified on a per mass basis so
heat and mass transfer coefficients are in units of L/t.

Definitions of the input parameters are as follows:

==================== ===============================================================
**KIN_LEAK**         Name of the boundary condition (<bc_name>).
**SS**               Type of boundary condition (<bc_type>), where **SS** denotes
                     side set in the EXODUS II database.
<bc_id>              The boundary flag identifier, an integer associated with
                     <bc_type> that identifies the boundary location (side set in
                     EXODUS II) in the problem domain.
<float1>             Mass transfer coefficient for bulk fluid (species *n* +1).
<float2>             Driving force concentration in external phase.
==================== ===============================================================

Please see Technical Discussion regarding the appropriate units for the mass transfer
coefficient and concentration in the external phase. For a pure liquid case, these inputs
are read directly from this card, while for a multi-component case these values are read
from *YFLUX* boundary conditions corresponding to each species that is needed. See
following examples.

------------
**Examples**
------------

Following are two sample input cards:

Pure Liquid Case
::

     BC = KIN_LEAK SS 3   0.1   0.

Two Component Case
::

     BC = KIN_LEAK SS 3   0.   0.

::

     BC = YFLUX SS 3   0   0.12   0.

Note, in the two component case, when *Goma* finds the *KIN_LEAK* card, it scans the
input deck to locate the applicable *YFLUX* conditions associated with side set 3 and
creates a linked list which is used by the applying function (*kin_bc_leak*). The
existence of this list is denoted in *Goma* by the addition of an integer into an unused
field of the BC structure for side set 3. The bulk fluid constitutes the second component
and is non-volatile so it requires no *YFLUX* card; a second volatile species would
require a second *YFLUX* input card.

-------------------------
**Technical Discussion**
-------------------------

Functionally, the *KIN_LEAK* boundary condition can be represented as the following:

.. figure:: /figures/045_goma_physics.png
	:align: center
	:width: 90%

where **EQUATION** is the vector velocity; **EQUATION** is the velocity of the boundary itself (not independent
from the mesh velocity); **EQUATION** is the normal vector to the surface; **EQUATION** is the concentration
of species *i*; **EQUATION** is the ambient concentration of species *i* at a distance from the surface
of interest and **EQUATION** is the mass transfer coefficient for species *i*. This function returns a
volume flux term to the equation assembly function.

*KIN_LEAK* is implemented through function *kin_bc_leak*; it sums the fluxes for all
species plus the bulk phase evaporation. These fluxes are computed via several other
function calls depending on the particular flux condition imposed on the boundary.
(See various *YFLUX* * cards for Mass Equations.) However, at the end of the
*kin_bc_leak* function, the accumulated flux value is assigned to variable *vnormal*, i.e.,
the velocity of fluid relative to the mesh. The apparent absence of a density factor here
to convert a volume flux to a mass flux is the crucial element in the proper usage of the
flux boundary conditions. The explanation is rooted in the formulation of the
convective-diffusion equation.

The convective-diffusion equation in *Goma* is given as

.. figure:: /figures/046_goma_physics.png
	:align: center
	:width: 90%

with mass being entirely left out of the expression. *J* is divided by density before
adding into the balance equation; this presumes that volume fraction and mass fraction
are equivalent. The users must be aware of this. This formulation is certainly
inconvenient for problems where volume fraction and mass fraction are not equal and
multicomponent molar fluxes are active elements of an analysis. However, *kin_bc_leak*
is entirely consistent with the convective-diffusion equation as a velocity is a volume
flux, and multiplied by a density gives a proper mass flux. If y\ :sub:`i` is a mass 
concentration,
and h\ :sub:`i` were in its typical velocity units, the result is a mass flux; if 
y\ :sub:`i` is a volume
fraction, then we have a volume flux. So *kin_bc_leak* is consistent.

The burden here lies with the user to be consistent with a chosen set of units. A
common approach is to build density into the mass transfer coefficient h\ :sub:`i` .


--------
**FAQs**
--------

1. See the FAQ pertaining to “Continuation Strategies for Free Surface Flows” on the
   *DISTNG* boundary condition card.

2. A question was raised regarding the use of volume flux in *Goma*; the following
   portion of the question and response elucidate this topic and the subject of units. Being from several emails exchanged during January 1998, the deficiencies or lack of clarity have since been remedied prior to *Goma* 4.0, but the discussions are relevant for each user of the code.

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

**Answer**: ... First let me state the only error in the manual that exists with regard to
the convection-diffusion equation (CDE) is the following:

	J\ :sub:`i` in the nomenclature table should be described as a volume flux with units 
	of L/t, i.e., D ⋅ ∇y\ :sub:`i`, where D is in L\ :sup:`2`/t units.

Now, this is actually stated correctly elsewhere, as it states the J\ :sub:`i` is a 
diffusion flux
(without being specific); to be more specific here, we should say it is a "volume flux of
species i." So, in this case D is in L ⋅ L/t units y\ :sub:`i`, is dimensionless and it 
is immaterial
that the CDE is multiplied by density or not, *as long as density is constant*.

Now, in *Goma* we actually code it with no densities anywhere for the FICKIAN diffusion model. For the HYDRO diffusion model, we actually compute a J\ :sub:`i`  /ρ in 
the code, and handle variable density changes through that ρ. In that case J\ :sub:`i` 
as computed in *Goma* is a mass flux vector, not a volume flux vector, but by dividing it 
by ρ and sending it back up to the CDE it changes back into a volume flux. i. e., 
everything is the same.

Concerning the units of the mass transfer coefficient on the YFLUX boundary
condition, the above discussion now sets those. *Goma* clearly needs the flux in the
following form:

.. figure:: /figures/047_goma_physics.png
	:align: center
	:width: 90%

and dimensionally for the left hand side

.. figure:: /figures/048_goma_physics.png
	:align: center
	:width: 90%

where D is in units L\ :sup:`2`/t, the gradient operator has units of 1/L so K *has* to 
be in units
of L/t (period!) because y\ :sub:`i` is a fraction.

So, if you want a formulation as follows:

.. figure:: /figures/049_goma_physics.png
	:align: center
	:width: 90%

then K’s units will have to accommodate for the relationship between p\ :sub:`i` and 
y\ :sub:`i` in the
liquid, hopefully a linear one as in Raoult’s law, i.e. if p\ :sub:`i` = PvV\ :sub:`i` 
where Pv is the vapor
pressure, then

.. figure:: /figures/050_goma_physics.png
	:align: center
	:width: 90%

and so K on the YFLUX command has to be KPv....and so on.

Finally, you will note, since we do not multiply through by density, you will have to
take care of that, i. e., in the Price paper he gives K in units of t/L. So, that must be
converted as follows:

.. figure:: /figures/051_goma_physics.png
	:align: center
	:width: 90%

This checks out!

--------------
**References**
--------------

Price, P. E., Jr., S. Wang, I. H. Romdhane, “Extracting Effective Diffusion Parameters
from Drying Experiments,” AIChE Journal, 43, 8, 1925-1934 (1997)

..
	 TODO - The picture on line 72 and 94 need to be exchanged with the equation. In lines 76-80, where it says "**EQUATION**" there is supposed to be something from the equation that needs to be written. 
************
**Porosity**
************

::

   Porosity = {model_name} <float> []

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the porosity, which is required for the
Brinkman or Darcy formulations for flow through porous media, viz. for
**POROUS_BRINKMAN, POROUS_TWO_PHASE, POROUS_SATURATED,**
and **POROUS_UNSATURATED** media types (see *Media Type* card).

Definitions of the {model_name} and <float> parameters are as follows:

+---------------+----------------------------------------------------------------------------+
|**CONSTANT**   |Name {model_name} of the constant porosity model.                           |
|               |                                                                            |
|               | * <float> - Value of porosity.                                             |
+---------------+----------------------------------------------------------------------------+
|**DEFORM**     |Name {model_name} of the model for a porosity that varies with deformation  |
|               |of the porous medium. A conservation balance is required for the solid      |
|               |material skeleton and is invoked in the equation specification section      |
|               |(see *EQ* section).                                                         |
|               |                                                                            |
|               | * <float> - Value of porosity (in the stress-free-state, i.e.,             |
|               |   undeformed state).                                                       |
+---------------+----------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Porosity = DEFORM 0.5

This model will result in a porosity of 0.5 (volume fraction of the interstitial space of a
porous skeleton) in the undeformed or stress-free state, but will allow the porosity to
vary affinely with the volume change invariant of the deformation gradient tensor (see
technical discussion). As mentioned above, the **DEFORM** model requires a field
equation for the mass-conservation of the solid matrix through the porous_deform
equation.

-------------------------
**Technical Discussion**
-------------------------

Porosity is a microstructural attribute of a porous medium which describes the fraction
of volume not occupied by the solid skeleton. For rigid porous media, it is a parameter
that weights the capacitance term (time-derivative term) of the Darcy flow equations
for liquid solvent and gas “solvent” concentrations. It often affects the Saturation
function (see *Saturation* card) and the permeability function (see *Permeability* card).
The references cited below elucidate the role of the porosity parameter in these
equations.

For deformable porous media, *Goma* uses the porosity as a measure of fraction solid
concentration, as a part of a mass balance for the solid skeleton. The reason this
equation is required is a result of the lack of an overall conservation law for the
mixture. Instead, we close the system by individual conservation equations for all
species components in the medium, including the solid; the liquid and gas phase
components are accounted for with individual Darcy flow equations. The conservation
law which governs the porosity assumes there is an affine deformation of the pores
with the overall deformation of the solid, and hence can be written as:

.. figure:: /figures/401_goma_physics.png
	:align: center
	:width: 90%

where F˜ is the deformation gradient tensor, φ0 is the initial porosity, and φ is the
porosity. This equation is invoked with the *porous_deform* option on the *EQ*
specifications.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s capabilities for partially saturated flow in porous media,
September 1, 2002, P. R. Schunk

SAND96-2149: Drying in Deformable Partially-Saturated Porous Media: Sol-Gel
Coatings, Cairncross, R. A., P. R. Schunk, K. S. Chen, S. S. Prakash, J. Samuel, A. J.
Hurd and C. Brinker (September 1996)
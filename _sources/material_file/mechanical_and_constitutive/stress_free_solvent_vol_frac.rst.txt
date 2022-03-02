****************************
Stress Free Solvent Vol Frac
****************************

::

   Stress Free Solvent Vol Frac = CONSTANT <float> []

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for the stress-free solvent volume
fraction, which is the volume fraction of solvents in the solid material in its stress-free
state. This card is used exclusively in materials of *LAGRANGIAN* or *TOTAL_ALE*
Mesh Motion types (see *Mesh Motion* card) which are being modeled as gelled solids
laden with solvent. At the gel-point, the solid is considered to be stress free, after which
a reduction of solvent leads to volume shrinkage and hence a rising stress state.
Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the stress-free solvent volume fraction.                                                                      |
+-----------------+------------------------------------------------------------------------------------------------------------------------------------+
|<float>          |The value of the stress-free solvent volume fraction; this value is unitless.                                                       |
+-----------------+------------------------------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card:

::

   Stress Free Solvent Vol Frac = CONSTANT 0.5

This specification sets the volume fraction of solvent in the material to 50 per cent.
That volume fraction is tantamount to the gel point of the material.

-------------------------
**Technical Discussion**
-------------------------

The stress free state volume fraction of solvent is basically the solvent fraction at which
a material gels, viz., the state at which the material solidifies from a liquid state. This
quantity is used in the continuity equation for incompressible solid materials, through
which is transported by a variety of diffusion models (see *Diffusivity* card). The
continuity equation, viz., *EQ = continuity*, is applied as follows:

.. figure:: /figures/369_goma_physics.png                                                          
   :align: center                                                                                  
   :width: 90%

where the dependent variable is the solid phase pressure (see *Solid Constitutive
Equation* card). Here *det* **F** is the determinant of the deformation gradient tensor, yi is
the volume fraction of component i (specified by the *EQ = species_bulk* card), and y0 is
the volume fraction of total solvents at the stress free state. Clearly, as the solvent
concentration decreases the local volume of solid decreases, creating a rising stress.



--------------
**References**
--------------

GT-001.4: GOMA and SEAMS tutorial for new users, February 18, 2002, P. R. Schunk
and D. A. Labreche

GT-019.1: Elastoviscoplastic (EVP) Constitutive Model in GOMA: Theory, Testing,
and Tutorial, P. R. Schunk, A. Sun, S. Y. Tam (Imation Corp.) and K. S. Chen, January
11, 2001

SAND96-2149: Drying in Deformable Partially-Saturated Porous Media: Sol-Gel
Coatings, Cairncross, R. A., P. R. Schunk, K. S. Chen, S. S. Prakash, J. Samuel, A. J.
Hurd and C. Brinker (September 1996)

.. 
	TODO - Line 50 is a photo that needs to be replaced with the correct equation. 


*******************
Plasticity Equation
*******************


::

   Plasticity Equation = {model_name} [ ]

-----------------------
**Description / Usage**
-----------------------

This optional card specifies the formulation for the potential yielding/plastic flow
regime during solid deformation. This card is not to be used in place of the *Solid
Constitutive Equation* card, but rather supplements that card to describe the constitutive
behavior during plastic deformation. Elastic deformation still proceeds according to the
model specified on the *Solid Constitutive Equation* card (i.e., for regimes that have not
yielded). The single input parameter is defined as

+-------------+---------------------------------------------------------------------------------------+
|{model_name} |Name of the plasticity model. This parameter can have one of the following values:     |
+-------------+---------------------------------------------------------------------------------------+
|**EVP_HYPER**|a constitutive equation that uses the elasticity portion specified on the *Solid       |
|             |Constitutive Equation* card for unyielding material and a complex hyperelastic         |
|             |plasticity equation for the yielding/flowing material as determined by the Von Mises   |
|             |yield criterion.                                                                       |
+-------------+---------------------------------------------------------------------------------------+
|**NO_MODEL** |this, or any value other than EVP_HYPER, will result in no plastic deformation.        |
+-------------+---------------------------------------------------------------------------------------+

Requirements for the use of this model are

* Transient problems only

* *LAGRANGIAN* mesh motion only; no *TALE*

* Continuous media only; no porous media (as specified on the *Media Type* card)

* Elastic Plane Strain models only (i.e., INCOMP_PSTRAIN in the *Solid Constitutive Equation* card)

* a *Plastic Viscosity* card and an *EVP Yield Stress* card must also be supplied.

------------
**Examples**
------------

Following is a sample card:
::

   Plasticity Equation = EVP_HYPER

which specifies hyperelastic elastoviscoplastic model is to be used for a solid phase
constitutive equation. In addition to the Lame coefficients that are still required as the
mechanical properties of the unyielded material, this model also requires a plastic
viscosity and a yield stress, viz.

::

   Plastic Viscosity = LINEAR 1.0 2.0

::

   EVP Yield Stress = CONSTANT 50.0

-------------------------
**Technical Discussion**
-------------------------

Detailed theoretical discussion, usage tutorials and troubleshooting tips for this model
are covered in the EVP tutorial (GT-019.1). Usage examples for four different strain
scenarios are given, including a solid yielding from an applied mechanical load and a
solid yielding from high shrinkage stress during drying.


--------
**FAQs**
--------

*Problem* – Trouble in continuing the first few time steps.


*Solution* – You may have a fast drying case with slow diffusion in the coating. Instead
of decreasing the time step size according to the normal procedure and intuition,
increase the time step size. With fast drying and slow diffusion, the initial
concentration gradient is very steep at the drying surface. This is a very difficult
numerical problem to solve. So when you increase the time step size, in effect, you are
relaxing the concentration gradient the program is solving, that will get you past the
initial numerical difficulty. However, even if the code can handle such a condition, the
concentration and stress profile may appear very wavy. This waviness only reflects the
degree of difficulty the code encountered and is not part of the real solution. In this
case, refining the mesh towards the drying surface will only increase the waviness of
the solution. Drawing from this observation, coarsening the mesh will also get you past
this initial numerical difficulty. Although this condition may pose numerical stability
problems initially, it does not affect subsequent solution. And most of the time, one is
not interested in the solution from the initial time steps.

*Problem* – Trouble in converging in the plastic region.

*Solution* – Reduce the time step size because viscoplasticity is in itself a time
dependent problem and elasticity in itself is not. Before the material yields, time
dependency is induced only through the drying process. The reduction in time step size
depends on the value of the plastic viscosity. The lower the viscosity, the small time
step should be used. Also, it takes more iterations to converge a time step in the
viscoplastic region than the elastic region, so increasing the maximum allowable
iterations per time step will help.

Other Cautions:
Always set the *MASS_FRACTION* in the input file to be the same as the *Stress
Free Solvent Vol Frac* in the material file.

The code has been tested for a wide range of initial solvent volume fractions
(up to 0.85). When using very high initial solvent volume fractions
(approaching 0.85 or beyond), use with caution.

--------------
**References**
--------------

GT-019.1: Elastoviscoplastic (EVP) Constitutive Model in GOMA: Theory, Testing,
and Tutorial, P. R. Schunk, A. Sun, S. Y. Tam (Imation Corp.) and K. S. Chen, January
11, 2001

S.Y. Tam’s thesis: “Stress Effects in Drying Coatings,” Ph.D Dissertation, University of
Minnesota, 1997

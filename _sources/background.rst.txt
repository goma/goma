==========================
Background Information
==========================

Program Features
####################


Free and Moving Boundary Capabilites
****************************************

*Goma* is a general purpose program designed for the solution of both steady and transient, two and
three-dimensional problems involving heat, mass, and momentum (solid and fluid) transport.
A unique feature is the treatment of all boundaries and interfaces as *free* (position unknown) or
*moving* (position unknown or prescribed, but variable). If the material domain of interest is a
solid, a Lagrangian formulation (i.e., the computational mesh follows the motion of material) of
the momentum equations naturally leads to mass conservation and a natural parameterization of
the boundaries and interfaces as material surfaces. If the material domain of interest is a fluid,
then an Arbitrary-Lagrangian-Eulerian (ALE) formulation allows the boundaries to respond to
constraint equations, hereafter referred to as *distinguishing conditions*. These conditions are
responsible for determining the location of all boundaries and interfaces, providing the necessary
mathematical closure of the system of equations governing the free boundary problem.
Distinguishing conditions available to the user fall into several classes, as described below.

Since publication of the *Goma* 2.0 manual in 1998 (and more recently the *Goma* 4.0 manual in
2002), the ALE formulation has been extended to solid-material regions (viz. the TALE
algorithm, Schunk, 2000) and purely Eulerian front tracking schemes based on the method of
level-sets have been incorporated for free surfaces with large deformations; moreover, both
algorithms have been implemented in a completely-coupled way. Of course Eulerian schemes are
inherently transient and less accurate in capturing interfacial physics, even though they are more
robust and even optimal for a certain class of problems. It is fair to say that of all the available
mechanics codes, *Goma* provides the greatest breadth of free and moving boundary tracking
formulations and options.

With regard to the ALE algorithms, the fully-implicit, pseudo-solid, unstructured mesh
deformation algorithm sets *Goma* apart from other finite element programs. All surfaces, internal
and external, together with other geometric features such as corners and junction points, are
permitted to move as part of the algorithm. The movement of boundaries, interfaces, and
geometric features is dictated by a weighted residual statement of the distinguishing conditions,
whether based on specific physical constraints or arbitrary conditions described by the analyst.
The internal mesh deforms as if it were embedded in a deforming elastic solid continuum; with
the mechanics of the solid governed by either infinitesimal (linear) or finite (nonlinear)
deformation theory. Through Newton’s method, the deformation of the mesh is determined
simultaneously with all of the other physics of the problem.

The key connection between the mesh deformation and the physics of interest is accomplishedthrough a library of distinguishing conditions. Currently, those conditions include (a) kinematic
(material surface of a fluid), (b) isotherm (phase transition temperature, such as melting), (c) isoconcentration
and (d) geometric (either smooth plane curves or fixed point specifications). As part
of the required input for *Goma*, the analyst specifies the associations between the particular
distinguishing conditions and corresponding sets of material points of the initial pseudo-solid
used to embody the mesh. Chapter 4 describes this process in more detail. Essentially, the
algorithm causes smooth boundaries of the pseudo-solid to slide tangentially in a “frictionless”
fashion. Further details of this algorithm and the corresponding equations can be found in several
references (e.g., Sackinger, Schunk, and Rao, 1995).

Coordinate Systems and Frames of Reference
**********************************************

Coordinate systems accessible through this version of *Goma* include two-dimensional and threedimensional
Cartesian coordinates, cylindrical coordinates for axisymmetric problems, spherical
coordinates, and a swirling option for two-dimensional axisymmetric problems with a (swirling)
velocity component in the third dimension. A limited framework has been built within *Goma* to
use arbitrary orthogonal curvilinear coordinate systems, but this has not yet been extensively
tested. As for frame of reference, all conservation equations are cast in an inertial frame (viz. nonaccelerating)
but with extensions to allow for arbitrary frame velocities that may or may not be
related to the material motion. Hereafter, when we refer to the frame/mesh motion type to be of
the *Eulerian* variety, we mean the mesh is fixed with respect to all material motion, which
basically means it is fixed in the laboratory frame. For now, we allow this frame of reference for
fluid systems and are researching ways to allow this frame for solid systems. The ALE frame of
reference, as mentioned above, allows for independent mesh motion in the interior of the domain,
but seeks to maintain a material frame of reference on the boundary. This means that the mesh
will move to accommodate material boundary motion. Currently, the ALE frame is allowed for all
classes of materials (cf. Schunk, 2000). Finally, a pure Lagrangian frame of reference implies that
our mesh moves with the material. This formulation is quite common in solid mechanics and is
one advocated here for truly solid regions.

Problem Physics and Thermophysical Properties
*************************************************

This brief section summarizes the physics capabilities in **Goma** and the thermophysical properties
and constitutive equations available to the user. The rest of the manual is designed to greatly
expand on all material parameter options, boundary condition options, and equation options;
perusing Chapter 4 and Chapter 5 is recommended to extract more detail.

The class of problems treated by **Goma** are those described by any one or a combination of the
incompressible form of the momentum conservation equation for generalized Newtonian fluids,
the momentum conservation and differential stress constitutive equations for viscoelastic fluids,
saturated and unsaturated flow equations cast for rigid or deformable porous media, the energy
conservation equation, the equations of quasi-static equilibrium of an elastic solid, and any number of additional or auxiliary species convection-diffusion-reaction equations. **Goma** has
been tested with the following types of fluid mechanics, solid mechanics, and heat transfer
problems: (a) mixed convection with mesh parameterization of an isotherm, (b) melting, with a
parameterization of the liquidus and solidus isotherms, (c) coating and related flows (slide
coating, curtain coating, etc.), (d) polymer processing (viscoelastic) flows (e.g. fountain flow,
planar and axisymmetric extrusion, simple mold filling, contraction flow), (e) neutral or charged
species transport in multicomponent concentrated systems, (f) partially saturated flow in
poroelastic systems, (g) suspension flows, (h) drying and shrinking of gelled polymer films (with
creep and elastic recovery), and (i) microfluidic systems with fluid-structure interaction (e.g.
MEMS device performance).

Thermophysical properties in the bulk for all equations may be taken as constant or variable, with
dependencies on any of the dependent and independent variables of the problem. General
property variation models of this sort can be implemented with a user-defined subroutine
capability. Moreover, a growing number of often-used standard models are supported within the
core routines. These include a Carreau-Yasuda model for the generalized Newtonian viscosity and
a Boussinesq source term for the fluid momentum equation that provides a means for simulating
flows with thermal and solutal buoyancy. A plethora of other constitutive models and properties
are available, including viscoelasticity, elastoviscoplasticity, nonFickian diffusivity, etc.

To enhance the capability for modeling problems in capillary hydrodynamics, e.g., coating flows,
a boundary condition expressing the normal stress balance for two-dimensional Cartesian and
axisymmetric problems has been implemented and tested. When capillary forces are activated, a
pressure jump term (proportional to the mean curvature) is added to the normal component of the
momentum flux balance at specified fluid material interfaces in a natural fashion. At three-phase
boundaries (points in two dimensions) a contact angle condition and a surface tangent force
condition may be applied. The former is used in place of a specified position on the mesh motion
equations and is best used to set static and dynamic contact angles, and the latter is an additional
endpoint force which is added to the momentum balance, necessitated because the curvature term
is integrated by parts. The current version of *Goma* also includes the ability to model tangential
shear forces along capillary surfaces, i.e., those originating from surface tension gradients caused,
for example, by variations in temperature or species concentration. To access this capability
requires a constitutive equation for the surface tension. A powerful low-level capability has been
implemented which allows the user to select which degree of freedom, or variable, is associated
with a particular boundary condition. Such a capability is useful at dynamic contact lines, where it
is often desirable to replace the liquid-phase momentum equations with auxiliary constraint
conditions.

Generalized interphase boundary conditions that allow for discontinuous field variables are
supported through a multiple degree-of-freedom capability. The prime targets for this capability
include flowing vapor-liquid equilibrium problems for which there are concentration and velocity
jumps between phases due to change in density and solute partitioning through the phase diagram
and multiphase/multicomponent corrosion problems. A series of boundary conditions which allow for the application of ideal and non-ideal vapor/liquid equilibrium (e.g. Raoult’s law and
Flory-Huggins theory), latent heat release/adsorption, and discontinuous velocity components due
to evaporation/condensation have been implemented. In the future this capability can be extended
to thermal contact resistance, which often involves a temperature jump at an interface.

Recently the solid mechanics module of *Goma*, which was originally installed as a part of the
pseudo-solid ALE mesh motion algorithm, has been exploited to solve problems in transport in
deformable porous media and other outstanding problems of elastohydrodynamics. For modeling
flow in non-deformable porous media, the Brinkman terms in the fluid momentum equations (cf.
Gartling, et. al., 1996) may be activated. Since *Goma* 2.0, generalized Darcy transport equations
for multiphase components (solid, liquid, gas) have been added and can be used for simulations of
deformable poroelastic media. For incompressible but deformable solids, a pressure term was
added to the solid momentum balance (e.g. rubber). In continuous shrinking or swelling solids,
the dilation is proportional to changes in solvent concentration. In deformable porous media, the
solid deformation is coupled to the pressure in the fluid-filled interstices of the porous matrix.
Several boundary conditions exist to apply normal tractions (i.e. compressive, tensile, or shear
boundary forces) to solid surfaces. To effectively simulate coupled fluid/solid interaction
problems, boundary conditions which balance the surface tractions exerted by the liquid and solid
phases at the common interface have been incorporated as have been the appropriate interface
impregnation/expulsion conditions at boundaries between porous and continuous media.

A complete rewrite of the species transport equations has been undertaken since the release of
*Goma* 2.0 that allows for generalized phase/species formulations on multimaterial problems.
Accommodating an arbitrary number of species, each of which can exist in an arbitrary number of
phases, was the goal of this development in order to model corrosion and charged species
transport.

Of course there are many more material property models and constitutive equations, specialized
boundary conditions, and more esoteric differential equations that can be solved for just about any
mechanics problem. Many of these capabilities are not cited in this manual because they were
under development at the time of publication. Interested readers should inquire about the status of
the following capabilities: generalized solid-model geometry features, wetting and spreading
models for Eulerian front tracking schemes, Eulerian/Eulerian fluid-structural interaction
capability, multiphase porous energy equation, Generalized surface and volume user-defined
Lagrange multiplier constraints, and much more.

Advanced Capabilities
*************************

Several developments in *Goma* that enable advanced engineering analysis of complex systems
have been completed since the last major release. These developments include a complete,
generalized capability of automated parameter continuation (zeroth-order, first-order, arclength,
multiparameter, user-defined parameter continuation, etc.) using the LOCA library (Salinger, et.
al., 2002), linear stability analysis of any dynamic system using normal modes, and augmenting condition capability. It is recommended that the user consult a separate manual (Gates et. al.,
2000; contact authors for a more recent version) for a complete user description of these features.
The input record sections required to activate these features are not covered in this document.

Numerical Methods
#####################

With over 150 different boundary conditions for 70 plus differential equation types, *Goma’s*
algorithms are very extensive for any brief discussion. In this section we simply point out the
foundation algorithms. A developer’s manual, advanced capabilities manual, and tutorial memos
can be consulted for more details (see *Goma* Document List in the Appendix for the citations.).

*Goma* is based primarily on the Galerkin/finite element method. The element library currently
includes (in two dimensions) 4- and 9-node isoparametric quadrilaterals (i.e., Q1 and Q2
interpolations) with available interpolations for linear discontinuous (P1) or piecewise constant
(P0) variables, and (in three dimensions) 8-node isoparametric hexahedral elements and 27-node
bricks, also available with piecewise constant interpolations. The overall solution algorithm
centers around a fully-coupled Newton-Raphson iterative scheme for solving the nonlinear
algebraic equations which results from the finite element discretization. That is, all active
equations and boundary conditions are solved simultaneously in a single matrix system at the
same time plane and during the same Newton iteration. The sparse matrix system is stored in a
central element-level matrix data structure that is injected into one of three sparse matrix formats
as dictated by the matrix solver chosen. The three formats are modified sparse row, MSR or
compressed row format (Hutchinson, et. al., 1995, Schunk and Shadid, 1992), the variable block
row, or VBR, format (see Heroux, 1992), or the frontal-solver element-level format (cf. Hood,
1976). If the matrix system is not too poorly conditioned, then iterative solvers of the generalized
preconditioned conjugate gradient-type can be used to solve the system (see Tuminaro, et. al.,
1999, Schunk and Shadid, 1992). A new matrix-services/solver-services library known as
TRILINOS (http://www.cs.sandia.gov/Trilinos), has been installed to handle all iterative solver
and preconditioner options. This package has greatly extended the robustness of iterative solvers
to the class of problems that *Goma* solves. Virtually all problems and all finite element
formulations are now solvable with these iterative schemes (see Schunk, et al., 2002). If all else
fails, *Goma* deploys a suite of direct solvers that, even though not always efficient for large threedimensional
problems, will always get a solution at the current Newton iteration. These solvers
are known as Sparse 1.3 (lu), a classical LU decomposition (Gaussian elimination) method, and
two frontal solvers, Umfpack (umf) and front; these are discussed in the next section.

The Galerkin least squares (GLS) method for pressure stabilization of Hughes and Franca (1987)
has also been added to *Goma*. The GLS method adds the momentum residual, weighted by the
gradient of the Galerkin weight function, to the standard Galerkin continuity equation, thus
providing a diagonal term for the pressure. This is a first-order convergent and consistent method
that enables the use of iterative solvers for incompressible equations over the entire range of
Reynold’s numbers.

The overall differential-algebraic system of equations may be advanced in time with implicit
time-integration techniques (simple backward Euler and Adams-Bashforth predictor, trapezoidal
corrector algorithms for fluid systems, species transport and energy transport; and Newmark-Beta
algorithms for solid dynamics). Time marching offers an alternative, albeit indirect, route to
attaining solutions to steady equations, as well as providing the capability of simulating process
transients directly. Automatic time step control based on current truncation error is also available.

Perhaps the most complicated part of the algorithm is the construction of the Jacobian sensitivity
matrix. Because the mesh point positions are actually unknowns in a free or moving boundary
problem, that matrix must include sensitivities of each weighted residual equation with respect to
each of the mesh variable unknowns that can affect the value of the residual. Unfortunately,
almost every term of the bulk equations and many boundary conditions contribute to this
sensitivity. This occurs mainly through gradient operators and surface normal and tangent vectors
(see Kistler and Scriven, 1983) and through dependencies on mesh position of the determinant of
the elemental Jacobian transformation matrix that maps between a fixed unit element and any
element in the computational domain. Great care has been taken to include analytical expressions
for all of these mesh sensitivities. However, some of this task inevitably falls to the user when
implementing user-defined boundary conditions, material property models, and constitutive
equations, particularly when any of these quantities depends directly on spatial position or spatial
gradients of other variables. In order to maintain the strong convergence properties of Newton’s
method, these sensitivities must be specified in those user-defined routines. To aid in this task, a
debugging option is available which computes a numerical finite-difference approximation of the
global Jacobian matrix and compares it with its analytical counterpart. This tool enables users and
developers to check the consistency of newly-created equations (whether bulk or boundary
constraints) with their corresponding analytic Jacobian contributions.

Portability, Software Library Infrastructure, and Code Accessibility
########################################################################

*Goma* is written in the C programming language (specifically Kernighan and Ritchie, 1988, C
with some ANSI extensions). It has been ported to a number of UNIX platforms including Solaris
and Linux, with the Linux Enterprise-4 version being the most actively maintained. Most recent
versions are aimed at Red-Hat RHEL5 and RHEL6 levels, almost exclusively. Many of the
machine dependencies in the program have been isolated using C preprocessor directives. Some
of the machine dependencies that occur in the I/O routines are insulated from the user by software
libraries. Building *Goma* requires EXODUS II v2.02 (Schoof and Yarberry, 1994), SPARSE 1.3
(cf. Kundert and Sangiovanni-Vincentelli, 1988), NetCDF v2.3.2 (Rew, et. al., 1993) libraries,
Umfpack direct solver libraries (Davis and Duff, 1997), and the TRILINOS 10.0 library
(Tuminaro, et. al., 1999; **http://software.sandia.gov/trilinos**). The first of these is part of the
SEACAS system at Sandia National Laboratories (Sjaardema, 1993); the latter two libraries are
available publicly. Parallel processing is enabled by OPEN-MPI. The user should consult the
build instructions for the most recent library revisitions. The most updated library needs are also made clear in the *Goma* makefile: Makefile. There are special versions of this makefile for
building for the test suite (Makefile_guts) and debug mode (Makefile_debug). These are
the most general makefiles that are deployed. Generally, pre- and post-processing is performed
outside of *Goma*, although some post-processing of results is available within the program. This
separation of the functionality permits the use of alternative solid-modeling and mesh-generation
software and visualization packages of choice, insofar as they may be interfaced with the
EXODUS II finite element data model.

Pre-processing options include mesh generation via CUBIT (**http://cubit.sandia.gov**), PATRAN
(PDA, 1990), and SolidWorks (**www.solidworks.com**). The latter two require special plug-ins.
These mesh generators currently support and will output a finite element database in the
EXODUS II format.

Post-processing options include BLOT (see the SEACAS distribution, Gilkey and Glick, 1989),
Paraview (**www.paraview.org**), and Ensight (**www.mscsoftware.com.au/products/software/cei/
ensight**).

Since *Goma* is built around the EXODUS II finite element data model, there are numerous
options available for communication with other analysis codes that also exchange data via the
same EXODUS II data model. Recent modifications to *Goma* permit not only the initialization of
unknown values from an EXODUS II file, but also the ability to incorporate field variables into
the analysis that are not unknowns. For example, the quasi-static and dynamic electromagnetic
fields from codes such as ALEGRA can be used to compute electric fields and current fluxes on a
specified finite element mesh that are input to *Goma* through the EXTERNAL FIELD data card.

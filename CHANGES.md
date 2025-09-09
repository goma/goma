## Changes in version 7.10.1

* 7.10.1 bugfix for Sacado disabled builds by @wortiz in https://github.com/goma/goma/pull/545

## Changes in version 7.10.0

### Highlights

* MUMPS solver is available directly
* Advanced capabilities manual is included in sphinx documentation
* Fluidity thixotropy model
* Work was done to support Trilinos removing Epetra based packages
* Logarithmic interpolation of level set properties is now available in general level set fields

### Included pull requests
* Gas egress by @rbsecor in https://github.com/goma/goma/pull/494
* Add logarithmic interpolation of properties across level set interface by @tjiptowi in https://github.com/goma/goma/pull/497
* Update github workflows and tpl packages by @wortiz in https://github.com/goma/goma/pull/498
* Add RPATH to CMAKE packages by default, fix no parmetis build by @wortiz in https://github.com/goma/goma/pull/499
* Add Fluidity model to species by @wortiz in https://github.com/goma/goma/pull/500
* Add Herschel Bulkley Papanastasiou model by @wortiz in https://github.com/goma/goma/pull/501
* Add some more options for Herschel Bulkley regularization by @wortiz in https://github.com/goma/goma/pull/502
* tiny fix for documentation examples by @wortiz in https://github.com/goma/goma/pull/503
* 2D Moving wall shear-thnning lubrication by @rbsecor in https://github.com/goma/goma/pull/496
* ls interface print for 3d by @wortiz in https://github.com/goma/goma/pull/504
* bugfix for printing wrong number of unknowns (tpetra/epetra) by @wortiz in https://github.com/goma/goma/pull/505
* Disable ls->Contact_Tolerance by default and bug fix for replacing Strong BC's instead of penalty by @wortiz in https://github.com/goma/goma/pull/506
* update for clang-format-18 by @wortiz in https://github.com/goma/goma/pull/507
* add missing POR_SINK_MASS check by @wortiz in https://github.com/goma/goma/pull/510
* Native MUMPS interface by @wortiz in https://github.com/goma/goma/pull/511
* add guard for no mumps solver by @wortiz in https://github.com/goma/goma/pull/512
* Compiler warning cleanup by @wortiz in https://github.com/goma/goma/pull/513
* some tpl updates by @wortiz in https://github.com/goma/goma/pull/508
* Read element variables from multiblock materials by @wortiz in https://github.com/goma/goma/pull/514
* small bugfix for stratimikos when using augmenting conditions by @wortiz in https://github.com/goma/goma/pull/516
* Lub_curv BC cleanup by @rbsecor in https://github.com/goma/goma/pull/517
* Update some TPLs and fixes for GCC 15 and building static libraries by @wortiz in https://github.com/goma/goma/pull/518
* Lubrication dynamic_contact_angle by @rbsecor in https://github.com/goma/goma/pull/519
* Improved facet based reinitialization of level sets by @wortiz in https://github.com/goma/goma/pull/509
* Add a cure shrinkage model to goma  by @wortiz in https://github.com/goma/goma/pull/521
* Making Aztec/Amesos/Epetra optional for upcoming releases of Trilinos by @wortiz in https://github.com/goma/goma/pull/488
* Autodiff version of porous shell saturation equations by @wortiz in https://github.com/goma/goma/pull/522
* Add a patch for mapvar IEEE_DENORMAL issues until next SEACAS release by @wortiz in https://github.com/goma/goma/pull/525
* Update mesh equation example in docs by @wortiz in https://github.com/goma/goma/pull/524
* force O1 optimization for SEACAS fortran by @wortiz in https://github.com/goma/goma/pull/527
* add DOUBLE_FILLET boundary condition by @wortiz in https://github.com/goma/goma/pull/528
* Trying to merge Chance's drying work by @wortiz in https://github.com/goma/goma/pull/515
* Add augmenting continuation cards to docs by @wortiz in https://github.com/goma/goma/pull/531
* Remove some set but unused variables by @wortiz in https://github.com/goma/goma/pull/533
* Add DOUBLE_FILLET_GEOM_BASED, move functions out of bc_colloc.c by @wortiz in https://github.com/goma/goma/pull/534
* Revert "Add DOUBLE_FILLET_GEOM_BASED, move functions out of bc_colloc.c" by @wortiz in https://github.com/goma/goma/pull/535
* Update SEACAS package, patch is upstreamed by @wortiz in https://github.com/goma/goma/pull/537
* Pixel to mesh fixes by @wortiz in https://github.com/goma/goma/pull/529
* DOUBLE_FILLET_GEOM and option to anneal every continuation print for remeshing scripts by @wortiz in https://github.com/goma/goma/pull/541
* More advanced capability manual being added to main docs by @wortiz in https://github.com/goma/goma/pull/542

## Changes in version 7.9.0
* Add Giesekus model to SQRT_CONF by @wortiz in https://github.com/goma/goma/pull/486
* Add skip ssl check to Goma tpl build script by @wortiz in https://github.com/goma/goma/pull/487
* Kelvin-Voigt solid model fix by @rbsecor in https://github.com/goma/goma/pull/485
* set LD_LIBRARY_PATH in configs when building TPLs with build script by @wortiz in https://github.com/goma/goma/pull/489
* SHELL_CONC_LS_BC by @rbsecor in https://github.com/goma/goma/pull/490
* Addition of Taylor Galerkin method to the moment population balance equations by @trkenne in https://github.com/goma/goma/pull/491
* Stratimikos bugfix by @wortiz in https://github.com/goma/goma/pull/492
* Update tpl packages for BLAS/LAPACK issues by @wortiz in https://github.com/goma/goma/pull/493

## Changes in version 7.8.1
- Goma v7.8.1 bugfix for viscoelastic flows by @wortiz in #483

## Changes in version 7.8.0
- fix readme links by @wortiz in #480
- Lubrication LS Heaviside Shift by @rbsecor in #481
- Goma 7.8.0: Add SHEAR_STRESS_APPLIED by @wortiz in #482

## Changes in version 7.7.0

Major new features are Autodiff and stabilization improvements and a new build script in `tpls/` folder

- Shell Energy Update by @rbsecor in https://github.com/goma/goma/pull/459
- write P1 element variables in segregated solve by @wortiz in https://github.com/goma/goma/pull/460
- Use Amesos2 Solver Package instead of just Mumps by @wortiz in https://github.com/goma/goma/pull/461
- GCC 14 warning cleanup, fix for xdot in level set segregated by @wortiz in https://github.com/goma/goma/pull/462
- Add Fluid Stress Post Processing, update Viscous Stress by @wortiz in https://github.com/goma/goma/pull/463
- tiny doc update by @wortiz in https://github.com/goma/goma/pull/464
- Shell Lubrication temperature and convection additions by @rbsecor in https://github.com/goma/goma/pull/465
- Add under relaxation and tolerances to subcycling. Fix some warnings. by @wortiz in https://github.com/goma/goma/pull/466
- Legacy build script changes for automation by @wortiz in https://github.com/goma/goma/pull/467
- Lubrication Logarithm Interpolation by @rbsecor in https://github.com/goma/goma/pull/468
- Add strumpack to legacy build script petsc, update flux docs by @wortiz in https://github.com/goma/goma/pull/471
- Updates for cmake to find non-system TPLs by @wortiz in https://github.com/goma/goma/pull/469
- Add a new TPL build script by @wortiz in https://github.com/goma/goma/pull/473
- AutoDiff and Stabilization and Line Search by @wortiz in https://github.com/goma/goma/pull/472
- Lub log ls2 by @rbsecor in https://github.com/goma/goma/pull/470
- Discontinuous Variable Initial Guess fix and Discontinuous Stress by @wortiz in https://github.com/goma/goma/pull/474
- Add current AC_VALUE to aprepro AC flag by @wortiz in https://github.com/goma/goma/pull/475
- Fix level set renormalization for quadratic tets by Harald Ziegelwanger in https://github.com/goma/goma/pull/477


## Changes in version 7.6.1

- bugfix for Trilinos 13 [goma#458](https://github.com/goma/goma/pull/458)

## Changes in version 7.6.0

- Tpetra and Amesos2 support, currently expects `TPETRA_INST_INT_INT` [goma#456](https://github.com/goma/goma/pull/456)
- `MULTI_CONTACT_LINE` Lame Mu model for multi CA problems (also 3D) [goma#454](https://github.com/goma/goma/pull/454)
- Output BC dup in 3D  [goma#452](https://github.com/goma/goma/pull/453)
- White Metzner type models [goma#450](https://github.com/goma/goma/pull/450)
- Support for Trilinos 14 and 15 [goma#442](https://github.com/goma/goma/pull/442)
- Elliptic mesh for single materials [goma#440](https://github.com/goma/goma/pull/440)

## Changes in version 7.5.0

- Addition of Rolie-Poly Viscoelastic Model [goma#433](https://github.com/goma/goma/pull/433)
- Spalart Allmaras turbulence model [goma#432](https://github.com/goma/goma/pull/432)
- Various cleanup [goma#431](https://github.com/goma/goma/pull/431)
  [goma#430](https://github.com/goma/goma/pull/430)

## Changes in version 7.4.5

- Bugfix for CA conditions [goma#429](https://github.com/goma/goma/pull/429)

## Changes in version 7.4.4

- Bugfix for parallel decomposition [goma#427](https://github.com/goma/goma/pull/427)

## Changes in version 7.4.3

- Element quality exits consistently with LOCA continuation [goma#424](https://github.com/goma/goma/pull/424)
- Various bugfixes/Cleanup [goma#424](https://github.com/goma/goma/pull/424)
  [goma#421](https://github.com/goma/goma/pull/421)
  [goma#420](https://github.com/goma/goma/pull/420)
  [goma#419](https://github.com/goma/goma/pull/419)

## Changes in version 7.4.2

- Online documentation fix, broken with newer Sphinx [goma#415](https://github.com/goma/goma/pull/415)

## Changes in version 7.4.1

- Bugfix for Goma builds without PETSc [goma#414](https://github.com/goma/goma/pull/414)

## Changes in version 7.4.0

- Add Nedelec Hcurl first order Basis functions [goma#413](https://github.com/goma/goma/pull/413)
- EM formulation based on Nedelec elements [goma#413](https://github.com/goma/goma/pull/413)
- Absorbing boundary condition [goma#413](https://github.com/goma/goma/pull/413)
- Post processing for EM equations [goma#413](https://github.com/goma/goma/pull/413)
- PETSc complex solver support [goma#413](https://github.com/goma/goma/pull/413)
- dpi fix for GCC 7 [goma#412](https://github.com/goma/goma/pull/412)
- KIN_DISPLACEMENT Jacobian entries [goma#411](https://github.com/goma/goma/pull/411)
- SUPG Jacobian entries [goma#410](https://github.com/goma/goma/pull/410)

## Changes in version 7.3.0

- Add SQRT form of Viscoelastic conformation tensor https://doi.org/10.1016/j.jnnfm.2011.02.008 for 2D, 3D, Axisymmetric OldroydB, PTT-EXPONENTIAL, PTT-LINEAR
- Add KDR model from @jtmcconnell https://doi.org/10.1103/PhysRevLett.126.218002
- Add `FLOW_GRADV_T` to impose $n\cdot\nabla v=0$ naturally (probably what `FLOW_GRADV` should've been)
- Update `STRESS_DEVELOPED` to use new forms
- Fix `LOG_CONF` numerical Jacobian
- $YZ\beta$ shock capturing for `SQRT_CONF`


## Changes in version 7.2.0

- Add basic aprepro library support, so aprepro isn't needed on path for most tasks
- Add Augmenting Condition on multple Boundary condition parameters using aprepro
  see [https://github.com/goma/goma/pull/406](https://github.com/goma/goma/pull/406)

## Changes in version 7.1.3

- Komplex support now compiles, untested

## Changes in version 7.1.2

- bugfix for LOCA parallel
- bugfix LOCA global variables exodus output
- add `ENABLE_KOMPLEX` to CMakeLists.txt

## Changes in version 7.1.1

- bugfix for linear PTT BC `STRESS_DEVELOPED`
- bugfix garbage malloc for `fix` program
- bugfix for RCB specification in input card


## Changes in version 7.1.0

- Add experimental standalone fix, usage: fix <exodusfile>.<num_proc>.00
  - example: `fix out.exoII.24.00`
- Linear PTT model, see PTT Form
- bugfix, CA conditions error checking
- bugfix, Output exodus properties, used with Nodeset associativity in CUBIT
- bugfix, iapply is always on on external surfaces for "Special" BCs
- build script contains fix for goma#346 issue
- bugfix, consistent NS and SS when using METIS decomposition

## Changes in version 7.0.5

- bugfix, fix anneal mesh not working
- bugfix, add `ENABLE_LOGGING` and fix `mm_eh.c` includes when `ENABLE_LOGGING`

## Changes in version 7.0.4

- bugfix for External Field Variables in Parallel

## Changes in version 7.0.3

- bugfix for VE Shock Capturing documentation and Parallel

## Changes in version 7.0.2

- Build script now builds `Omega_h` and `PETSc`
- Various `CMake` changes to address consistency of setting third party libraries
- Fix for Numerical Jacobian log conformation with multi-material

## Changes in version 7.0.1

- Bugfix for Ubuntu OpenMPI and using build script on Ubuntu

## Changes in version 7.0

- Partial Segregated solver support
    - Some developer notes are located `docs/segregated_notes.org`
    - Small memo on segregated solves is available on website
    - Mostly tested the following equations
        - Momentum + Continuity
        - Level Set
        - Energy
        - Species
        - Viscoelastic flow
        - Moments
        - Mesh
    - If you have need for segregated solver support on other problems or find a bug
      open an issue
- Automatic rotations for normal conditions in 3D
    - Requires that no ROT conditions are specified and only NORMAL rotated BCs are used
    - Works for Momentum and Mesh
    - Does not work for TANGENT or EDGE BCs, or if ROT cards are used
- Quadratic tetrahedrons (TET10 in Cubit)
- Partial PETSc support for matrix solvers
    - Can be used with a `petscrc` or `-petsc "petsc options"` on the command line
    - See manual `Solution Algorithm` for more information
- Various new equations, including:
    - Population balance equations using Moments (see Ortiz et al. 2022 AIChE)
    - Time-Harmonic Maxwell Equations with Lagrange Finite Elements and a penalty stabilization
- Adaptive meshing using `Omega_h` library for first order simplices
    - Requires first order Tet or Tri
    - Currently mostly used for level set problems see `Level Set Adapt` cards
    - ALE adaptivity of uniform ("iso") mesh size, see `ALE Adapt` cards
    - Currently limited to single materials
    - Nodesets and Sidesets must contain the same geometry and have the same number
- Brk/Fix support has been removed
    - Now compatible with SEACAS tools using elemental decomposition
    - epu from SEACAS replaces fix
    - brkfiles are no longer supported
    - See SEACAS decomp / `nem_slice`, `nem_spread` for more controlled decomposition
    - Automatic fix is still included, Fix Frequency is still available though not recommended
    - New output format is compatible with Paraview without `fix`ing
    - Automatic basic decomposition is available when configured with METIS
        - See cards `Decomposition Type`, `External Decomposition`
        - Command line flags `-e`, `-external_decomp`, `-kway`, `-rcb`
        - Defaults to `rcb` when less than 8 processors and `kway` when 8 or more
- Various bug fixes and feature updates, see merged Pull Requests
- SUPG formulation for Viscoelastic `EVSS_F`/`LOG_CONF` and Species has been modified
- PSPG has an additional formulation `pspg_shakib`
- Level set has additional SUPG formulations `SUPG_GP` and `SUPG_SHAKIB`

## Changes in version 6.2

- Many bug fixes
- Stratimikos support
- Log conformation tensor stress model
- Hysing and Denner surface tension models for level set
- Suspension balance updates
- Updated SUPG for species
- Quadratic triangles
- And more...

## Changes in version 6.1

- Automatic creation of brk files (see `docs/parallel_integration.md`)
- Bug fixes (mostly MPI related bugs)
- Epetra is now a supported matrix format (see `docs/epetra_matrix_format.md`)
- MUMPS is now supported through Amesos
- An experimental build script to build library dependencies is available under `scripts/`

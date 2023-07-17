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

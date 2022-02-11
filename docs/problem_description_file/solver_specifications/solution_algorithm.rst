**********************
Solution Algorithm
**********************

::

	Solution Algorithm = {char_string}

-----------------------
Description / Usage
-----------------------

This required card selects an algorithm for the solution of the linear matrix system that
arises at each Newton iteration (either for a steady-state solution or for the solution at
each discrete time). Please note that at the time of this writing, new solver capabilities
were being generated; although the following information was complete and accurate,
it will likely be out of date by the time of publishing. Users should consult the CD
version of this document in the Goma Documentation System for up to date options.

There are three major matrix solver packages accessible in *Goma*, two direct
factorization collections and an iterative solver package. The first collection of direct
factorization methods in *Goma* include the Sparse1.3 package (Kundert and
Sangiovanni-Vincentelli, 1988) and Y12M direct factorization technique (Zlatev,
Wasniewski and Schaumburg, 1981) accessible via the Aztec linear solver package.
The second collection of direct factorization methods include two frontal solvers,
SNL_MPFRONT, an adaptation of R. Benner’s implementation of Hood’s (1976)
frontal method, and UMFPACK (Davis and Duff, 1997). SNL_MPFRONT is a
traditional frontal method while UMFPACK is a multi-frontal solver.

The Aztec 2.x linear solver package (Tuminaro, et. al., 1999) is the iterative solver
component of *Goma*. A successor to the krysolve 1.0 package (Schunk and Shadid,
1992) and the Aztec 1.0 package (Hutchinson, Shadid and Tuminaro, 1995), Aztec 2.x
includes support for distributed memory architectures and for matrices in either a
modified sparse row (MSR) format or a variable block row (VBR) format, as well as
their distributed memory extensions. Generally, convergence of these iterative methods can be accelerated by judicious use of a preconditioner (which many of the other *Solver
Specifications* cards address).

The options for this input card are listed below, but additional usage comments are
included as part of the Technical Discussion section of this card. These comments
provide assistance in choosing the *Solution Algorithm* for your problem.

Valid options for {char_string} are as follows:

lu
    Direct factorization via Gaussian elimination using Sparse 1.3. This solver
    is robust even for poorly conditioned matrix systems. It is unavailable
    when running *Goma* on multiple processors.
front
    Direct factorization based on Benner’s SNL_MPFRONT that eliminates
    equations and variables as the fully assembled rows of the matrix are
    acquired. This is the latest solver installed within *Goma* and users are
    encouraged to report their successes and failures with this option as part
    of testing. It is unavailable when running *Goma* on multiple processors.
umf/umff
    Direct factorization using UMFPACK. This multi-frontal solver has been
    hardwired to perform elimination only upon complete assembly. The **umff**
    option forces a full factorization every time, whereas **umf** does not. It
    is unavailable when running *Goma* on multiple processors.
y12m
    Direct factorization using the Y12M package. This package is accessible
    through the Aztec matrix solver interface and cannot be used for multiple
    processor computations. Other direct solvers are recommended against this
    one.
gmres
    Iterative solver from the Aztec package using the restarted generalized
    minimum residual method. Iterative solver options are important to
    convergence of this method, e.g. *Preconditioner, Size of Krylov subspace,
    Matrix,* etc.
cg
    Iterative solver from the Aztec package using the conjugate gradient
    method. Like other iterative solvers, the successful convergence of the
    conjugate gradient method for a linear system depends on preconditioners
    and other cards in the *Solver Specifications* section.
cgs
    Iterative solver from the Aztec package using the conjugate gradient
    squared method. Convergence of this method is frequently contingent on the
    linear system and on the choice of other cards in the *Solver
    Specifications* section.
tfqmr
    Iterative solver from the Aztec package using the transposefree
    quasi-minimum residual method. Convergence of this method is frequently
    contingent on the linear system and on the choice of other cards in the
    *Solver Specifications* section.
bicgstab 
    Iterative solver from the Aztec package using the biconjugate gradient with
    stabilization. Convergence of this method is frequently contingent on the
    linear system and on the choice of other cards in the Solver Specifications
    section.
amesos
    Allows access to direct solver options implemented in parallel. Please see
    the user-notes below for Goma build options that must be exercised. This
    package is part of the Trilinos 6.0 framework. With this option, you must
    add an additional input card to specify the parallel direct solvers:

    ::

        Amesos Solver Package = {superlu | mumps | klu | umfpack}
							  
    Of these four options, we currently recommend **mumps**.
    All options can be run in parallel.
stratimikos
    Interface to Trilino's Stratimikos package
    requires:

    ::

        Matrix storage format = epetra
    
    Allows block solvers, see also ref:`Stratimikos File`
petsc
    PETSc solver and preconditioner, will use `petscrc` file or `-petsc`
    command line, see Technical Discussion for more information


------------
Examples
------------

Following is a sample card:
::

	Solution Algorithm = lu

Another example (two cards) shows how to invoke a parallel direct solver:
::

	Solution Algorithm = amesos

::

	Amesos Solver Package = superlu

-------------------------
Technical Discussion
-------------------------

The direct factorization options are the most robust but consume the most
computational resources (CPU time and memory, particularly for large and 3D
problems). The iterative methods consume less resources but may take some
experimentation to obtain convergence to the solution of the linear system. For
example, a poorly conditioned linear system may require a lot of preconditioning. The
conjugate gradient method may not be very useful on linear systems that are not symmetric 
positive definite. Although the following guidelines are useful, selection of
the “right” linear solver requires experience, understanding and sometimes, luck.

* **lu** - The Sparse1.3 direct solver, is the most robust solver in *Goma* in terms of
  obtaining successful convergence for even poorly conditioned matrix systems. A significant 
  disadvantage, however, is that it can be computationally expensive for
  large problems. Not only do the memory and CPU requirements grow with
  problem size, but the initial symbolic factorization that seeks optimal reordering
  also consumes greater CPU resources with larger problem sizes. For example, a
  problem with 70,000 degrees of freedom that required 22 hours of CPU for the
  initial factorization required only 1/2 hour for subsequent factorizations.
  Furthermore, this solver is unavailable when *Goma* is run on multiple processors.
  Its robustness makes it an excellent choice for small- and medium-sized problems.

* **front** - This solver is an adaptation for *Goma* of R. Benner’s frontal solver, which
  itself includes considerable improvements compared to the pioneering frontal
  solvers (Irons, 1970; Hood, 1976). The SNL_MPFRONT library is compiled and linked into *Goma* 
  only by choice. Direct factorization is done as the fully
  assembled rows of the matrix are acquired. The frontal solver consumes CPU time
  roughly comparable to Sparse 1.3, with the noted advantage of eliminating intraelement
  fully summed equations as they are encountered and only keeping the
  active working matrix in-core, thereby reducing memory requirements and
  possible storage of matrix components to disk.

* **umf/umff** - UMFPACK 2.0d is a powerful direct solver that is generally faster
  than Sparse 1.3a, though it might lack the robustness of the latter on infrequent
  occasions. The implementation of UMFPACK within *Goma* is only barebones, i.e.
  the multi-frontal solver has been hardwired to perform elimination only upon
  complete assembly. Finally, usage of UMFPACK is governed by a license that limits usage to 
  educational, research and benchmarking purposes by nonprofit
  organizations and the U.S. government. Please refer to the license statement
  contained in the UMFPACK distribution for exact details. This solver was
  implemented prior to **front** so it was the only direct solver alternative to lu for a
  period of time. User’s should now evaluate performance of this solver against **front** on a case by case basis.

* **gmres, cg, cgs, tfqmr, bicgstab** - The convergence of each of these iterative
  solvers is highly influenced by the kind of preconditioning selected. Often, the
  method(s) will not converge at all without an appropriate level of preconditioning.
  GMRES is considered one of the best iterative methods available, although there
  are instances where each of the others is superior. It is a Krylov-based method and has an 
  additional input card, *Size of Krylov subspace*. As mentioned earlier, CG
  should only be used on systems that are symmetric positive definite. See the *Matrix 
  subdomain* 
  solver card, and other *Solver Specifications* cards for guidance
  on appropriate use of preconditioners; also consult Schunk, et. al. (2002).

* **amesos**: superlu, klu, umfpack - These solvers are all direct (not iterative, but
  based on Gaussian elimination) and can be run in parallel with mpi. We
  recommend these solvers when robustness is required over iterative solvers and
  when the matrix assembly time is excessive, which is often the case when
  overloaded equations like species diffusion, porous media equations, etc. are used.
  This option also performs well for three-dimensional problems of small to
  moderate size. 

* **stratimikos**: mostly used for interfacing with Trilino's `Teko` but can also call
  full solver suite that is supported in Trilinos through xml files

* **petsc**: There are quite a lot of linear solvers and preconditioners available through
  PETSc and most are configured through either command line arguments using `-petsc` or 
  using a `petscrc` file in your goma problem directory specifying petsc options

  Options are specified using the usual `ksp_type` and `pc_type` etc

  ::
    
    -ksp_type gmres
    -pc_type asm
    ... etc

  When in a segregated solve `ksp` and `pc` options should be prefixed with a 0-indexed `-sys#`
  corresponding to each matrix

  ::
    
    -sys0_ksp_type gmres
    -sys0_pc_type asm
    -sys1_ksp_type gmres
    -sys1_pc_type hypre
    ... etc


--------------
**References**
--------------

SAND2001-3512J: Iterative Solvers and Preconditioners for Fully-coupled Finite
Element Formulations of Incompressible Fluid Mechanics and Related Transport
Problems, P. R. Schunk, M. A. Heroux, R. R. Rao, T. A. Baer, S. R. Subia and A. C.
Sun, March 2002.

G. H. Golub and C. F. V. Loan, Matrix Computations, Johns Hopkins University Press,
Baltimore, MD 3rd ed. (1996)

For all other references, please see *References* at the end of this manual.

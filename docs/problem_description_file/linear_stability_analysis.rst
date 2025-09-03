Linear Stability Analysis
============================

Stability analysis addresses the following question: After disturbing a given steady base state with 
random disturbances, do the disturbances grow or decay in time? If the disturbances vanish and 
the base state persists, then the system is stable. Otherwise the disturbances grow and the state 
changes to a new state that can be physically undesirable.

Stability analysis tests the stability of a steady-state flow solution by adding disturbances to it. By 
invoking the method of small disturbances and representing a generic disturbance as a 
superposition of normal modes, the stability problem takes the form of an asymmetric and 
singular generalized eigenvalue problem. This eigenvalue problem characterizes the stability of 
the base state subjected to the imposed disturbance. The disturbance modes are the eigenvectors 
and the growth or decay rates are the eigenvalues. The sign of the real parts of the eigenvalues 
(growth rate if positive; decay rate if negative) tell whether the flow is either stable or not. If one 
or more of the real parts of the eigenvalues are positive, the base flow is unstable; otherwise the 
base flow is stable.

The stability problem arising from incompressible fluid mechanics contains a few troublesome 
features that must be taken into account. The first concerns the properties of the two, usually 
large, matrices involved in linear stability analysis, The Jacobian matrix (sensitivities of the 
residuals to the unknowns), is in general asymmetric, as are shifted combinations with the mass 
matrix (sensitivities of the residuals to time-dependent unknowns). This means that fast 
symmetric eigenvalue extraction methods are not suitable. The second is that the mass matrix is 
singular. This is because both Dirichlet boundary conditions and the incompressibility constraint 
contain no time dependence. This gives rise to 'infinite' growth (or decay) rates, since the speed 
of sound is infinite in incompressible media (Christodoulou, 1988). The 'infinite magnitude' 
eigenvalues must be filtered out so that the eigenvalue extractor finds the physically relevant, 
most dangerous eigenvalues. A simple manipulation that maps 'infinite' eigenvalues to zero is the 
shift-and-invert transformation. The 'shift' is usually set near to the anticipated values of the 
leading eigenvalues. Most often, the shifts are set to the values of the leading eigenvalues 
available at a 'close' parameter value. More details about the nature of the matrix stability 
problem are provided in Coyle (1984), Christodoulou (1988), and Gates (1998).

The method used in Goma to extract eigenvalues is Arnoldi's method (Saad, 1992). It makes use 
of projections of the shifted matrix onto a subspace spanned by iterates of the power method. The 
method belongs to the class of Krylov subspace methods. Saad (1992) provides a very thorough 
discourse on these methods. Arnoldi (1951) presented the procedure as a means for reducing 
dense matrices to Hessenberg form. The interested user can consult Saad (1992), Christodoulou 
(1988), and Gates (1998) for more details on the algorithm implemented in Goma.

There are now two eigensolvers available for Goma. Eggroll, which has been available for some 
time, is included in the source code, and ARPACK (Lehoucq and Sorensen, 1996, Lehoucq and 
Scott, 1997) can be linked in as an external library. Both eigensolvers use Arnoldi-based Krylov 
subspace methods to extract the first few (n) eigenvalues of interest; however, the LOCA driver 
for ARPACK uses either a one-parameter shift-and-invert algorithm or a two-parameter Cayley 
transformation algorithm, whereas eggroll only uses shift-and-invert. Moreover, eggroll cannot be 
accessed during continuation, and requires UMFPACK to perform the necessary solves of the 
shifted linear eigensystem (which precludes running in parallel). ARPACK can use any linear 
solver supported by Goma (except FRONT), and with the provided sub-package PARPACK can 
be used for parallel eigensolves. The interface to ARPACK is contained in LOCA, which enables 
continuation to be combined with stability analysis at specified step intervals

The sections of this chapter present an explanation of the input specifications (4.1 and 2.2), an 
example to illustrate the use of LSA and to highlight the output and it's meaning (4.4), and a 
presentation of several considerations relevant to the use of LSA (4.5). Where necessary, 
notations will be made in the remainder of this chapter about the applicability of inputs, features, 
etc. to either one or both eigensolvers.

Required Specifications in the Goma Input File
===================================================

This section describes how to carry out linear stability analysis with Goma. The input information 
that Goma requires is supplied in the usual Goma input file under the Eigensolver Specifications 
section. To activate the stability software, the Linear Stability keyword in the Solver 
Specifications section must be set equal to yes, that is:

::

    ----------------------------------------------------------
    Solver Specifications
    ----------------------------------------------------------
    Linear Stability = yes

A sample of the information/parameters that must be supplied in the Eigensolver Specifications 
section is listed below. The eigensolver applicability of each input is indicated here by <A> for 
ARPACK only, <E> for eggroll only, or <AE> for both.

::

    ----------------------------------------------------------
    Eigensolver Specifications
    ----------------------------------------------------------
    Eigen Algorithm = cayley <A>
    Eigen Number of modes = 10 <AE>
    Eigen Record modes = 5 <AE>
    Eigen Size of Krylov subspace = 30 <AE>
    Eigen Maximum Iterations = 2 <E>
    Eigen Number of Filter Steps = 1 <E>
    Eigen Recycle = no <E>
    Eigen Tolerance = 1.0e-06 <E>
    Eigen Matrix Output= no <AE>
    Eigen Initial Vector Weight = 0.1 <E>
    Eigen Initial Shifts = -50.0 -51.0 -52.0 -54.0 <E>
    Eigen Cayley Sigma= 100.0 <A>
    Eigen Cayley Mu= 1000.0 <A>
    Eigen Relative tolerance= 1.0e-6 <A>
    Eigen Linear Solver tolerance= 1.0e-6 <A>
    Eigenvalue output frequency= 1 <A>
    Eigenvector output frequency= 1 <A>
    Eigenvector output file= eigen.exoII <A>
    Eigen Wave Numbers= 0.0 0.1 0.2 1.0 <AE>

Eigensolver Specifications
==========================

The ability to solve for the stability of a base flow is a very powerful tool. Often, the important 
characteristics of a flow can be summarized in the answer to the question "is the flow stable?". 
Although the following cards are in active use at the time of this writing, sweeping changes are 
coming to the eigensolver sections of Goma. In particular, the old code (called "eggroll") is being 
replaced with newer methods (in the ARPACK library), as well as being coupled to the 
continuation and tracking algorithms (in the LOCA library).

Input specifications for this section of input records is discussed in a separate, comprehensive 
manual (Gates, et. al., 2000); an update to this manual will be completed during the fall of 2002 
(Labreche, et. al., 2002). Either of these manuals contains a thorough discussion of how to 
successfully compute the stability and interesting modes of an underlying base flow.

.. toctree::
   :maxdepth: 1

   linear_stability/eigen_algorithm
   linear_stability/eigen_number_of_modes
   linear_stability/eigen_record_modes
   linear_stability/eigen_size_of_krylov_subspace
   linear_stability/eigen_maximum_iterations
   linear_stability/eigen_number_of_filter_steps
   linear_stability/eigen_recycle
   linear_stability/eigen_tolerance
   linear_stability/eigen_matrix_output
   linear_stability/eigen_initial_vector_weight
   linear_stability/eigen_initial_shifts
   linear_stability/eigen_wave_numbers
   linear_stability/eigen_cayley_sigma
   linear_stability/eigen_cayley_mu
   linear_stability/eigen_relative_tolerance
   linear_stability/eigen_linear_solver_tolerance
   linear_stability/eigenvalue_output_frequency
   linear_stability/eigenvector_output_frequency
   linear_stability/eigenvector_output_file

3D of 2D Stability Analysis
================================

Theory
-------------

For some 2D flow problems, it may be possible that the steady state flow is stable to all 
disturbances in either flow direction but unstable to disturbances in the third (transverse) 
direction. This situation could be analyzed by a full 3D simulation with stability analysis, but this 
would have considerably larger computational expense. An alternate method is to use normal 
mode analysis (e.g. Coyle 1984, Christodoulou 1988, Gates 1998) to model the transverse 
disturbances. This approach considers such disturbances as the sum of a series of Fourier normal 
modes of different wavelengths, or different wavenumbers - in this case, the actual transverse 
stability is that of the "most dangerous" normal mode -- i.e. the one which yields the largest 
leading eigenvalue (real part if complex). These authors have developed algorithms to include 
Fourier factors in the base 2D equations for residual and Jacobian assembly which account for the 
effects of these normal modes, such that the 2D finite element mesh can still be used and only a 
third velocity component needs to be added to the eigensystem.

This is the general approach taken in Goma for 3D of 2D stability analysis, which makes this 
method accessible to a wide range of physical problems. Basically, when 3D of 2D stability is 
requested for a given list of wavenumbers, the steady state problem is first solved as usual, then 
the augmented eigensystem is assembled using a two-pass approach for both the Jacobian and the 
mass matrix. All of the necessary factors are included with the basis and weight functions, to 
minimize intrusion on the assembly equations.

Usage
------------

The 3D of 2D stability algorithm is accessed by the user similarly to basic stability analysis, and 
requires only a few changes to an input file which is already set up for LSA (using either eggroll 
or ARPACK). One requirement of this method is that both momentum equations of the 2D 
problem must be included and must use Q2 interpolation for basis and weight functions; a third 
momentum equation will then be added as shown in Example 1.4.2.

For a file which is set up for LSA, the following changes must be made:

• Add the "momentum3" equation, using "Q2_LSA" interpolation type.
• Update (add 1 to) number of equations.
• Add BC's for third velocity (W)
• Change Coordinate System to "PROJECTED_CARTESIAN"
• Change Linear Stability choice to "3D"
• Add the Eigen Wave Numbers card with a list of (1 or more) wavenumbers.

When used with either eigensolver, a list of converged eigenvalues will be generated on screen for 
each requested wave number, but print formats will differ. If the Eigen Record Modes card is set 
to n (1 or more), then the eigenvectors for the first n modes will be written to separate ExodusII 
files. In both cases, the file names are appended with mode and wave numbers, but the base names 
will differ. Also, eigenvector output files generated with ARPACK will have zero-based mode 
numbers and the eigenvalue itself is used as the time stamp. For continuation runs with ARPACK, 
each of these files will contain eigenvectors for the same mode/wavenumber combination at each 
continuation step. If the run was done in parallel, it will be necessary to run fix for any file you 
wish to view.

Examples
=============

Stability of the lid driven cavity problem at Re = 1
-----------------------------------------------------------

The relevant sections of the input file to carry out linear stability analysis are as follows:

::

    ----------------------------------------------------------
    Solver Specifications
    ----------------------------------------------------------
    .
    .
    Linear Stability = yes
    ------------------------------------------------------------
    Eigensolver Specifications
    ------------------------------------------------------------
    Eigen Number of modes = 10
    Eigen Record modes = 5
    Eigen Size of Krylov subspace = 30
    Eigen Maximum Iterations = 2
    Eigen Number of Filter Steps = 1
    Eigen Recycle = no
    Eigen Tolerance = 1.0e-06
    Eigen Initial Vector Weight = 0.1
    Eigen Initial Shifts = -50.0 -51.0 -52.0 -54.0

To run the stability software, the Linear Stability keyword in the Solver 
Specifications section was set to yes. Given that this problem is small (dimension is 1182), 
the dimension of the Krylov subspace was taken to be 30. The number of modes wanted was 
chosen to be 10 with the five most dangerous modes being written to file. The shifts were chosen 
to be close to the leading eigenvalue found below. This was found by trial and error (refer to 
discussion in Section 4.5). The initial vector weight was taken to be 0.1. In this problem, changing 
the initial vector weight did not help speed up the eigenvalue extraction.

Goma converges to a solution for the cavity problem; the linear stability analysis is performed on 
the converged solution. The Goma output for both analysis steps is:

::

    R e s i d u a l            C o r r e c t i o n
     ToD    itn    L_oo    L_1     L_2     L_oo    L_1     L_2   lis asm/slv (sec)
    -------- --- ------- ------- ------- ------- ------- ------- --- ---------------
    11:56:10 [0] 2.6e-12 2.2e-11 5.3e-12 1.9e-10 1.3e-08 1.2e-09   1 5.0e-02/6.0e-02
    scaled solution norms  3.386794e+01  4.267119e-01  2.028370e+00

    WARNING: Linear Stability Analysis is not (and cannot be) compatible
     with all user-supplied boundary conditions (i.e., BC cards with
     the term USER in them). If you are using some, and you want LSA,
     then you must modify your boundary conditions accordingly.

    Assembling J...
    Assembling B...
     Initializing variables and allocating space...
     Arnoldi Eigenvalue Extractor
     Allocated work vectors
     Allocated work martices
     Pass 1 New Shift = -5.000000e+01 +0.000000e+00 i
     Pass 2
     De-allocated work storage
     ARN-IT done.
    -------------------------------------------------------------------------------
     Eigensolver required 85 iterations.
     Found 8 converged eigenvalues.
     Leading Eigenvalue = -5.235398e+01-0.000000e+00 i RES = 4.163666e-47
     Real            Imag            RES
     -5.235398e+01  -0.000000e+00 i  4.163666e-47
     -9.219072e+01  +3.614838e-01 i  7.087330e-29
     -9.219072e+01  -3.614838e-01 i  7.087330e-29
     -1.283342e+02  -0.000000e+00 i  9.032547e-20
     -1.544776e+02  -0.000000e+00 i  5.007987e-16
     -1.675191e+02  -0.000000e+00 i  1.794834e-13
     -1.899878e+02  -8.560254e-01 i  2.109350e-10
     -1.899878e+02  +8.560254e-01 i  2.109350e-10
     push_mode                       =           5
     Writing modes to file ...
     Mode     0 ... recorded.
     Mode     1 ... recorded.
     Mode     2 ... recorded.
     Mode     3 ... recorded.
     Mode     4 ... recorded.
    Deallocating memory ... done.
    -done

    Proc 0 runtime:         0.02 Minutes.

The Goma steady state solution for this problem was written to the output exodusII file named 
out.exoII. The names of the files containing the five eigenvectors requested for writing (by the 
Eigen Record modes input) are:

LSA_1_of_5_out.exoII

LSA_2_of_5_out.exoII

LSA_3_of_5_out.exoII

LSA_4_of_5_out.exoII

LSA_5_of_5_out.exoII

Each file is an exodusII file that can be read and displayed by BLOT. These eigenvectors 
correspond to the first five eigenvalues. Note that this means that LSA_2_of_5_out.exoII
corresponds to -9.219072e+01 +3.614838e-01 i whereas 
LSA_3_of_5_out.exoII corresponds to -9.219072e+01 -3.614838e-01 i in this 
example. A velocity vector plot of the leading mode is shown in Figure 4.1.

3D of 2D stability of the lid driven cavity problem at Re=1
------------------------------------------------------------------

To access the 3D of 2D stability algorithm for this problem, the previous input file is modified 
as indicated in Section 2.3.2. The affected cards are:

::

    Linear Stability= 3D
    Eigen Initial Shifts= -40.0 -51.0 -52.0 -54.0
    Eigen Wave Numbers= 0.0 1.0 2.0
    BC= W NS 3 0.0
    BC= W NS 1 0.0
    BC= W NS 2 0.0
    BC= W NS 4 0.0
    Coordinate System= PROJECTED_CARTESIAN
    Number of EQ = 4
    EQ= momentum3 Q2_LSA U3 Q2_LSA 0. 1. 1. 1. 0. 0.

Note that the new momentum3 EQ card must be inserted just below the existing momentum2 
card -- the three cards must remain together. Also, the first Eigen Initial Shift was changed to 
accommodate the higher wavenumbers.

The output from Goma shows leading eigenvalue (real part) increases with wavenumber over this range, 
indicating that some transverse perturbations are less stable than those in the other two directions 
(although not unstable in this case). This run generated fifteen eigenvector output files (one for 
each requested mode/wavenumber combination) with names appended accordingly. These names 
are of the type "LSA_5_of_5_wn=2_out.exoII"; those generated with ARPACK contain this 
information in a slightly different format. One check on this algorithm is that when the 
wavenumber is zero (or infinite wavelength), the standard 2D stability results should be recovered 
-- comparison of the first set of eigenvalues with those found in Example 4.4.1 reveals that this is 
the case.

Stability of the lid driven cavity with continuation in Re
----------------------------------------------------------------

To perform continuation with stepwise stability analysis, it is necessary to switch from eggroll to 
ARPACK; this is achieved when LOCA continuation is specified. The input file is modified by 
removing the Hunting Specifications section (including its header) and adding or changing the 
following cards:

::

    (in Continuation Specifications:)
    Continuation= loca
    Continuation Type= BC
    Boundary condition ID= 4
    Boundary condition data float tag= 0
    Material id= 1
    Material property tag= 1700
    Material property tag subindex= 0
    Initial parameter value= 1.0
    Final parameter value= 5.0
    delta_s= 2.0
    Maximum number of path steps= 3
    Minimum path step= 10.0
    Maximum path step= 50.0
    Continuation Printing Frequency= 1
    LOCA method = zero

    (in Eigensolver Specifications:)
    Eigen Algorithm= si
    Eigen Cayley Sigma= -50.0

These cards set up a continuation run in the lid speed (here, taken to be Re) from 1 to 5 in constant 
steps of 2 by changing the float on BC #4 (fifth from the top) with a call to ARPACK (shift-andinvert method with shift of -50) after each step is completed. The Goma output for this run shows 
eigenvalue calculations at each continuation step, with eigenvalue and eigenvector frequencies controlled 
by the respective output frequency cards.

While ARPACK supports an option to specify a base name for eigenvector output files, one was 
not selected here and the default "LSA.exoII" was applied. Five eigenvector files were created with 
names of the type "LSA_mode0.exoII" Note that these modes are zero-based whereas they are 
one-based with eggroll. Each of these files contain three "time steps" which correspond to the 
continuation steps taken, and the time stamp for each step is the corresponding eigenvalue. Note 
also that the computed eigenvalues for the first step agree with those found in Example 4.4.1 
(with eggroll).

User Guidance
==================

Linear stability analysis of the solution of a previously unstudied process can be quite difficult 
without some fortuitous choices in the LSA setup. Experience will erase the uncertainty but the 
user requires a jump off point. To this end, the following guidance is provided; included are some 
do's and don'ts in problem setup.

User Boundary Conditions
--------------------------------

In general, a LSA should not be performed if user-specified boundary conditions have been used 
in the Goma analysis, since boundary conditions are the main source of difficulties in doing an 
LSA. To solve for the requested eigenvalues, a matrix is constructed from the time derivatives of 
the solution variables. The boundary conditions from the associated transient problem are 
considered, and only those boundary conditions that contain time derivatives are inserted into this 
matrix. This can lead to unexpected complications: one would not normally consider the timedependent boundary conditions when solving a steady-state problem. User-supplied boundary 
conditions should not be used, unless the user is aware of, and can implement, the appropriate 
time derivatives in the user function necessary to construct the correct matrix.

Leading Eigenvalue
-------------------------

The LSA reports a Leading Eigenvalue and a list of values up to the Number of Modes (Section 
2.2.1) selected in the analysis. The Leading Eigenvalue is simply the eigenvalue with the largest 
real part from the list. If the shift has been chosen well (i.e., near the actual leading eigenvalue), 
then this will be the eigenvalue with the largest real part. However, if the shift was chosen poorly, 
or a phantom eigenvalue is present (the eigensolver may report an eigenvalue converged when it 
is not an actual eigenvalue), then the leading Eigenvalue may be misleading. In this case, repeated 
runs should be performed and the list consulted (ignoring the Leading Eigenvalue reported). This 
is true to a greater extent with eggroll than ARPACK, as the LOCA routines that check for 
convergence are more effective in discarding spurious modes, thereby detecting fewer phantom 
eigenvalues for similar conditions. Thus, the choices of ARPACK shift parameters appears to have 
a smaller effect on the reliability of the results. But beware that ARPACK may generate an error if 
it detects an eigenvalue to the right of (larger than) the primary shift parameter sigma; when this 
happens, just increase sigma and try again. Another way to decrease the effect of spurious modes 
(for either method) is to increase the Krylov subspace size.

Selecting Initial Shifts (eggroll)
-----------------------------------------

The eggroll eigensolver is very sensitive to the user-supplied shifts (Section 2.2.9). The guidance 
is to pick a set of initial shifts which has at least one shift to the right (in the positive real 
direction) of the leading eigenvalue. In the example in Section 4.4, the initial shifts selected were 
-50, -51, -52, -54 while the leading eigenvalue was determined to be -52.35. With time this 
portion of the LSA module will become more robust, but for the present a trial and error approach 
is necessary for new models/processes. Following is some assistance to carry out this trial and 
error procedure.

The first shifts can be selected to span several orders of magnitude, or successive analyses can be 
used to span several orders of magnitude, e.g., -1 to -1000. (Recall, we are looking for stable 
solutions, so the eigenvalues will be negative.) By expanding this range, reducing the range, or 
refining within a selected range, the user should look for consistency in the leading eigenvalue, 
i.e., the leading eigenvalue should emerge from all the stable modes as the largest real value. Then 
a final set of shifts can be chosen which is in keeping with the "... at least one shift to the right of 
the leading eigenvalue" approach.

Alternatively, the user could apply "analytic intuition." The eigenvalues being sought imply a 
time constant, τi, for the problem of interest, which is simply

.. math::

    \lambda_i = \frac{1}{\tau_i}                                                    (4-1)

The decay (or growth in unstable solutions) is given by :math:`e^{\lambda_i t}`. If a user has some insight into the 
"time constant" for the process being modeled, i.e., seconds or milliseconds, for example, an 
estimate of the leading eigenvalue can be made to guide the trial and error search referred to in the 
above paragraph.

Selecting Initial Shifts (ARPACK)
----------------------------------------

The two Cayley shift parameters, σ and μ, are specified using the Eigen Cayley Sigma and 
Eigen Cayley Mu input cards, respectively. With the LOCA drivers for ARPACK there are 3 
transformations implemented. The Shift-and-Invert behaves the same as the one in eggroll and 
one can follow the same strategy in choosing the shift parameter (in fact it was added to LOCA 
expressly for eggroll users), and in this case μ does not apply. The two Cayley transformation 
options available through LOCA for driving ARPACK (A and B as discussed below) have 
different strengths depending on the type of instability that is expected. For reliability in 
determining the stability of a solution, the LOCA developers recommend Cayley version B. 
Method A is invoked by choosing σ < μ, and Method B is invoked by choosing μ < σ.

The Cayley version A has λ < σ < μ, where λ is the expected leading eigenvalue, and is good for 
finding real eigenvalues or those with relatively small imaginary parts. Since eigenvalues are 
dimensional quantities of units inverse time, it is not possible to pick good default values for the 
transformation parameters. The paper by Lehoucq and Salinger (2001) presents some guidance on 
choosing the Cayley parameters and discusses the trade-offs. If one wants to converge the first n 
eigenvalues and expects them to have real parts in the range -10 < λ < 0 and imgainary parts less 
then 10i, a reasonable choice would be σ=5 and μ = 25 n (Remember that σ must be to the right of 
all eigenvalues. Choosing μ to be 5 times further away from the eigenvalues than σ works well. 
Typical choices for the Krylov subspace size would be 20-50.

The Cayley version B has μ < σ and usually has μ = -∞. This transformation is very robust for 
finding Hopf instabilities where there is a relatively large imaginary component to the eigenvalue. 
There are some details of this transformation in the Sandia report by Burroughs et al.(2001). If 
one expects an instability with pure imaginary λ = iω (i.e. with period 2π/ω) the a good choice 
would be σ = ω and μ = -∞. This transformation corresponds directly to trapezoid rule time 
integration, with σ=2/Δt, so one can use intuition on choosing time steps to choose the Cayley 
parameters. Typical values would be σ=100, μ=-100. With this transformation, we often need 
larger Krylov subspaces than the other transformations, often 50-200; yet when using an iterative 
solver this expense is in part mitigated by the fact that the matrices are diagonally dominant and 
solve very quickly. Also, this method is by far the most robust in locating all eigenvalues with 
positive real part, irrespective of the imaginary part.
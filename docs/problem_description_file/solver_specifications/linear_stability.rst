********************
Linear Stability
********************

::

	Linear Stability = {char_list}

-----------------------
Description / Usage
-----------------------

This optional card indicates whether or not linear stability analysis should be
performed, as well as what kind.

The valid options for {char_list} are:

no
    Do not perform any kind of linear stability analysis.
yes
    Perform regular linear stability analysis. If your problem was 2D, then 2D
    analysis is performed. If your problem was 3D, then 3D analysis is
    performed.
inline
    Same as yes, perform regular linear stability analysis.
3D
    Subject the 2D flow to 3D linear stability analysis by normal mode
    expansion for the modes specified with the *Eigen Wave Numbers* card.
file
    Set up the problem as in **yes** or **inline,** but output the matrices
    involved instead of determining stability.
3Dfile
    Set up the problem as in **3D,** but output the matrices involved instead
    of determining stability.

The default value is **no.**

------------
Examples
------------

Here is a sample card:
::

	Linear Stability = yes

-------------------------
Technical Discussion
-------------------------

When linear stability analysis is performed, a steady-state solution is first acquired, and
then the eigenvalue/eigenvector spectrum is computed subject to the choices made in
the *Eigensolver Specifications* section. In the case of **file** or **3Dfile,** the steady-state
solution is acquired and then the matrices that would have been used to compute the
spectrum are exported to file and no spectrum is actually computed. Refer to the
Advanced Capabilities (Gates, et. al., 2001) document for a more thorough description.

The name of the output files when **file** is specified are:

* LSA_mass_coo.out         for the mass matrix, B or M,
* LSA_jac_coo.out          for the jacobian matrix, J,
* LSA_vars.out             for variable names associated with unknowns.

When **3Dfile** is specified, the names are:

* LSA_mass_coo-<f>.out,    for the mass matrix, B or M,
* LSA_jac_coo-<f>.out,     for the jacobian matrix, J,
* LSA_vars.out,            for variable names associated with unknowns.

where <f> is the value of the requested normal mode (see the *Eigen Wave Numbers*
card). The *Eigen Matrix Output* card must be set to **yes** in order to create and write
these files.

When computing the 3D stability of a base 2D flow, other modifications need to be
made (see the 3D stability of 2D flow memo).

See the Advanced Capabilities document (Gates, et. al., 2001), or it’s replacement
(Labreche, et. al., 2002).

--------------
References
--------------

SAND2000-2465: Advanced Capabilities in Goma 3.0 - Augmenting Conditions,
Automatic Continuation, and Linear Stability Analysis, I. D. Gates, D. A. Labreche and
M. M. Hopkins (January 2001)

SAND2002-xxxx: Advanced Capabilities in Goma 4.0 - Augmenting Conditions,
Automatic Continuation, and Linear Stability Analysis, Labreche, D. A., Wilkes, E. D.,
Hopkins, M. M. and Sun, A. C., (in preparation)

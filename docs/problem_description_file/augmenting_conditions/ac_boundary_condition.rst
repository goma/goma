***********************
AC (Boundary Condition)
***********************

::

    AC = BC <bc_id> <integer> <float_list>

-----------------------
Description / Usage
-----------------------

This type of augmenting condition allows the user to connect an additional constraint 
equation (user-supplied) to a specific boundary condition float parameter which is 
allowed to vary as an additional degree of freedom during the solution process. The 
constraint equation is added to the function user_aug_cond_residuals in the file 
user_ac.c. A discussion of this function specification and examples is supplied in later 
sections of this manual.

A description of the card syntax follows:

BC
    A mandatory string indicating that the augmenting 
    condition is attached to a boundary condition parameter.

<bc_id>
    An integer parameter identifying the boundary condition 
    index whose parameter is going to be used as the additional 
    unknown. The first index is zero, starting from the first read 
    boundary condition in the input file, and proceeding 
    sequentially upwards. The Technical Discussion describes a 
    method for automatically determining <bc_id>.

<integer>
    A parameter identifying the index of the float parameter that 
    is to be varied on the boundary condition specification 
    identified by <bc_id>. The leftmost float value is assigned 
    index zero and the index increases sequentially left to right.

<float_list>
    A list of float parameters that can be used in 
    user_aug_cond_residuals to evaluate the augmenting 
    conditions. They are stored in sequence in the array 
    augc[i].DataFlt.

------------
Examples
------------

For the following set of boundary condition cards:

::

    Number of BC = -1
    BC = U NS 10 0.0
    BC = V NS 10 0.0
    BC = W NS 10 0.0 1.0

an example augmenting condition card of type BC is:

::

    AC = BC 2 0 1.0 0.25

The augmenting condition card attaches the (unspecified, in this example) constraint to 
the first float parameter on the BC = W NS 10 card. Note that there are two float 
parameters on this card - the first is the specified z-velocity component on node set 10 
and the second parameter is required when a Dirichlet condition is attached to an 
augmenting equation. Ordinarily, Dirichlet conditions are applied by direct 
substitution, but when the second parameter is present, the condition is included as a 
residual equation along with all the other residual equations. Naturally, it is this latter 
form that should be used when an augmenting condition is attached to the boundary 
condition. The two float parameters supplied on the card can be used in 
user_aug_cond_residuals as the variables augc[0].DataFlt[0] = 1.0, 
augc[0].DataFlt[1] = 0.25.

-------------------------
Technical Discussion
-------------------------

- The function user_aug_cond_residuals is passed an integer index iAC. This 
  is the index of each AC card in order of its appearance in the input deck (with the 
  zero being assigned to the first AC card that appears). The arbitrary float 
  parameters supplied on the card in <float_list> can be accessed in this function in 
  the array augc[iAC].DataFlt. The variable DataFlt[0] corresponds to the 
  first float parameter in <float_list> and so on.

- Often the augmenting condition constraint written in function 
  user_aug_cond_residuals must make use of the values of nodal unknowns. 
  These can be determined using the function Index_Solution as follows:

  The boundary condition index of any boundary condition can be 
  found with the following feature. The NS or SS designator string on 
  the boundary condition in question should be changed 
  (temporarily) to NC or SC, respectively. Then Goma should be run 
  with the arguments as follows:

  ::

      goma -i input_file -a -bc_list

  This invokes a code fork that does not run the problem but analyses 
  the input deck and prints out the index of the boundary conditions 
  flagged by the previous procedure. The values of these indices can 
  be used directly in the augmenting condition cards for <bc_id>

- The augmenting condition constraints are additional residual equations solved 
  simultaneously with the other residual equations associated with the nodal degrees 
  of freedom. The unknown degrees of freedom associated with these equations are 
  the boundary condition float values identified on each augmenting condition card. 
  The new equations also introduce a different structure to the matrix being solved. 
  The formerly sparse, banded Jacobian matrix has now been augmented by rows 
  and columns on its periphery that are potentially populated. This requires a 
  bordered matrix algorithm for the solution of this type of matrix. A consequence 
  of this algorithm is that the residual and update norms of the augmenting 
  conditions are computed separately. The user will see these norms appear as a 
  second row of output at each and every iteration. If these are missing, the 
  augmented system is not being solved.

- Note that the value of the augmenting condition parameter values can be displayed 
  after each iteration by setting Debug = 1 in the input deck.


*********
**Debug**
*********

::

	Debug = <integer>

-----------------------
**Description / Usage**
-----------------------

This optional card specifies the level of information output to files stdout and
stderr. The permissible values for <integer> are **-3** through **4**, depending on the level
of informational (debugging) output desired; higher values of the output level will
produce more diagnostic information output on the stdout and stderr output
channels. The default level is 0. Specific results produced for each level are
summarized below. The user should exercise caution in using values other than the
default for problems with large numbers of unknowns as the volume of information
increases very quickly.

.. tabularcolumns:: |l|L|

==============  ===============================================================
Level           Results Output
==============  ===============================================================
0               No output (default).
1               Logs activity as the code does problem setup, including setting
                addresses and sizes, boundary-condition (BC) conflictresolution
                information, and identification of the rotation
                conditions at every node with a boundary flag. Prints out surface
                boundary integral setup information. Lists matrix and solver
                information for each solution step.
2               Prints same information as level 1, plus provides a summary of
                BC type information for each BC and logs the beginning and
                end of matrix fill operations.
3               Prints same information as level 2, but also prints a list of
                variables/unknowns at each node.
4               Prints same information as for level 3.
-1              Logs activity as the code does problem setup and prints out
                surface boundary integral setup information as is done for mode
                1. Triggers a comparison of the analytical Jacobian and the
                numerical Jacobian in un-scaled form, which can be used to
                check the compatibility of the analytical residual equations and
                Jacobian. Prints results only if the analytical and numerical
                Jacobian are different. Does not solve any equations; terminates
                after Jacobian print out.
-2              Same initial information as for level -1. Triggers a comparison
                of the analytical Jacobian and the numerical Jacobian scaled by
                the sum of each row of the analytical Jacobian (this helps
                suppress small errors in large Jacobian entries). Prints results
                only if the analytical and numerical Jacobian are different.
-3              Similar to level -2 except each row is scaled by the diagonal
                value which is usually the largest. Prints results only if the
                analytical and numerical Jacobian are different.
==============  ===============================================================

------------
**Examples**
------------

Following is a sample card:
::

	Debug = -2

-------------------------
**Technical Discussion**
-------------------------

For options -1, -2, -3, viz. numerical Jacobian checking, the user must take care when
interpreting the cited differences in the numerical and analytical Jacobian. The comparison is made by perturbing each variable and comparing the numerical Jacobian
computed between the perturbed and unperturbed states to the analytical Jacobians at
the two states. A difference is deemed significant if the numerical Jacobian falls outside
the band between the two analytical values with an additional allowance for roundoff
error. It is the roundoff error in the residual that is the most difficult for the Jacobian
checker to estimate. This is particularly true for problems with zero initial conditions
since it is impossible to determine the scale of a velocity, for example, if all the values
of velocity are zero. For this reason, it is often better to use a nonzero initial condition
or a scaled problem (with values order unity) when using the Jacobian checker.

Currently, there are two parameters output by the Jacobian checker that can help the
user decide on the significance of the entry. The first is the relative change in the
analytical residual. This quantity, labeled d\ :sub:`aj`, is the percentage of the acceptance band
that comes from changes in the analytical Jacobian from the unperturbed to perturbed
states. For a non-linear dependency, the difference between the analytical Jacobians
will be significant and it is reasonable to expect that the numerical Jacobian should fall
within the band. If the analytical Jacobian is nearly constant over the perturbation, the
accuracy of the check becomes increasingly dependent on knowing the roundoff error
in the residual. So, as d\ :sub:`aj` gets closer to unity, the user can have more confidence that
the entry is significant.

The second parameter is a confidence measure that is the deviation between the
numerical jacobian and analytical values divided by the expected value of the deviation
based on roundoff error. Since the roundoff error is only known approximately, this
value, called *conf*, is only a qualitative measure of the confidence. A *conf* value of 100
means that the deviation between the numerical jacobian and the analytical values is
100 times larger than the expected deviation based on roundoff error.

Here is a sample of output from a convective heat transfer problem, using the -2 option

:: 

    Eqdof=92 T_0 n=31 Vardof=95 T_0 n=32 x=0
    dx=0.0001 aj=-0.008188 nj=-0.008126 aj_1=-0.008188 d_aj=0
    conf=1.889e+06
    
      >>> QCONV on SSID=1

This entry can be read as follows: The sensitivity of global equation number 92, which
happens to be the T_0 energy equation at node 31, with respect to the temperature
variable at node 32 (variable global degree of freedom number 95) has an analytical
Jacobian of -0.008188 at the unperturbed state and a computed numerical Jacobian of
0.008126. The analytical jacobian at the perturbed state is -0.008188. For this problem
the change in the analytical Jacobian is zero between the unperturbed and perturbed
states, so d\ :sub:`aj` is zero. But even though the difference is small between the analytical and
numerical values, it is huge relative to the expected roundoff error, with the deviation
being 1.889e+6 times the deviation attributable to roundoff error.

For each node where a deviation is found, the side boundary conditions applied at the
node are printed, as shown above. If one of these boundary conditions are applied to the
equation that shows an error and have the same dependency that is showing the error,
this boundary condition is flagged as shown for the QCONV boundary condition
above.

Before the user/developer concludes that there is a discrepancy in the analytical
Jacobian, a few things should be tried:

* Giving the problem a nonzero initial guess, either by reading in a STEADY state
  solution, if one exists, or on transient problems using the “one” option on the
  Initial Guess card. Sometimes this will make many differences disappear.

* Checking whether the nodes cited in the difference outputs are boundary nodes.
  Specifically, if they are boundary nodes on which Dirichlet boundary conditions
  are specified, artificial errors can occur.

* Also, if you are in doubt that there are not reported errors, put one in by a 10
  percent perturbation to the residual. The Jacobian checker should hit on those
  errors and report them to you.

* Check the settings in *mm_numjac.h*.


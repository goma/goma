Augmenting Conditions Specifications
########################################

The augmenting condition capability in Goma enables the addition of supplemental constraints
and auxiliary conditions to the problem via standardized and user-defined function(s) (in files
mm_augc_util.c and user_ac.c respectively). These supplemental constraints (equations) are used
to define relationships between unknowns, boundary condition data, material properties, and
virtually any other extracted quantity from Goma post processing routines. Examples include:

**Goma unknowns:**
  Specify a relationship between two or more unknowns, e.g. the distance between 
  two nodes remains fixed at a specified value.

**Material properties:**
  Solve for a material property value such that another condition is met, e.g., find the 
  liquid viscosity at which the velocity gradient at a wall becomes zero, or to maintain a relationship between the temperature, pressure, and density.

**Boundary condition floats:**
  Constrain the problem so that a specified relationship always exists between 
  boundary condition data floats (which may not necessarily be Dirichlet conditions), i.e., the parameters of a boundary condition specification, or between parameters of different boundary conditions.

**Postprocessing constraints:**
  Postprocessing of solutions can define auxiliary variables that allow feedback of 
  integrated results into the solution as constraints. Examples include integral constraints on surface forces, heat flux through a wall, volume flux and species flux.

**Other constraints:**
  Examples include constraints on volume, mass and component mass.

**Mixed constraints:**
  Augmenting conditions can invoke any combination of the above types in a single 
  auxiliary constraint. Thus, a relationship can be specified between any of the 
  Goma unknowns, material properties, boundary condition floats, and postprocessed variables.

The Goma augmenting condition capability handles all of the above types. Providing the variable 
of interest is passed to the augmenting condition user-definition routine (described below) via a 
data structure or argument, the augmenting condition can be applied. The number of augmenting 
conditions is unlimited but must be set by the user. An important consideration in the addition of 
augmenting conditions however is increased computational cost (time and storage) as the number 
of augmenting conditions increases. There are no limits on the density (number of non-zeroes) of 
the augmenting condition residual equation or the level of interaction between the unknowns.

Two pieces of information must be specified to add an augmenting condition to a Goma analysis: 
identification of the AC and its parameters (Section 2.3) and the addition of new unknowns to the 
Goma system of equations (Section 2.4). This latter step consists of creating a function which 
specifies the relationship between the variables of interest. The extra unknowns must be equal in 
number to the augmenting conditions and be one of four types: boundary condition data floats, 
material properties, volume constraints or flux constraints. Volume and flux augmenting condition 
standardized relationships have already been defined in function std_aug_cond (file 
mm_augc_util.c); the user must specify the functional relationship for boundary conditions and 
material properties in the user_aug_cond_residuals function (file user_ac.c).


Augmented System Solution
*************************

Adding extra equations to the original Goma system of equations yields an augmented set of 
equations. In Goma, the augmented system of equations is solved by a bordering algorithm (Chan 
and Resasco, 1986). This method is essentially a block elimination of the augmented equations.

The original Goma system has the advantage that it is sparse, usually with small bandwidth. This 
leads to computational benefits (speed and storage) that can be realized when using direct frontal 
and iterative solvers. The bordered algorithm is suitable for the augmented system of equations 
because the decomposition of the original Goma set of equations is already required; one 
additional, but much smaller, decomposition for the augmented equations adds little calculational 
cost compared to the original system decomposition. If there are N augmenting conditions, the 
bordered algorithm solves the augmented system of equations with one decomposition and N+1 
back-substitutions of the Goma original system, plus one NxN system solution. One further 
advantage of the bordered algorithm is that the augmenting equations can have much wider 
bandwidth than the original system without substantial growth in the computer time and memory 
for decomposition.

In the case where the original system is singular but the augmented system is not, the bordered 
algorithm will fail to solve the system. An instance in which this occurs is the tracking of a 
turning point in a parameter. Here, the original Goma set of equations is singular while the 
augmented system is not.

Required Specifications in the Goma Input File
**********************************************

A new section must be added to the Goma input file to identify the new unknowns and other 
needed data. This section is required only when augmenting conditions are present in the 
numerical model and has the following form:

.. code-block:: none

    -------------------------------------
    Augmenting Condition Specifications
    -------------------------------------
    Number of augmenting conditions  = -1
    AC = BC 29 1
    AC = MT 4 1300
    END OF AC



Augmenting Condition Specifications
***********************************

The following cards are used to specify augmenting conditions in Goma.

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/augmenting_conditions_initial_guess.rst

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/number_of_augmenting_conditions.rst

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/ac_boundary_condition.rst

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/ac_material_property.rst

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/ac_flux_condition.rst

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/ac_volume_constraint.rst

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/ac_overlapping_grid_boundary_conditions.rst

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/ac_phase_velocity_level_set.rst

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/ac_periodic_boundary_condition.rst

---------------------------------------------------------------------------------------

.. include:: augmenting_conditions/end_of_ac.rst

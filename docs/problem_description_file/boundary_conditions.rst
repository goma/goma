Boundary Condition Specifications
#####################################

The broad range of mechanics capabilities that has been built into *Goma* necessitates an equally
broad range of boundary conditions (BCs) to provide all boundary condition information that the
differential equations specified in the *Problem Description* section will require for a well-posed
system. The BCs for *Goma* have been categorized according to the differential equation set to
which they apply. First are listed those boundary conditions which can be applied to any equation
followed by BCs for mesh, real solid, fluid momentum, energy, mass, continuity, porous, stress,
gradient, shear rate, fill and potential equations. Each boundary condition (BC) card follows a
general syntax as follows:

::

	BC = <bc_name> <bc_type> <bc_id> {integer_list}/{float_list}

The <bc_name> identifies the desired control of the physics/mechanics at the boundary as
identified by the <bc_type> and its associated <bc_id>. The <bc_type> is either nodeset, NS
(NODEBC or POINBC in EXODUS II) or sideset, SS (ELEMBC in EXODUS II) depending on
the <bc_name> and can be located in the problem domain by means of its flag or <bc_id> number
(set in EXODUS II). The {integer_list} and/or {float_list} specify parameters of the boundary
condition. Within each equation category are Dirichlet nodeset boundary conditions (i.e. T, U, V,
W, DX, DY, DZ, Y, S11, S12, S13, S22, S23, S33, G11, G12, G13, G21, G22, G23, G31, G32,
G33) that can be handled (i.e., processed) in two ways in *Goma*. The first way is application of the
BC as a “hard-set” on the primitive variable, and the second as a residual equation; differences in
these methods are discussed below. The cards belonging to this category have the following
general syntax:

::

	BC = <bc_name> <bc_type> <bc_id> <float1> <float2>

where <float2> flags whether a hard-set or residual equation is to be used.

Prior to introducing individual boundary conditions and their parameters, some general comments
regarding the first category of BCs, boundary condition types and the resolution of boundary
condition conflicts will be made.

**Any Equation Boundary Conditions** - There are several boundary condition types that are
not necessarily best binned with a specific equation type. The FIX, GD_* and TABLE boundary
condition types are general and can be applied to any equation type. A general description of these
types (called Category 1) is given below.

**Boundary Condition Types** - Beyond the generalized boundary conditions types and the Dirichlet
types, *Goma* has strong-collocated, weak form, and several others that are intrinsic to the
Galerkin finite element method; these are applied in a variety of ways. Because of this, boundary
conditions at a given node might interact in ways that produce unexpected results. For this reason,
it is important to understand the differing methods of application that occur in *Goma* and how
each affects the other. In addition, by cleverly mixing boundary conditions, the analyst is often
able to achieve a desired result, but only if the nature of each boundary condition is understood.
Toward this end, the user will find a *special label* assigned to each boundary condition, which,
with the ensuing explanation below, will provide each user with an understanding of how that BC
is applied within *Goma*.

On each boundary condition card, the boundary condition type appears in the **Description/Usage**
section. These are the following boundary condition types that will be found here:

* **DIRICHLET (DC)**
* **STRONGLY INTEGRATED (SIC)**
* **STRONGLY INTEGRATED EDGE (SIC_EDGE)**
* **COLLOCATED (PCC)**
* **COLLOCATED EDGE (PCC_EDGE)**
* **WEAKLY INTEGRATED (WIC)**

The following sections discuss the method of application of each boundary condition type along
with the implications of using each.

**DIRICHLET (DC):**

In the hierarchy of boundary conditions, Dirichlet conditions are at the top.
Nothing trumps a Dirichlet conditions. A Dirichlet condition is applied by
discarding all mechanics information related to a particular field variable that has
been accumulated at a given node and replacing it with a direct assignment of the
nodal unknown of that field with a fixed *a priori* value. Algorithmically, applying
a Dirichlet condition on a degree of freedom at a node involves zeroing the entire equation row, inserting a unity value on the diagonal element of the Jacobian
matrix, inserting a zero value at the appropriate place in the residual vector, and
inserting the known boundary condition value at the appropriate place in the
solution vector. This is referred to in many places as the “*hard set*” method. An
alternate formulation imposes the boundary condition by replacing the mechanics
equation at a node with the simple *residual equation*, EQUATION , where φ and φ0
are the nodal unknown field and its assigned value, respectively.The sensitivities
of this residual equation are entered into the Jacobian appropriately and solution
takes place normally.

Dirichlet conditions are strictly node-based. Neighbor nodes and shared elements
have no influence on them. For this reason, all Dirichlet conditions are applied to
nodesets. Furthermore, Dirichlet conditions are assigned the highest precedence in
terms of boundary conditions. If a Dirichlet condition appears at a node, it will be
applied. Any other boundary condition that could be applied will be discarded (at
that node).

Dirichlet conditions are limited, however in that they can only affect the nodal
value of a degree of freedom. Derived quantities cannot be set with a Dirichlet
condition. You will never see a Dirichlet condition being applied to a heat flux for
example.

**STRONGLY INTEGRATED (SIC):**

The next class of boundary condition is referred to within *Goma* as the strongly
integrated boundary conditions. These boundary conditions replace the mechanics
equation at the i\ :sup:`th` node with a surface integral of some derived quantity. The
general form of these conditions is:

.. math::

   \int_{S} \phi_i \, g(\mathbf{x}) \, \mathrm{d}S = 0

where in this case :math:`g(\mathbf{x})` is not a residual equation but some derived quantity. Unlike
strong constraints, this term is not multiplied by a penalizing factor before it is
added to the accumulated mechanics equation at node *i*. Consequently, it
represents boundary contributions to the mechanics at that node. Note also that
since these conditions only make additions to the boundary mechanics, if a
strongly enforced condition (SIC or PCC) is also present at the node, the weakly
integrated constraint will be clobbered along with the rest of the mechanics. As an
example, a CAPILLARY boundary condition that is applied to the same sideset as
a VELO_NORMAL condition will have no effect in the final answer.

Weakly integrated boundary conditions are also very much a consequence of the
“natural” boundary conditions that emerge from the finite element formulation. As
anyone familiar with the finite element method knows, these are the ghostly
boundary terms that enforce zero boundary fluxes or forces as a convenient
default. The weakly integrated boundary condition step into the space afforded by
the natural boundary conditions and allow the user to specify values for these
boundary fluxes or forces as functions of conditions on those boundaries.

In addition, to the various classes of boundary conditions detailed above, there are special cases
that arise when applying boundary conditions to the “vector” degrees of freedom. Currently, the
only “vector” degrees of freedom are the mesh displacement and fluid velocity unknowns. When
a boundary condition is applied to these degrees of freedom, it may be ROTATED, VECTOR or
SCALAR. These labels appear in the boundary condition documentation along with the class of
the condition.

**ROTATED:**

When a boundary condition is designated as “ROTATED,” the vector components
of the appropriate equations for the surface nodes are projected into a new
coordinate system that is locally based on the surface normal vector and tangent
vectors. It is the presence of the “ROTATED” boundary condition that prompts
this process. Usually, only one of these rotated components is then affected by the
boundary condition constraint and in this sense ROTATED conditions are
SCALAR conditions (see below). Also generally speaking, ROTATED boundary
conditions are strongly enforced as described above.

**VECTOR:**

When a boundary condition is designated as a “VECTOR” condition, the
implication is that a vector quantity will be added to the vector components of the
original mechanics equations. “VECTOR” boundary conditions are generally
always applied weakly.

**SCALAR:**

When a boundary condition is designated a “SCALAR” condition, only a single
mechanics equation is going to be influenced by the boundary condition. In the
case of the vector degrees of freedom, only a single component would be affected
by the boundary condition. Boundary conditions that apply to degrees of freedom
that are naturally scalars, for instance temperature and species, are by default
SCALAR conditions.

An example of these special labels for the *VELO_NORMAL_EDGE* condition (found on the line
with the **Description/Usage** section header) is **PCC-EDGE/ROTATED MOMENTUM** indicating
a rotated collocated edge condition applied to the fluid momentum equation. Given this labeling
convention, boundary conditions which are not specified to be rotated or vector conditions can
be presumed to be unrotated scalar conditions. Boundary conditions that may be applied to any
equation are labeled “varied.”

The user will not find “periodic boundary conditions” discussed in this manual. Those interested
in such conditions should consult the Advanced Capabilities Manual (SAND2006-7304).

**Resolving Conflicts between Boundary Conditions** - In *Goma*, the bulk equations and
boundary conditions are evaluated on an element-by-element basis. After the residual and Jacobian
entries for the bulk equations have been calculated, the boundary conditions are used to modify
or replace the bulk entries where necessary. Often the selection of boundary conditions from the
input deck may cause two boundary conditions to be applied to the same equation (equation associated
with a nodal point); this is especially true at junction points. Frequently the multiple boundary
conditions perform the same function (i.e. duplicates) but in some important instances they
are different (i.e. conflicts). In *Goma*, a decision making process was developed for determining
which boundary conditions have priority. The flow chart for this decision-making is shown in Figure
4. While this process resolves boundary-condition conflicts, it **does not** eliminate the possibility
of setting boundary conditions that are incompatible and lead to errors in solving the problem.
However, this method should clarify how BC’s are chosen from the input deck and should enable
the user to determine why a given combination of boundary conditions does not work.

The flow chart in Figure 3 shows the procedure for resolving what boundary conditions get applied
to a given equation at a given node. The starting point assumes that a list of all the potential
boundary conditions for the equation are known. Boundary conditions in *Goma* fall into several
classes: *Dirichlet, Pointwise Collocation, Strong Integrated, Weak Integrated and Special conditions*,
in order of priority. For boundary conditions applied to vector equations (mesh or momentum),
a boundary condition can cause the bulk equations to be rotated prior to applying the
boundary condition; in conflicts between boundary conditions, conditions which do not rotate the
bulk equations (*unrotated conditions*) have priority over conditions which rotate the bulk equations
(*rotated conditions*). In certain cases (e.g. two PLANE conditions which intersect at a point),
conflicting boundary conditions can be checked to determine if they are duplicates, in which case
only the first of the duplicates in the input deck is applied. Most boundary conditions are designed
to apply by themselves, but a special class of boundary conditions, the generalized dirichlet (GD) conditions, are designed so that multiple GD conditions can apply along the same boundary and to
the same equation.

While running, *Goma* prints the results of conflict resolution for every node at which it found at
least two boundary conditions being applied to the same equation. The results indicate the node
number, equation type, boundary conditions chosen by *Goma*, and the side-set or node-set numbers
to which the boundary conditions apply. Thus to determine what boundary conditions are actually
used by *Goma*, carefully check the output from conflict resolution. Setting the *Debug_Flag*
= 1 causes *Goma* to print out more information regarding which boundary conditions apply and
which do not. Despite the complexity of the logic built into *Goma* to resolve conflicts between
boundary conditions, there are several combinations of boundary conditions that do not have a
clear resolution. It is up to the user to resolve the final conflicts.

And finally, the first (*Number of BC*) and last (*END OF BC*) boundary condition cards are a pair
and stand alone; the remaining cards belong to the categories of conditions discussed above. The
ordering of input cards within this collection of BC input records (i.e., section) is sequential and
some sections of interspersed comments accompany each boundary condition category.

.. toctree::
   :maxdepth: 1

   boundary_conditions/number_of_bc
   boundary_conditions/any
   boundary_conditions/mesh
   boundary_conditions/real_solid
   boundary_conditions/fluid_momentum
   boundary_conditions/energy
   boundary_conditions/mass
   boundary_conditions/continuity
   boundary_conditions/porous
   boundary_conditions/stress
   boundary_conditions/gradient
   boundary_conditions/shear_rate
   boundary_conditions/fill
   boundary_conditions/potential
   boundary_conditions/fluid_solid_interaction
   boundary_conditions/level_set
   boundary_conditions/shell
   boundary_conditions/acoustic
   boundary_conditions/turbulence
   boundary_conditions/end_of_bc

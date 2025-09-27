*******************
AC (Flux Condition)
*******************

::

    AC = FC <mat_id> <bc_id> <data_float_index> <FC_type> [species_id] <sideset> <flux_value>

-----------------------
Description / Usage
-----------------------

This card attaches a boundary condition float parameter to an integrated flux constraint 
on a boundary sideset. The flux constraint consists of requiring a specific value for 
some integrated quantity over the sideset, for example, force, heat flux or fluid 
flowrate. A number of standard flux constraints have been provided so that a user-
defined flux is not necessary. During the solution process the boundary condition 
parameter is allowed to vary as the additional degree of freedom associated with the 
flux constraint.

Definitions of the input parameters are as follows:

FC
    Mandatory string indicating that this augmenting 
    condition is of variety Flux Condition.

<mat_id>
    An integer parameter identifying the identification 
    number of the material on which the flux integration 
    will be performed. Material properties needed in 
    evaluating the flux constraint will be obtained from this 
    material and only those elements belonging to this 
    material will contribute to the flux integral

<bc_id>
    An integer parameter giving the index of the boundary 
    condition card whose parameter is being used as the 
    new degree of freedom. Numbering begins with zero, 
    starting with the first boundary condition in the input 
    file and proceeding sequentially upward with each read 
    boundary condition card.

<data_float_index>
    An integer parameter that identifies the boundary 
    condition parameter that will be varied. It is an index 
    that starts at zero with the leftmost float value on the 
    <bc_id> boundary condition card and increments 
    upward from right to left.

<FC_type>
    A string parameter specifying one of the standard flux
    constraint equations. Choices are (all in capitals):

    - FORCE_X
    - FORCE_Y
    - FORCE_Z
    - FORCE_NORMAL
    - FORCE_TANGENT1
    - FORCE_TANGENT2
    - HEAT_FLUX
    - VOLUME_FLUX
    - SPECIES_FLUX
    - CHARGED_SPECIES_FLUX
    - CURRENT
    - PORE_LIQ_FLUX
    - TORQUE
    - AVERAGE_CONC
    - SURF_DISSIP
    - AREA
    - VOL_REVOLUTION
    - CURRENT_FICKIAN
    - NEG_LS_FLUX
    - POS_LS_FLUX
    - N_DOT_X
    - ELEC_FORCE_NORMAL
    - ELEC_FORCE_TANGENT1
    - ELEC_FORCE_TANGENT2
    - ELEC_FORCE_X
    - ELEC_FORCE_Y
    - ELEC_FORCE_Z
    - NET_SURF_CHARGE
    - DELTA
    - ACOUSTIC_FLUX_NORMAL
    - ACOUSTIC_FLUX_TANGENT1
    - ACOUSTIC_FLUX_TANGENT2
    - ACOUSTIC_FLUX_X
    - ACOUSTIC_FLUX_Y
    - ACOUSTIC_FLUX_Z
    - ACOUSTIC_INTENSITY
    - LS_DCA

    See below for a detailed description of each of these 
    constraints.

<side set>
    An integer parameter that identifies the side set over 
    which the flux condition will be integrated. Those 
    elements that belong to <mat_id> and have faces 
    included in <side set> will be evaluated in the 
    integration.

[species_id]
    This is an integer parameter that identifies the species 
    component that is evaluated in the SPECIES_FLUX 
    constraint. When other flux constraints are used this 
    parameter should not appear.

<flux_value>
    A float parameter that specifies the value that is 
    assigned to the flux integral. At the end of the solution 
    procedure, the flux integral will be equal to this float 
    value.

------------
Examples
------------

The following are two examples of this AC type:

::

    Number of augmenting conditions = -1
    AC = FC 1 0 0 VOLUME_FLUX 4 {-PI}
    AC = FC 1 1 0 FORCE_TANGENT1 3 {3*5.0}
    END OF AC

-------------------------
Technical Discussion
-------------------------

The following list describes precisely some of the preceding FC_type constraint 
equations. In the following, :math:`\Phi` is the <flux_value> parameter and :math:`\Gamma` is the side set 
specified above.

For material blocks with ARBITRARY mesh motion:

**FORCE_X:**

.. math::

    \Phi = \int_\Gamma (\mathbf{n} \cdot \boldsymbol{\Pi} - \rho(\mathbf{u} - \mathbf{u}_m)(\mathbf{n} \cdot \mathbf{u})) \cdot \mathbf{e}_x d\Gamma

**FORCE_Y:**

.. math::

    \Phi = \int_\Gamma (\mathbf{n} \cdot \boldsymbol{\Pi} - \rho(\mathbf{u} - \mathbf{u}_m)(\mathbf{n} \cdot \mathbf{u})) \cdot \mathbf{e}_y d\Gamma

**FORCE_Z:**

.. math::

    \Phi = \int_\Gamma (\mathbf{n} \cdot \boldsymbol{\Pi} - \rho(\mathbf{u} - \mathbf{u}_m)(\mathbf{n} \cdot \mathbf{u})) \cdot \mathbf{e}_z d\Gamma

**FORCE_NORMAL:**

.. math::

    \Phi = \int_\Gamma (\mathbf{n} \cdot \boldsymbol{\Pi} - \rho(\mathbf{u} - \mathbf{u}_m)(\mathbf{n} \cdot \mathbf{u})) \cdot \mathbf{n} d\Gamma

**FORCE_TANGENT1:**

.. math::

    \Phi = \int_\Gamma (\mathbf{n} \cdot \boldsymbol{\Pi} - \rho(\mathbf{u} - \mathbf{u}_m)(\mathbf{n} \cdot \mathbf{u})) \cdot \mathbf{t}_1 d\Gamma

**FORCE_TANGENT2:**

.. math::

    \Phi = \int_\Gamma (\mathbf{n} \cdot \boldsymbol{\Pi} - \rho(\mathbf{u} - \mathbf{u}_m)(\mathbf{n} \cdot \mathbf{u})) \cdot \mathbf{t}_2 d\Gamma

**HEAT_FLUX:**

.. math::

    \Phi = \int_\Gamma (\mathbf{n} \cdot \mathbf{q} - \rho C(\mathbf{u} - \mathbf{u}_m) \cdot \mathbf{n}) d\Gamma

**VOLUME_FLUX:**

.. math::

    \Phi = \int_\Gamma (\mathbf{u} - \mathbf{u}_m) \cdot \mathbf{n} d\Gamma

**SPECIES_FLUX:**

.. math::

    \Phi = \int_\Gamma (\mathbf{n} \cdot \mathbf{j}_i - c_i(\mathbf{u} - \mathbf{u}_m) \cdot \mathbf{n}) d\Gamma

where:

- :math:`\boldsymbol{\Pi}` = fluid stress tensor
- :math:`\rho` = fluid phase density
- :math:`\mathbf{u}` = fluid velocity vector
- :math:`\mathbf{u}_m` = mesh velocity vector
- :math:`\mathbf{e}_x`, :math:`\mathbf{e}_y`, :math:`\mathbf{e}_z` = cartesian unit bases
- :math:`\mathbf{n}` = normal vector to side set
- :math:`\mathbf{t}_1` = first tangent vector to side set
- :math:`\mathbf{t}_2` = second tangent vector side set (3D only)
- :math:`\mathbf{q}` = diffusive heat flux
- :math:`C` = heat capacity per unit mass
- :math:`\mathbf{j}_i` = diffusive flux of species i
- :math:`c_i` = concentration of species i.

In the case of a LAGRANGIAN material, only the FORCE flux constraints are 
relevant. They have the same form as the preceding except that there are no 
convective contributions to the flux and Π represents the solid stress tensor.

For additional technical notes regarding this augmenting condition type see the AC 
(Boundary Condition) section.


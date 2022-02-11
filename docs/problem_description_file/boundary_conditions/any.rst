~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Category 1: Any Equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This category includes a set of cards that are used to provide all boundary condition information
for a generalized dirichlet (GD) boundary condition. The condition is applied as a pointwise
collocation along a given node set. The general syntax for the GD_cards is as follows:

::

	BC = <bc_name> <bc_type> <bc_id> <equation_name> <integer1> <variable_name> <integer2> {float_list}

The current allowable definitions and/or values for < bc_name>, <bc_type>, <bc_id>,
<integer1>, <integer2> and {float_list} are provided in the individual cards. As a general note,
<integer1> and <integer2> are the species number of the mass transport equation and
concentration variable, respectively; they should be zero for other equation and variable types.
Currently these conditions assume that the variable is defined at all the nodes at which the
equation is defined (no subparametric mapping).

However, the values for <equation_name> and <variable_name>, which apply generally to all
cards in this category (except as subsequently noted), are given here:

<equation_name>
    A character string indicating the equation to which this
    boundary condition is applied, which can be

    * ``R_MOMENTUM1 R_MOMENTUM2 R_MOMENTUM3``
    * ``R_MESH1 R_MESH2 R_MESH3``
    * ``R_MASS R_ENERGY R_MASS_SURF``
    * ``R_PRESSURE`` 
    * ``R_STRESS11 R_STRESS12 R_STRESS13``
    * ``R_STRESS22 R_STRESS23 R_STRESS33``
    * ``R_GRADIENT11 R_GRADIENT12``
    * ``R_GRADIENT13 R_GRADIENT21``
    * ``R_GRADIENT22 R_GRADIENT23``
    * ``R_GRADIENT31 R_GRADIENT32``
    * ``R_GRADIENT33 R_POTENTIAL R_FILL``
    * ``R_SHEAR_RATE R_MESH_NORMAL`` (rotate mesh equations and apply
      this condition to normal component)
    * ``R_MESH_TANG1 R_MESH_TANG2``
    * ``R_MOM_NORMAL`` (rotate momentum equations and apply this condition
      to normal component)
    * ``R_MOM_TANG1 R_MOM_TANG2``
    * ``R_POR_LIQ_PRESS R_POR_GAS_PRESS``
    * ``R_POR_POROSITY R_POR_SATURATION``
    * ``R_POR_ENERGY R_POR_LAST``
    * ``R_POR_SINK_MASS R_VORT_DIR1``
    * ``R_VORT_DIR2 R_VORT_DIR3 R_VORT_LAMBDA``
    * ``R_CURVATURE R_LAGR_MULT1``
    * ``R_LAGR_MULT2 R_LAGR_MULT3``
    * ``R_BOND_EVOLUTION R_SURF_CHARGE``
    * ``R_EXT_VELOCITY R_EFIELD1 R_EFIELD2``
    * ``R_EFIELD3 R_ENORM R_NORMAL1``
    * ``R_NORMAL2 R_NORMAL3 R_ _CURVATURE``
    * ``R_SHELL_TENSION R_SHELL_X R_SHELL_Y``
    * ``R_SHELL_USER R_PHASE1 R_PHASE2``
    * ``R_PHASE3 R_PHASE4 R_PHASE5``
    * ``R_SHELL_ANGLE1 R_SHELL_ANGLE2``
    * ``R_SHELL_SURF_DIV_V R_SHELL_SURF_CURV``
    * ``R_N_DOT_CURL_V R_GRAD_S_V_DOT_N1``
    * ``R_GRAD_S_V_DOT_N2 R_GRAD_S_V_DOT_N3``
    * ``R_ACOUS_PREAL R_ACOUS_PIMAG``
    * ``R_SHELL_DIFF_FLUX``
    * ``R_SHELL_DIFF_CURVATURE``
    * ``R_SHELL_NORMAL1 R_SHELL_NORMAL2``
    * ``R_ACOUS_REYN_STRESS R_SHELL_BDYVELO``
    * ``R_SHELL_LUBP R_LUBP R_SHELL_FILMP``
    * ``R_SHELL_FILMH R_SHELL_PARTC``
    * ``R_SHELL_SAT_CLOSED R_SHELL_SAT_OPEN``
    * ``R_SHELL_ENERGY R_SHELL_DELTAH``
    * ``R_SHELL_LUB_CURV R_SHELL_SAT_GASN``
    * ``R_SHELL_SHEAR_TOP R_SHELL_SHEAR_BOT``
    * ``R_SHELL_CROSS_SHEAR R_MAX_STRAIN``
    * ``R_CUR_STRAIN  R_LUBP_2``
    * ``R_SHELL_SAT_OPEN_2 or``
    * ``R_SHELL_LUB_CURV_2``

<variable_name>
    A character string indicating the variable which should be
    fixed, which can be

    * ``VELOCITY1 VELOCITY2 VELOCITY3``
    * ``MESH_DISPLACEMENT1``
    * ``MESH_DISPLACEMENT2``
    * ``MESH_DISPLACEMENT3 MESH_POSITION1``
    * ``MESH_POSITION2 MESH_POSITION3``
    * ``MASS_FRACTION SURFACE TEMPERATURE or``
    * ``PRESSURE`` (pressure will have no effect if not using Q1 or Q2 basis functions)
    * ``POLYMER_STRESS11``
    * ``POLYMER_STRESS12 POLYMER_STRESS13``
    * ``POLYMER_STRESS22 POLYMER_STRESS23``
    * ``POLYMER_STRESS33 VOLTAGE FILL``
    * ``SHEAR_RATE VEL_NORM D_VEL1_DT``
    * ``D_VEL2_DT D_VEL3_DT D_T_DT D_C_DT``
    * ``D_X1_DT D_X2_DT D_X3_DT D_S_DT D_P_DT``
    * ``VELOCITY_GRADIENT11``
    * ``VELOCITY_GRADIENT12``
    * ``VELOCITY_GRADIENT13``
    * ``VELOCITY_GRADIENT21``
    * ``VELOCITY_GRADIENT22``
    * ``VELOCITY_GRADIENT23``
    * ``VELOCITY_GRADIENT31``
    * ``VELOCITY_GRADIENT32``
    * ``VELOCITY_GRADIENT33 POR_LIQ_PRESS``
    * ``POR_GAS_PRESS POR_POROSITY``
    * ``POR_POROSITY POR_TEMP  POR_SATURATION``
    * ``POR_LAST MAX_POROUS_NUM``
    * ``POR_SINK_MASS VORT_DIR1 VORT_DIR2``
    * ``VORT_DIR3 VORT_LAMBDA CURVATURE``
    * ``LAGR_MULT1 LAGR_MULT2 LAGR_MULT3``
    * ``BOND_EVOLUTION SURF_CHARGE``
    * ``EXT_VELOCITY EFIELD1 EFIELD2 EFIELD3``
    * ``ENORM NORMAL1 NORMAL2 NORMAL3``
    * ``SHELL_CURVATURE SHELL_TENSION``
    * ``SHELL_X SHELL_Y SHELL_USER PHASE1``
    * ``PHASE2 PHASE3 PHASE4 PHASE5``
    * ``SHELL_ANGLE1 SHELL_ANGLE2``
    * ``SHELL_SURF_DIV_V SHELL_SURF_CURV``
    * ``N_DOT_CURL_V GRAD_S_V_DOT_N1``
    * ``GRAD_S_V_DOT_N2 GRAD_S_V_DOT_N3``
    * ``ACOUS_PREAL ACOUS_PIMAG``
    * ``SHELL_DIFF_FLUX SHELL_DIFF_CURVATURE``
    * ``SHELL_NORMAL1 SHELL_NORMAL2``
    * ``ACOUS_REYN_STRESS SHELL_BDYVELO``
    * ``SHELL_LUBP LUBP SHELL_FILMP``
    * ``SHELL_FILMH SHELL_PARTC``
    * ``SHELL_SAT_CLOSED SHELL_PRESS_OPEN``
    * ``SHELL_TEMPERATURE SHELL_DELTAH``
    * ``SHELL_LUB_CURV SHELL_SAT_GASN``
    * ``SHELL_SHEAR_TOP SHELL_SHEAR_BOT``
    * ``SHELL_CROSS_SHEAR MAX_STRAIN``
    * ``CUR_STRAIN LUBP_2 SHELL_PRESS_OPEN2``
    * ``SHELL_LUB_CURV_2``

EXCEPTIONS to the above parameter definitions: For the *GD_TIME* card, the <variable_names>
of **LINEAR, EXPONENTIAL**, or **SINUSOIDAL** are acceptable (see examples below). There
are also differences in the use of the *GD_TABLE* card, which are explained in the description of
that card below.

A GD boundary condition can be applied multiple times to the same side set and equation to build
up a general multiparameter condition. When this is done, the function is built by expanding the
equations sequentially in the order specified in the BC list.

Descriptions of the GD cards are given next. An insert entitled “Usage Notes on the GD Cards”
follows the descriptions, explaining how the cards are used together in various combinations.

.. include:: /problem_description_file/boundary_conditions/any/fix.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/gd_const.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/gd_linear.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/gd_parab.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/gd_polyn.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/gd_time.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/gd_circ.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/gd_table.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/gd_usage.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/table_wicv.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/table_wics.rst

---------------------------------------------------------------------------------------

.. include:: /problem_description_file/boundary_conditions/any/table.rst


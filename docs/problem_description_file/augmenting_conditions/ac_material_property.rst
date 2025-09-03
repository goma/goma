**********************
AC (Material property)
**********************

::

    AC = MT <mat_id> <material_prop_tag> <float_data_list>

-----------------------
Description / Usage
-----------------------

This augmenting condition type attaches a material property (viscosity, thermal 
conductivity, surface tension, etc.) to an augmenting constraint. The value of the 
material property value is permitted to vary as any other degree of freedom in order that 
the user-supplied augmenting condition is satisfied.

Definitions of the input parameters are as follows:

**MT**
    A string designator identifying the AC card as a material 
    property type.

**<mat_id>**
    An integer parameter identifying the material whose 
    property is being varied.

**<material_prop_tag>**
    An integer parameter that associates an integer tag with 
    a type of material property. The material property value 
    designated by this tag will be varied in the solution 
    process. The following tables show the association 
    between material property and integer tag.

**<float_data_list>**
    A list of float parameters that can be used in
    user_aug_cond_residuals to evaluate the 
    augmenting conditions. They are stored in sequence in 
    the array augc[i].DataFlt.

**General Physical Properties**

=============================  ====
Property                       TAG
=============================  ====
THERMAL_CONDUCTIVITY           1100
ELECTRICAL_CONDUCTIVITY        1200
VISCOSITY                      1300
SURFACE_TENSION                1400
HEAT_CAPACITY                  1500
VOLUME_EXPANSION               1600
DENSITY                        1700
POROSITY                       1800
PERMEABILITY                   1900
REL_GAS_PERM                   2000
REL_LIQ_PERM                   2100
SATURATION                     2200
MELTING_POINT_LIQUIDUS         2500
MELTING_POINT_SOLIDUS          2600
FLOWINGLIQUID_VISCOSITY        2700
=============================  ====

**Generalized Newtonian Models**

===============  ====
Property         TAG
===============  ====
MU0              4000
NEXP             4100
MUINF            4200
LAM              4300
AEXP             4400
ATEXP            4500
WLFC2            4550
TAU_Y            4600
FEXP             4610
MAXPACK          4800
FICKDIFF_X       4810
FICKDIFF_Y       4820
===============  ====

**Viscoelastic Models**

===============  ========
Property         TAG
===============  ========
TIME_CONST       5000\ :sup:`a`
WT_FUNC          5100\ :sup:`a`
ALPHA            5200\ :sup:`a`
PTT_XI           5300\ :sup:`a`
PTT_EPS          5400\ :sup:`a`
===============  ========

\ :sup:`a` This is the TAG_ID for the first viscoelastic mode; subsequent mode 
TAG_IDs are incremented by 1, e.g., TAGC_TIME_CONST(mode i) 
= TAGC_TIME_CONST + i

**Solid Elastic Models**

===============================  ====
Property                         TAG
===============================  ====
LAME_MU                          6000
LAME_MU_CONTACT_LINE_G0          6001
LAME_MU_CONTACT_LINE_G1          6002
LAME_MU_CONTACT_LINE_R0          6003
LAME_LAMBDA                      6100
CONV_LAG_VELX                    6201
CONV_LAG_VELY                    6202
CONV_LAG_VELZ                    6203
CONV_LAG_ROTRATE                 6221
CONV_LAG_ROT_X0                  6222
CONV_LAG_ROT_Y0                  6223
CONV_LAG_ROT_Z0                  6224
RS_LAME_MU                       6300
RS_LAME_LAMBDA                   6400
POISSON                          6600
STRSS_FR_SOL_VOL_FRAC            6610
===============================  ====

This is only a partial list of the available material property tags. To find the latest tag 
definitions (those active in your version of Goma), look in the Goma header file 
mm_mp_const.h. The initial value of the extra unknowns is taken from those values 
given in the input, material definition, and geometry definition (if any) files.

------------
Examples
------------

The following is an example this AC card:

::

    AC = MT 1 1400

This card indicates that this augmenting condition constraint is associated with the 
surface tension value found in the material file designated for material 1.

Note that the value of surface tension originally supplied in the material final will be 
used as a starting guess. However, at the end of the computation a different value for 
surface tension will have been determined and save in the memory location for this 
material property. If an ASCII output file is requested (SOLN file = ), the value of the 
surface tension will appear at the end of the file following the output of the nodal 
unknown vector.

-------------------------
Technical Discussion
-------------------------

See the technical discussion appearing in the documentation for the AC(Boundary 
Condition) card.

----------
Theory
----------

No Theory.

-------
FAQs
-------

No FAQs.

--------------
References
--------------

Gates, I.D., D.A. Labreche, and M.M. Hopkins, "Advanced Capabilities in GOMA 3.0 
- Augmenting Conditions, Automatic Continuation, and Linear Stability Analysis," 
SAND2000-2465, Albuquerque, NM, (2001).
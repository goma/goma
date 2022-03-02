********************************
**Shell Energy Source External**
********************************

::

   Shell Energy Source Viscous External = {model_name} <float_list>

-----------------------
**Description / Usage**
-----------------------

This card activates a heat source (or sink, as it were) in the shell_energy equation
which corresponds to a user-supplied or constant value. Two models are available.

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model invokes a constant heat source term (heat sink if negative) in units of   |
|                          |energy per area per time.                                                            |
|                          |                                                                                     |
|                          | * <float1> - Value of heat source.                                                  |
+--------------------------+-------------------------------------------------------------------------------------+
|**JOULE**                 |This model invokes a constant energy source which is determined by an external       |
|                          |current density field of the form                                                    |
|                          |                                                                                     |
|                          |Qjoule = hξ J ⋅ J. Here J is the current density, h is the gap, and is the electrical|
|                          |resistivity, or the inverse conductivity. Both h and ξ are determined from other     |
|                          |models in the material file. J is brought in as an external field variable from      |
|                          |another exodusII file (see discussion below).                                        |
|                          |                                                                                     |
|                          | * <float1> - Scale factor, usually set to 1.0.                                      |
+--------------------------+-------------------------------------------------------------------------------------+
|JOULE_LS                  |This model differs from JOULE only in that the electrical conductivity is pulled out |
|                          |and must be specified with a LEVEL_SET model. This model is not well tested (PRS     |
|                          |12/14/2012)                                                                          |
|                          |                                                                                     |
|                          | * <float1> - Scale factor, usually set to 1.0.                                      |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Shell Energy Source External = JOULE {scale=1.0}

-------------------------
**Technical Discussion**
-------------------------

To bring in an external field of the appropriate form, see the main Goma user’s manual
and refer to the External Field card. As an example, you might consider solving
a simple electrostatic problem using the EQ = V (voltage) equation and output the
magnitude of the current density vector. In Goma, this is done with the post
processing
::

   Electric Field Magnitude = yes

Card. This card outputs this J-magnitude as the exodusII variable EE. You then bring
it in as follows in the input script:

::

   External Field = EE Q1 current_dens_out.exoII

With the JOULE model, this field is used to compute the Joule heating term.



--------------
**References**
--------------

No References.
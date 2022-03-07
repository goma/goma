*******************************
**Lubrication Momentum Source**
*******************************

::

   Lubrication Momentum Source = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card sets a fluid “body force per unit volume” source term in the lub_p
equation. This capability can be used to specify a force field over the entire shell area
(over the shell material).. Currently two models {model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |THIS MODEL NOT IMPLEMENTED AS OF 11/11/2010. This model is used to set a constant    |
|                          |fluid momentum source in units of force per unit volume. Only one floating point     |
|                          |value is required.                                                                   |
|                          |                                                                                     |
|                          | * <float1> is the fluid momentum source in F/L^3.                                   |
+--------------------------+-------------------------------------------------------------------------------------+
|**JXB**                   |This model is used to set fluid momentum source in units of force per unit volume    |
|                          |which comes from externally supplied current density J field and magnetic B fields.  |
|                          |These fields are suppled with the external field capability in Goma in a component   |
|                          |wise fashion. Please consult the technical discussion below.                         |
|                          |                                                                                     |
|                          | * <float1> is scale factor which may be used for nondimensionalization. Typically   |
|                          |   this is set to 1.0.                                                               |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Lubrication Momentum Source = JXB 1.

-------------------------
**Technical Discussion**
-------------------------

The two vector fields J, the current flux, and B, the magnetic induction, must be
supplied to Goma in order to activate this option. At present, these fields must be
supplied with the External Field cards, which provide the specific names of nodal
variable fields in the EXODUS II files from which the fields are read. The three
components of the J field must be called JX_REAL, JY_REAL, and JZ_REAL.
Likewise the B field components must be called BX_REAL, BY_REAL, and
BZ_REAL. These names are the default names coming from the electromagnetics code
like Alegra. Because of the different coordinate convention when using cylindrical
components, the fields have been made compatible with those arising from TORO II. It
is the interface with TORO that also makes the Lorentz scaling (lsf) necessary so that
the fixed set of units in TORO (MKS) can be adjusted to the user-selected units in
Goma.



--------------
**References**
--------------

No References.
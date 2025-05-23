********************
**Default Database**
********************

::

   Default Database = {GOMA_MAT | CHEMKIN_MAT}

-----------------------
**Description / Usage**
-----------------------

This card sets the default material database type. The default for this card is
GOMA_MAT. In that case, all material properties for the current material are obtained
from the current material file being read. If the default database is Chemkin, then the
Chemkin 3 linking files are read in, and initialization of most of the methods and data
for thermodynamic function evaluation, the stoichiometry and names of species and
elements, the homogeneous and heterogeneous source terms for chemical reactions and
their coupling into the energy equation, and transport property evaluations occurs.
Many fields in the materials database file that were required now are optional. After
Chemkin initialization, the rest of the materials database file is then read in. At that
time, some fields containing methods and data that were initialized with Chemkin
methods and data may be overwritten with methods and data specified by the material
file. Other fields not initialized or even handled by Chemkin (such as surface tension)
must be initialized for the first time by the materials file. Thus, the use of Chemkin
materials database doesn’t mitigate the need for a Goma materials file.

------------
**Examples**
------------

Following is a sample card:
::

   Default Database = CHEMKIN_MAT

-------------------------
**Technical Discussion**
-------------------------

Chemkin includes its own rigorous treatment of ideal gas thermodynamics and
transport property evaluations, providing it with a solid foundation on which to build
kinetics mechanisms and a rigorous treatment of gas phase transport property
evaluation. In order to maintain internal consistency, the new treatment must be used in
its entirety.



--------------
**References**
--------------

No References.
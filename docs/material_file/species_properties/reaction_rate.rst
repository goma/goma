*****************
**Reaction Rate**
*****************

::

   Reaction Rate = <model_name> <float1> <float2>

-----------------------
**Description / Usage**
-----------------------

This card is used to specify rates of species electrochemical reactions in the anode and
cathode regions in a LiSi/LiCl-KCl/FeS2 thermal battery cell using Butler-Volmer
kinetics.

This property currently allows for a single {model_name} which has two parameters:

+----------------------+-------------------------------------------------------------------------------------+
|**ELECTRODE_KINETICS**|the name of reaction rate model                                                      |
|                      |                                                                                     |
|                      | * <float1> - Anodic direction transfer coefficient                                  |
|                      | * <float2> - Cathodic direction transfer coefficient                                |
+----------------------+-------------------------------------------------------------------------------------+

Two companion cards, *THERMODYNAMIC POTENTIAL* and *INTERFACIAL AREA*
are required to complete the specification of parameters present in the Butler-Volmer
kinetic model of current density.

------------
**Examples**
------------

The following are two sample cards:

::

   Reaction Rate = ELECTRODE_KINETICS 0.5 0.5

::

   Reaction Rate = ELECTRODE_KINETICS 1.0 1.0

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.
**********************
**Chemical Potential**
**********************

::

   Chemical Potential = {IDEAL_SOLUTION | STOICHIOMETRIC_PHASE}

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the formulation of the chemical potential for the phase. It is
currently unconnected to *Goma’s* functionality. Two values are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**IDEAL_SOLUTION**        |Ideal solution thermodynamics                                                        |
+--------------------------+-------------------------------------------------------------------------------------+
|**STOICHIOMETRIC_PHASE**  |Phase consists of fixed set of molecular composition                                 |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Chemical Potential = IDEAL_SOLUTION

-------------------------
**Technical Discussion**
-------------------------

The chemical potential of species k in an ideal solution is given by the expression,
[Denbigh, p. 249],

.. figure:: /figures/457_goma_physics.png
	:align: center
	:width: 90%

where μk*(T, P) is defined as the chemical potential of species k in its pure state (or a
hypothetical pure state if a real pure state doesn’t exist) at temperature T and pressure
P. μk*(T, P) is related to the standard state of species k in the phase, μk, o(T), which
is independent of pressure, through specification of the pressure dependence of the
pure species k. Xk is the mole fraction of species k in the phase.

The chemical potential of species k (actually there is only one species!) in a
stoichiometric phase is equal to

.. figure:: /figures/458_goma_physics.png
	:align: center
	:width: 90%



--------------
**References**
--------------

Denbigh, K., The Principles of Chemical Equilibrium, 4th Ed., Cambridge University
Press, 1981
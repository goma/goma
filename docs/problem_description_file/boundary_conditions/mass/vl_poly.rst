***********
**VL_POLY**
***********

::

	BC = VL_POLY SS <bc_id> {char_string} <integer_list> <float>

-----------------------
**Description / Usage**
-----------------------

**(SIC/MASS)**

This boundary condition card enforces vapor-liquid equilibrium between a gas phase
and a liquid phase using Flory-Huggins activity expression to describe polymer-solvent
mixtures. The condition only applies to interphase mass, heat, and momentum transfer
problems with discontinuous (or multivalued) variables at an interface, and it must be
invoked on fields that employ the **Q1_D** or **Q2_D** interpolation functions to “tie”
together or constrain the extra degrees of freedom at the interface in question.

There are three input values in the <integer_list>; definitions of the input parameters
are as follows:

+-------------+---------------------------------------------------------------+
|**VL_POLY**  | Name of the boundary condition (<bc_name>).                   |
+-------------+---------------------------------------------------------------+
|**SS**       | Type of boundary condition (<bc_type>), where **SS** denotes  |
|             | side set in the EXODUS II database.                           |
+-------------+---------------------------------------------------------------+
|<bc_id>      | The boundary flag identifier, an integer associated with      |
|             | <bc_type> that identifies the boundary location (side set in  |
|             | EXODUS II) in the problem domain.                             |
+-------------+---------------------------------------------------------------+
|{char_string}| the concentration basis; two options exist:                   |
|             |                                                               |
|             |    * **MASS** - the concentration variable in *Goma* is       |
|             |      equivalent to mass fractions.                            |
|             |    * **VOLUME** - the concentration variable in *Goma* is     |
|             |      based on volume fractions for all species.               |
|             |                                                               |
+-------------+---------------------------------------------------------------+
|<integer1>   | Species number of concentration.                              |
+-------------+---------------------------------------------------------------+
|<integer2>   | Element block id that identifies the liquid phase.            |
+-------------+---------------------------------------------------------------+
|<integer3>   | Element block id that identifies the vapor phase.             |
+-------------+---------------------------------------------------------------+
|<float>      | Total pressure of the system.                                 |
+-------------+---------------------------------------------------------------+

------------
**Examples**
------------

This is a sample input card for this boundary condition:
::

   BC = VL_POLY SS 7 MASS 0 1 2 1.e+05

-------------------------
**Technical Discussion**
-------------------------

One of the simplest forms of the equilibrium relation is the Raoult’s law, where the
mole fraction of a species is equal to its mole fraction in the liquid multiplied by the ratio of its pure component vapor pressure to the total pressure in the system.

.. math::

   y_iP^{total}=\phi_i\gamma_iP_i^v\quad\forall i

.. math::

   \sum_{i=1}^N\phi_i=1

:math:`\gamma_i` is defined as the activity coefficient of species *i* and is considered a departure
function from the Raoult’s law. The fugacity in the liquid is reformulated in terms of
volume fraction :math:`\phi_i` for polymer mixtures to avoid referencing the molecular weight of polymer (Patterson, et. al., 1971).

Based on an energetic analysis of excluded volume imposed by the polymer, the
activity coefficient model of Flory-Huggins is widely used for polymer-solvent
mixtures (Flory, 1953). The general form of the Flory-Huggins model for
multicomponent mixtures is a summation of binary interactions terms; i.e.,

.. math::

   \ln\gamma_i=\sum_{k=1}^N\left[\left(\frac{\tilde{v}_i}{\tilde{v}_k}\right)(\delta_{ki}-\phi_i)\right]+\sum_{k=2}^N\sum_{j=1}^{k-1}\left[\delta_{ij}\phi_k\chi_{jk}+\left(\frac{\tilde{v}_i}{\tilde{v}_j}\right)\phi_j(\delta_{ki}-\phi_k)\chi_{jk}\right]

:math:`v_i` is the molar volume of component *i* (or the average-number molar volume if *i* is a polymer). :math:`\zeta_{ki}` is the Dirac delta. :math:`\chi_{jk}` is known as the Flory-Huggins interaction
parameter between components *j* and *k*, and is obtainable by fitting the solubility data
to the above model. For a simple binary pair (solvent (1)-polymer (2)) and assuming
:math:`v_2` » :math:`v_1`, the above model reduces to a simpler form.

.. math::

   \ln\gamma_1=\phi_2+\chi_{12}\phi_2^2



--------------
**References**
--------------

Flory, P., Principles of Polymer Chemistry, Cornell University Press, New York (1953)

Patterson, D., Y.B. Tewari, H.P. Schreiber, and J.E. Guillet, "Application of Gas-Liquid
Chromatography to the Thermodynamics of Polymer Solutions,"Macromolecules, 4, 3,
356-358 (1971)

GTM-007.1: New Multicomponent Vapor-Liquid Equilibrium Capabilities in GOMA,
December 10, 1998, A. C. Sun

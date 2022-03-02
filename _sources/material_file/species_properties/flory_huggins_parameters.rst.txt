****************************
**Flory-Huggins parameters**
****************************

::

   Flory-Huggins parameters = CONSTANT <integer1> <integer2> <float>

-----------------------
**Description / Usage**
-----------------------

This card specifies the Flory-Huggins binary interaction parameters. It is assumed that
the binary parameters are symmetric; i.e.,

.. figure:: /figures/462_goma_physics.png
	:align: center
	:width: 90%

Therefore, one set of i-j coefficients is sufficient to describe the binary interaction
coefficients.

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |Model for constant Flory-Huggins parameters.                                         |
|                          |                                                                                     |
|                          | * <integer1> - first species number.                                                |
|                          | * <integer2> - second species number.                                               |
|                          | * <float> - Flory-Huggins binary interaction coefficient.                           |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is an example set of cards for a three-species mixture:

::

   Flory-Huggins parameters = CONSTANT 0 1 0.3
   Flory-Huggins parameters = CONSTANT 0 2 0.3
   Flory-Huggins parameters = CONSTANT 1 2 0.3

This example shows that two species are solved in the *Goma* problem explicitly:
species 0 and species 1.

-------------------------
**Technical Discussion**
-------------------------

No discussion; see Sun (1998).



--------------
**References**
--------------

GTM-007.1: New Multicomponent Vapor-Liquid Equilibrium Capabilities in GOMA,
December 10, 1998, A. C. Sun
************************
**LATENT_HEAT_INTERNAL**
************************

::

	BC = LATENT_HEAT_INTERNAL SS {char_string} <integer_list> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition card is used for latent heat release/adsorption at an internal
interface. See usage comments in the Technical Discussion.

The <integer_list> requires two values be specified; definitions of the input 
parameters are as follows:

+------------------------+----------------------------------------------------------+
|**LATENT_HEAT_INTERNAL**| Name of the boundary condition (<bc_name>).              |
+------------------------+----------------------------------------------------------+
|**SS**                  | Type of boundary condition (<bc_type>), where **SS**     |
|                        | denotes side set in the EXODUS II database.              |
+------------------------+----------------------------------------------------------+
|<bc_id>                 | The boundary flag identifier, an integer associated with |
|                        | <bc_type> that identifies the boundary location (side set|
|                        | in EXODUS II) in the problem domain.                     |
+------------------------+----------------------------------------------------------+
|{char_string}           | Variable name with the following permissible values:     |
|                        |                                                          |
|                        |   * **LIQUID_VAPOR**                                     |
|                        |   * **SOLID_LIQUID**                                     |
|                        |                                                          |
+------------------------+----------------------------------------------------------+
|<integer1>              | NOT ACTIVE. Any integer will do.                         |
+------------------------+----------------------------------------------------------+
|<integer2>              | NOT ACTIVE. Any integer will do.                         |
+------------------------+----------------------------------------------------------+
|<float3>                | Value of latent heat of vaporization/fusion for a pure   |
|                        | material case, in units of Energy/mass.                  |
+------------------------+----------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:
::

   BC = LATENT_HEAT_INTERNAL SS 40 SOLID_LIQUID 1 2 2.6e5

-------------------------
**Technical Discussion**
-------------------------

The *LATENT_HEAT_INTERNAL* card should be used for internal surfaces, or
interfaces, at which transfer is governed by actual physics being modeled as a part of
the problem. See *LATENT_HEAT* card for further information.

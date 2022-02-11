*******************************************
**Shell Energy Source Viscous Dissipation**
*******************************************

::

   Shell Energy Source Viscous Dissipation = {model_name} <float_list>

-----------------------
**Description / Usage**
-----------------------

This card activates a heat source (or sink, as it were) in the shell_energy
equation resulting from viscous dissipation due to shear combined Couette and
pressure-driven flow in the Reynolds lubrication equation (lubp equation).

+--------------------------+-------------------------------------------------------------------------------------+
|**LUBRICATION**           |This model invokes a viscous dissipation model simplified for the lubrication        |
|                          |approximation.                                                                       |
|                          |                                                                                     |
|                          | * <float1> - Scale factor for the term, typically taken as 1.0                      |
+--------------------------+-------------------------------------------------------------------------------------+
|**LUBRICATION_FRICTION**  |This model invokes the same viscous dissipation source term as in the LUBRICATION    |
|                          |model, but adds on an additional linear friction model of the form μ*Pload vslider.  |
|                          |Commensurately there are seven floats required:                                      |
|                          |                                                                                     |
|                          | * <float1> - coefficient of friction, μ.                                            |
|                          | * <float2> - Applied external load Pload                                            |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Shell Energy Source Viscous Dissipation= LUBRICATION_FRICTION {load=5e8} {coeff=0.9}




--------------
**References**
--------------

No References.
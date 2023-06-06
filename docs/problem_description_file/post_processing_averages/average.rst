**************
AVERAGE
**************

::

   AVERAGE = {average_type} <species_no>

-----------------------
Description / Usage
-----------------------

The AVERAGE card calculates average values at nodes and outputs to the Exodus file.
Definitions of the input parameters are as follows:

<average_type>
   Several choices are available:

   * **EM** Total electromagnetic wave, 6 fields (real and imaginary of 3 field vector)
   * **EMSCAT** Scattered electromagnetic wave (Total - Incident), 6 fields (real and imaginary of 3 field vector)
   * **EMINC** Incident electromagnetic wave, 6 fields (real and imaginary of 3 field vector)
   * **EM_MAG** Total electromagnetic wave magnitude, 1 field
   * **EM_SCAT_MAG** Scattered electromagnetic wave magnitude, 1 field
   * **EM_INC_MAG** Incident electromagnetic wave magnitude, 1 field
   * **TEMPERATURE** Temperature
   * **VISCOSITY** Viscosity

<species_no>
   The species number for **SPECIES_MASS**, usually 0 

------------
**Examples**
------------

::

   Post Processing Averages =
   AVERAGE = EM 0
   AVERAGE = EMSCAT 0
   AVERAGE = EMINC 0
   AVERAGE = EM_INC_MAG 0
   AVERAGE = EM_MAG 0
   AVERAGE = EM_SCAT_MAG 0
   END OF AVERAGES

-------------------------
**Technical Discussion**
-------------------------


--------------
**References**
--------------

No References.
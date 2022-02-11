**************************************************
Restart Time Integration After Renormalization
**************************************************

::

	Restart Time Integration After Renormalization = {yes | no}

-----------------------
Description / Usage
-----------------------

This card is used to specify whether or not to restart time integration each time Goma
renormalizes the level set function during the course of the computation. When time
integration is restarted, the time step is reset to its initial size and held at this step size
for the following 3 time steps. If this card is not present, the default is yes (time
integration will be restarted after each renormalization). The syntax of this card is as
follows:

{yes | no}
    Indicates the specified choice. {yes | on | true}can all be used to specify
    restarting of time integration. {no | off | false} can all be used to
    specify no restart.

------------
Examples
------------

This is a sample renormalization method input card:
::

	Restart Time Integration After Renormalization = no


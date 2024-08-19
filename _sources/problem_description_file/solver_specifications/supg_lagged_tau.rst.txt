***************
 Lagged Tau
***************

::

    SUPG Lagged Tau = {yes | no}

-----------------------
Description / Usage
-----------------------

Lag the  SUPG :math:`\tau`  (se the value fromt he previous step values).


yes
    Lag the :math:`\tau` computation
no
    Use the current step to compute :math:`\tau`

The default value is **no**.

------------
Examples
------------

The following is a sample card:
::

    SUPG Lagged Tau = yes

-------------------------
Technical Discussion
-------------------------



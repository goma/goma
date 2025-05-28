***************************************
Pressure Stabilization Disable Tau Sens
***************************************

::

    Pressure Stabilization Disable Tau Sens = {yes | no}

-----------------------
Description / Usage
-----------------------

Disable the sensitivities of PSPG :math:`\tau` scaling term in the Jacobian.


yes
    Omit the sensitivities of :math:`\tau`
    from the Jacobian
no
    Form the complete Jacobian.

The default value is **no**.

------------
Examples
------------

The following is a sample card:
::

    Pressure Stabilization Disable Tau Sens = yes

-------------------------
Technical Discussion
-------------------------



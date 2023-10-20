******************************************
Strong Boundary Condition Replace Equation
******************************************

::

    Strong Boundary Condition Replace Equation = {yes | no}

-----------------------
Description / Usage
-----------------------

This optional card  will replace the boundary residual contributions rather than use a penalty approach

yes
    Replace the boundary residual contributions.
no
    Use penalty method (Default)


The default value is **no**.

------------
Examples
------------

The following is a sample card:
::

	Strong Boundary Condition Replace Equation = yes

-------------------------
Technical Discussion
-------------------------

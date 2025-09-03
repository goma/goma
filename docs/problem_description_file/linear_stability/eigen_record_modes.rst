Eigen Record modes
==================

.. code-block:: none

    Eigen Record modes = <int>

Description/Usage
-----------------

This optional card determines how many of the converged modes should have their 
eigenvectors written to exodus II files. The valid input is:

0<=n<=N
    This requests that n exodus II outputs be written. n can be no larger than N, the value supplied in the Eigen Number of modes card.

The default value of n is 0.

Examples
--------

Here is a sample card:

::

    Eigen Record modes = 5

Technical Discussion
--------------------

This card is especially important when the 3D stability of a 2D flow is being computed. 
Each of the requested normal mode wave numbers receives n outputs, so the number of 
files written can become quite large if the user isn't careful.

This card is applicable to both eggroll and ARPACK eigensolvers. When eggroll is 
used, the name of the output file is LSA_<i>_of_n_<out_name>, where <i> is the ith 
mode (of n), and <out_name> is the name of the regular Exodus II output file (from the 
Output Exodus II file card), including any extension if there was one. In the case of 
3D stability of a 2D flow, the output file is LSA_<i>_of_n_wn=<f>_<out_name>, 
where <f> corresponds to the value of the normal mode (see the Eigen Wave Numbers
card). When ARPACK is used, these files will have the base name specified by the 
Eigenvector Output File card, similarly augmented with mode and wave numbers.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.

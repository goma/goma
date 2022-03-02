************************
**Q Tensor Diffusivity**
************************

::

   Q Tensor Diffusivity = <integer> <float_list>

-----------------------
**Description / Usage**
-----------------------

This card specifies the coefficients for use in the Q-tensor suspension rheology model.
The <float_list> has three values, one for each direction, so the input parameters are
defined as follows:

+----------------------+-------------------------------------------------------------------------------------+
|<integer>             |Species number for suspension volume fraction.                                       |
+----------------------+-------------------------------------------------------------------------------------+
|<float1>              |Coefficient of eigenvectors in the flow direction.                                   |
+----------------------+-------------------------------------------------------------------------------------+
|<float2>              |Coefficient of eigenvectors in the gradient direction.                               |
+----------------------+-------------------------------------------------------------------------------------+
|<float3>              |Coefficient of eigenvector in the vorticity direction.                               |
+----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The current best selection of coefficients is given by:

::

   Q Tensor Diffusivity = 0 1.0 1.0 0.5

-------------------------
**Technical Discussion**
-------------------------

The three directions (1, 2, 3) are often called the (flow, gradient, vorticity) directions.
Here, vorticity is *not* curl(u), but defined (along with the other three) for a particular set
of circumstances: steady simple shear flow. Their analogous definitions in other
regimes, as well as the selection of the coefficients, is an active area of research. The
interested reader should review the references listed below.



--------------
**References**
--------------

Brady, J. F. and Morris J. F., “Microstructure of strongly sheared suspensions and its
impact on rheology and diffusion,” J. of Fluid Mechanics, v. 348 pp.103-139, Oct 10,
1997.

Fang, Z. W., Mammoli, A. A., Brady, J.F., Ingber, M.S., Mondy, L.A. and Graham,
A.L., “Flow-aligned tensor models for suspension flows,” Int. J. of Multiphase Flow, v.
28(#1) pp. 137-166, January 2002.

Hopkins, M. M., Mondy, L. A., Rao, R. R., Altobelli, S. A., Fang, Z., Mammoli, A. A.
and Ingber, M. S., 2001. “Three-Dimensional Modeling of Suspension Flows with a
Flow-Aligned Tensor Model”, The 3rd Pacific Rim Conference on Rheology, July 8-
13, 2001, Vancouver, B.C., Canada.

Morris, J. F. and Boulay, F., “Curvilinear flows of noncolloidal suspensions: The role of
normal stresses,” J. of Rheology, v. 43(#5) pp. 1213-1237 Sep-Oct 1999.
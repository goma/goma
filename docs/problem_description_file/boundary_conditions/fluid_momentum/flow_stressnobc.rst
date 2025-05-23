*******************
**FLOW_STRESSNOBC**
*******************

::

	BC = FLOW_STRESSNOBC SS <bc_id> <float> [integer]

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition card applies the *free* outflow boundary condition developed
by Papanastasiou, et.al. (1992) on the fluid momentum with the option of setting the
pressure level. It is appropriate only for outflow boundaries where it is inappropriate to
use natural boundary conditions or *FLOW_PRESSURE*-class boundary conditions. It is
only supported for generalized Newtonian fluid constitutive equations.

Definitions of the input parameters are as follows:

+-------------------+------------------------------------------------------------+
|**FLOW_STRESSNOBC**| Name of the boundary condition.                            |
+-------------------+------------------------------------------------------------+
|**SS**             | Type of boundary condition (<bc_type>), where **SS**       | 
|                   | denotes side set in the EXODUS II database.                |
+-------------------+------------------------------------------------------------+
|<bc_id>            | The boundary flag identifier, an integer associated with   |
|                   | <bc_type> that identifies the boundary location (side set  |
|                   | in EXODUS II) in the problem domain.                       |
+-------------------+------------------------------------------------------------+
|<float>            | :math:`P_{applied}`, the applied pressure.                 |
+-------------------+------------------------------------------------------------+
|[integer]          | An optional parameter.                                     |
|                   |                                                            |
|                   |   * blank/-1 the pressure in the normal stress is replaced |
|                   |     by :math:`P_{applied}`.                                |
|                   |   * ≠ –1 the pressure in the solution vector is            |
|                   |     retained in the normal stress.                         |
+-------------------+------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:
::

     BC = FLOW_STRESSNOBC SS 10   1.0 -1

Here the boundary condition is applied to sideset 10 with a constant pressure of 1.0
required.

-------------------------
**Technical Discussion**
-------------------------

* The finite element formulation of the fluid momentum equations generates
  boundary integrals of the form:

.. figure:: /figures/096_goma_physics.png
	:align: center
	:width: 90%

where *P* is the isotropic pressure and :math:`\tau` the viscous stress. Often this boundary term
is left off entirely on a particular boundary with the result that a zero normal force
is applied implicitly. These are referred to as imposing a “natural” boundary
conditions. Alternatively, this integral might be included but with the integrand
replaced by an known value of force. This is the concept behind the
*FLOW_PRESSURE* and *FLOW_HYDROSTATIC* boundary conditions.

However, both types of boundary conditions imply that something is known about
the stress and, by association, the velocity field on the boundary. It is often the case
that outflow boundaries are present where it is difficult to provide this information.
A prime example is the outflow of a fluid jet accelerating downward due to
gravity. In this case, the downward velocity field is still developing at this
boundary so it is problematic to specify a stress value. Other examples include
imposing conditions at a “truncated” outflow where the exiting fluid is still
developing.

The *FLOW_STRESSNOBC* seeks to remedy this problem. Formulationally, the
boundary term as written above is included as just another term dependent upon
solution degrees of freedom. This permits the pressure and velocity gradients on
the boundary to float as needed so that one does not need to say anything about the
stress or pressure on the boundary.

Now strictly speaking, the ellipticity of the viscous flow equations suggests that
this operation should result in an ill-posed problems. Elliptic equations by their
very nature require that something be said about every boundary in the problem.
However, in the case of outflow boundaries it appears that this restriction can be
relaxed in certain circumstances with good results. Papanastasiou, et.al., (1997),
Renary (1997), Griffiths (1997) and Sani and Gresho (1994) discuss this.

* The boundary condition does permit that the pressure value be fixed while the
  viscous stress is allowed to float. This is done by setting the optional parameter to
  -1 and supplying the pressure value as :math:`P_{applied}`. When this is done depends upon
  circumstance. Note that this is distinctly different from setting a normal stress
  component using *FLOW_PRESSURE*.

* As noted above, this boundary condition is currently implemented only for
  generalized Newtonian fluid models. Polymeric fluid models will not work with it.



--------------
**References**
--------------

Griffiths, D.F., “The ‘no boundary condition’ outflow boundary condition,” IJNMF, 24,
393-411, (1997)

Papanastasiou, T. C., N. Malamataris, and K. Ellwood, “A New Outflow Boundary
Condition”, IJNMF, 14, 587-608, (1992).

Renardy, M., “Imposing ‘NO’ boundary conditions at outflow: Why does this work?”
IJNMF, 24, 413-417, (1997).

Sani, R.L., and P.M. Gresho, “Resume and remarks on the open boundary condition
minisymposium,” IJNMF, 18, 983-1008, (1994).

.. TODO - Line 62 contains a photo that needs to be exchanged for the equation.
**************
**CURRENT_BV**
**************

::

	BC = CURRENT_BV SS <bc_id> <integer> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(WIC/POTENTIAL)**

The *CURRENT_BV* card enables the specification of variable electrical current density
as given by Butler-Volmer kinetics and the Faraday’s law at the specified boundary
(namely, an electrode surface).

The <floatlist> has seven parameters for this boundary condition; definitions of the
input parameters are as follows:

=============== =================================================================
**CURRENT_BV**  Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS**
                denotes side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set
                in EXODUS II) in the problem domain.
<integer>       Species number of concentration.
<float1>        Stoichiometric coefficient.
<float2>        Kinetic rate constant.
<float3>        Reaction order.
<float4>        Anodic direction transfer coefficient.
<float5>        Cathodic direction transfer coefficient.
<float6>        Electrode potential or applied voltage.
<float7>        Theoretical open-circuit potential.
=============== =================================================================

------------
**Examples**
------------

An example input card:
::

   BC = CURRENT_BV SS 1 0   -1.0 0.000002 1.0 0.21 0.21 -0.65   -0.22

-------------------------
**Technical Discussion**
-------------------------

Users are referred to Chen (2000) for details of the Butler-Volmer model and also
Newman (1991), particularly Equations 8.6 and 8.10 and Chapter 8, pp. 188-189 in the
latter.



--------------
**References**
--------------

GTM-025.0: Modeling diffusion and migration transport of charged species in dilute
electrolyte solutions: GOMA implementation and sample computed predictions from a
case study of electroplating, K. S. Chen, September 21, 2000


J. S. Newman, "Electrochemical Systems", Second Edition, Prentice-Hall, Inc. (1991).



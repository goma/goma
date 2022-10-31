**************
FLOW_GRADV_T 
**************

::

	BC = FLOW_GRADV_T SS <bc_id> <float> [integer]

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition card stipulates a vanishing normal velocity gradient on a
boundary with the option of setting the pressure level.

Definitions of the input parameters are as follows:

FLOW_GRADV_T 
  Name of the boundary condition.

SS         
  Type of boundary condition (<bc_type>), where **SS** denotes    
  side set in the EXODUS II database.                             

<bc_id>
  The boundary flag identifier, an integer associated with        
  <bc_type> that identifies the boundary location (side set in    
  EXODUS II) in the problem domain.                               

<float>
  :math:`P_{applied}`, the applied pressure.                      

[integer]      
  An optional parameter.                                          

  *blank/-1*  the pressure in the normal stress is replaced      
  by :math:`P_{applied}`.                            

  *≠ –1*     the pressure in the solution vector is              
  retained in the normal stress.                     

------------
**Examples**
------------

The sample input card:
::

     BC = FLOW_GRADV_T SS 15   0.0

sets the transpose of gradient of velocity normal to sideset 15 to zero. A pressure value of zero is
used in the boundary condition.

::

    BC = FLOW_GRADV_T SS 15   0.0   1.0

In the preceding example, the pressure value used is obtained from the solution itself.

-------------------------
**Technical Discussion**
-------------------------

* This boundary condition is related in form and formulation to the
  *FLOW_STRESSNOBC* boundary condition in that it includes terms for the
  boundary integrals that appear in the momentum equation after application of
  integration by parts and the divergence theorem. In this boundary condition, the
  following integral is included with the momentum equation (for Newtonian):

  
.. math::

   n \cdot T = n \cdot (-pI + \mu \nabla v^T)

  where :math:`\mu` is the viscosity of a Newtonian or generalized Newtonian fluid. As in the
  case of the *FLOW_STRESSNOBC* condition the preceding integral appears as a
  function of pressure and velocity unknowns as any other term.
  
  This imposes a natural BC of :math:`n\cdot\nabla v = 0`

* The pressure term in the preceding may be replaced by a fixed, imposed pressure
  value. This is done by setting the optional input integer to -1 and providing the
  imposed value in :math:`P_{applied}` ; otherwise, the value set in :math:`P_{applied}` is ignored.

* Modifications are done to work for Viscoelastic flow with DEVSS-G stabilization




.. TODO - Line 68 contains a photo that needs to be exchanged for the equation.
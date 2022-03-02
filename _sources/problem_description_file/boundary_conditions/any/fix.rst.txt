*******
FIX
*******

::

	BC = FIX NS <bc_id> {char_string} <integer1>

-----------------------
Description / Usage
-----------------------

**(DC/VARIED)**

This boundary condition card is used to fix the value of a nodal variable along a node
set to the value it receives from an initial guess file (viz. either from the neutral file
specified by the *Initial Guess* card or an input EXODUS II file as also specified by the
*read_exoII_file* option on the *Initial Guess* card). The boundary condition is applied as
a Dirichlet condition (see technical discussion below).

Definitions of the input parameters are as follows:

FIX
    Name of the boundary condition (<bc_name>).                             

NS
    Type of boundary condition (<bc_type>), where **NS** denotes            
    node set in the EXODUS II database.                                     

<bc_id>          
    The boundary flag identifier, an integer associated with                
    <bc_type> that identifies the boundary location (node set in            
    EXODUS II) in the problem domain.                                       
{char_string}    
    Variable name that is to be fixed. This parameter can have              
    the following permissible values:                                       
                                                                                         
    * ``VELOCITY1``                                                       
    * ``VELOCITY2 VELOCITY3``                                            
    * ``MESH_DISPLACEMENT1``                                              
    * ``MESH_DISPLACEMENT2``                                              
    * ``MESH_DISPLACEMENT3``                                              
    * ``SOLID_DISPLACEMENT1``                                             
    * ``SOLID_DISPLACEMENT2``                                             
    * ``SOLID_DISPLACEMENT3 MASS_FRACTION``                              
    * ``TEMPERATURE PRESSURE VOLTAGE FILL``                            
    * ``POLYMER_STRESS11 POLYMER_STRESS12``                              
    * ``POLYMER_STRESS13 POLYMER_STRESS22``                              
    * ``POLYMER_STRESS23 POLYMER_STRESS33``                              
    * ``VELOCITY_GRADIENT11``                                             
    * ``VELOCITY_GRADIENT12``                                             
    * ``VELOCITY_GRADIENT13``                                             
    * ``VELOCITY_GRADIENT21``                                             
    * ``VELOCITY_GRADIENT22``                                             
    * ``VELOCITY_GRADIENT23``                                             
    * ``VELOCITY_GRADIENT31``                                             
    * ``VELOCITY_GRADIENT32``                                             
    * ``VELOCITY_GRADIENT33 POR_LIQ_PRES``                               
    * ``POR_GAS_PRES POR_POROSITY``                                      
    * ``POR_POROSITY POR_TEMP POR_SATURATION``                          
    * ``POR_LAST MAX_POROUS_NUM``                                        
    * ``POR_SINK_MASS VORT_DIR1 VORT_DIR2``                             
    * ``VORT_DIR3 VORT_LAMBDA CURVATURE``                               
    * ``LAGR_MULT1 LAGR_MULT2 LAGR_MULT3``                              
    * ``BOND_EVOLUTION SURF_CHARGE``                                     
    * ``EXT_VELOCITY EFIELD1 EFIELD2 EFIELD3``                         
    * ``ENORM NORMAL1 NORMAL2 NORMAL3``                                
    * ``SHELL_CURVATURE SHELL_TENSION``                                  
    * ``SHELL_X SHELL_Y SHELL_USER PHASE1``                            
    * ``PHASE2 PHASE3 PHASE4 PHASE5``                                  
    * ``SHELL_ANGLE1 SHELL_ANGLE2``                                      
    * ``SHELL_SURF_DIV_V SHELL_SURF_CURV``                               
    * ``N_DOT_CURL_V GRAD_S_V_DOT_N1``                                   
    * ``GRAD_S_V_DOT_N2 GRAD_S_V_DOT_N3``                                
    * ``ACOUS_PREAL ACOUS_PIMAG``                                        
    * ``SHELL_DIFF_FLUX SHELL_DIFF_CURVATURE``                           
    * ``SHELL_NORMAL1 SHELL_NORMAL2``                                    
    * ``ACOUS_REYN_STRESS SHELL_BDYVELO``                                
    * ``SHELL_LUBP LUBP SHELL_FILMP``                                   
    * ``SHELL_FILMH SHELL_PARTC``                                        
    * ``SHELL_SAT_CLOSED SHELL_PRESS_OPEN``                              
    * ``SHELL_TEMPERATURE SHELL_DELTAH``                                 
    * ``SHELL_LUB_CURV SHELL_SAT_GASN``                                  
    * ``SHELL_SHEAR_TOP SHELL_SHEAR_BOT``                                
    * ``SHELL_CROSS_SHEAR MAX_STRAIN``                                   
    * ``CUR_STRAIN LUBP_2 SHELL_PRESS_OPEN2``                           
    * ``SHELL_LUB_CURV_2``                                                 

<integer1>
    Species number of concentration, or zero if variable is not             
    concentration.                                                          

------------
Examples
------------

The following is an example of using this card to set the mesh displacement
components in a 2-D problem:
::

	BC =    FIX    NS    4    MESH_DISPLACEMENT1    0

::

	BC =    FIX    NS    4    MESH_DISPLACEMENT2    0

In this example, several continuation steps were taken to deform part of an elastic block
of material. The displacements on boundary node set 4 were then held constant while
moving another boundary (because the current displacements were not known, FIX
was a convenient tool).

-------------------------
Technical Discussion
-------------------------

This boundary condition capability is indispensable for moving-mesh problems when
the dependent variable is the mesh displacement from a stress free state. If one were to
try to use the DX/DY/DZ type Dirichlet condition to suddenly freeze a mesh along a
node set after a parameter continuation or transient problem restart, then they would be
faced with figuring out the displacement of each node and defining individual node sets
for each node for boundary condition application. This capability is also beneficial
when using previous simulation results to generate boundary conditions for more
complex analysis. We have on occasion used this boundary condition for most of the
variable types shown.





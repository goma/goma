************
moment
************

::

	EQ = moment{0|1|2|3} {Galerkin_wt} {MOM0|MOM1|MOM2|MOM3} {Interpol_fnc} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for one component of
a population balance equation

All four moments are expected to be enabled at the same time


moment0 | moment1 | moment2 | moment3
   Name of the equation to be solved, where the 
   and 3 components correspond to one of the       
   principal coordinate directions, e.g. X, Y and Z
   for Cartesian geometry.                         
{Galerkin_wt}                        
   Two- or four-character value that defines the    
                                     type of weighting function for this equation,    
                                     where:                                           
                                                                                      
                                      * **Q1**-Linear                                 
                                      * **Q2**-Quadratic                              
MOM0 | MOM1 | MOM2 | MOM3
   Name of the variable associated with the moment equation

{Interpol_fnc}                       
   Two- or four-character value that defines the    
                                     interpolation function used
                                                                                      
                                      * **Q1**-Linear Continuous                      
                                      * **Q2**-Quadratic Continuous                   
<float1>                             
   Multiplier on mass matrix term ( d ⁄dt ).        
<float2>                             
   Multiplier on advective term.                                   
<float4>                             
   Multiplier on diffusion term.                    
<float5>                             
   Multiplier on source term.                       
<float6>                             
   Multiplier on divergence term

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

::

          EQ = moment0 Q1 MOM0 Q1 1. 1. 1e-6  1. 1.
          EQ = moment1 Q1 MOM1 Q1 1. 1. 1e-6  1. 1.
          EQ = moment2 Q1 MOM2 Q1 1. 1. 1e-6  1. 1. 
          EQ = moment3 Q1 MOM3 Q1 1. 1. 1e-6  1. 1.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

Ortiz, Weston, Lisa Mondy, Christine Roberts, and Rekha Rao. "Population balance modeling of polyurethane foam formation with pressure‐dependent growth kernel." AIChE Journal 68, no. 3 (2022): e17529.
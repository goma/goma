**************
**shear_rate**
**************

::

	EQ = eddy_visc {Galerkin_wt} eddy_nu {Interpol_fnc} <float1> <float2> <float3> <float4> <float5>

-----------------------
**Description / Usage**
-----------------------

This card provides information for solving a differential equation for Spalart
Allmaras turbulence model Definitions of the input parameters are defined below.
Note that <float1> through <float5> define the constant multipliers in front of
each type of term in the equation. The Galerkin weight and the interpolation
function must be the same for the code to work properly.

eddy_visc
    Name of the equation to be solved

{Galerkin_wt}
   Two-character value that defines the type of weighting function for this equation,
   where:
   Q1-Linear
   Q2-Quadratic

eddy_nu
    Name of the variable associated with this equation

{Interpol_fnc}
   Two-character value that defines the interpolation function used to represent the
   variable eddy_nu, where:
   Q1-Linear
   Q2-Quadratic

<float1>
   Multiplier on the mass term

<float2>
   Multiplier on the advective term

<float3>
   Multiplier on the boundary terms

<float4>
   Multiplier on the diffusion terms

<float5> 
   Multiplier on the source term


Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The following is a sample card that uses quadratic continuous interpolation for the
species equation and turns on all the term multipliers:
::

   EQ = eddy_visc Q2 eddy_nu Q2 1. 1. 1. 1. 1.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.


--------------
**References**
--------------

No References.
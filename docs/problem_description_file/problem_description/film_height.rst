***********
film_height
***********

::

	EQ = film_height {Galerkin_wt} FILM_HEIGHT {Interpol_fnc} <floatlist>

-------------------
Description / Usage
-------------------

This is the equation card for the film height 2.5D equations for film casting problems.
It is used in conjunction with momentum1 and momentum2 cards as well as mesh1 and mesh2 cards
to solve film casting problems.

film_height
    Name of the equation
{Galerkin_wt}
    Q1 Linear
    Q2 Quadratic
FILM_HEIGHT
    Variable associated with the film_height equation
{Interpol_fnc}
    Q1 Linear
    Q2 Quadratic
<float1>
    Mass term
<float2>
    advection term
<float3>
    boundary term
<float4>
    diffusion term
<float5>
    source term

Note: These multipliers are intended to provide a means of activating or deactivating
terms of an equation, and hence should be set to zero or one. If a multiplier is zero, the
section of code that evaluates the corresponding term will be skipped.

------------
**Examples**
------------

The literature suggests treating the film_height equation like pressure and using Taylor-Hood like elements for 
velocity and film_height combinations.

::

   EQ = momentum1 Q2 U1 Q2 0. 1. 1. 1. 1. 0.
   EQ = momentum1 Q2 U1 Q2 0. 1. 1. 1. 1. 0.
   EQ = film_height Q1 FILM_HEIGHT Q1 1. 1. 1. 1. 1.
   EQ = mesh1 Q2 D1 Q2 0. 0. 1. 1. 1.
   EQ = mesh2 Q2 D2 Q2 0. 0. 1. 1. 1.

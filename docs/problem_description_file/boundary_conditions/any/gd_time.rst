***********
GD_TIME
***********

::

	BC = BC = GD_TIME SS <bc_id> <equation_name> <integer1> <time_func_name> <integer2> <float1> <float2> [float3]

-----------------------
Description / Usage
-----------------------

**(PCC/VARIED)**

This boundary condition card is actually a multiplicative building block that can be
used to impose a multiplicative time modulation of a specified functional form on any
set of *GD_** conditions. *NOTE: unlike the other* *GD_** *cards which are additive*, 
*this card is multiplicative*. This condition must be placed after any single or set of *GD_**
cards for which the user wishes to modulate (viz. *GD_LINEAR, GD_PARAB*, etc.). The
card can be used as many times as needed to construct the desired function. The
examples below will clarify its use. Definitions of the input parameters are as follows:

GD_TIME
    Name of the boundary condition (<bc_name>).
SS
    Type of boundary condition (<bc_type>), where **SS** denotes side set in
    the EXODUS II database.
<bc_id> 
    The boundary flag identifier, an integer associated with <bc_type> that
    identifies the boundary location (side set in EXODUS II) in the problem
    domain.
<equation_name>
    A character string indicating the equation to which this boundary condition
    is applied (see the list of permissible values in the discussion above for
    *Category 1*).
<integer1>
    Species number of the mass transport equation. The value should be 0 unless
    the <equation_name> is of type R_MASS.
<time_func_name>
    Keyword to identify the functional form of the time modulation. Permissible
    values for this parameter are LINEAR, EXPONENTIAL, and SINUSOIDAL.
<integer2>
    Set this required but unused parameter to zero.
<float1>
    :math:`C_0` model parameter
<float2>
    :math:`C_1` model parameter
[float3]
    Optional parameter to add a maximum time to be applied if :math:`t > t_{max}` then
    :math:`t` is set to :math:`t_{max}`
    
The functional form of each time-modulation model is as follows:

.. tabularcolumns:: |l|L|

===================== ===============================================================
LINEAR:               :math:`f(t) = C_0 + C_1t`
EXPONENTIAL:          :math:`f(t) = \exp \left(C_0 + C_1t \right)`
SINUSOIDAL:           :math:`f(t) = \sin \left(C_0 + C_1t \right)`  
===================== ===============================================================

------------
Examples
------------

Following is a sample card set:
::

	BC = GD_LINEAR SS 1 R_MESH_NORMAL 0 MESH_DISPLACEMENT1 0 1. 0.
	BC = GD_TIME SS 1 R_MESH_NORMAL 0 SINUSOIDAL 0 10. 2.
	BC = GD_LINEAR SS 1 R_MESH_NORMAL 0 MESH_POSITION1 0 0. -1.

This set of cards leads to the application of :math:`x = sin(10.0 + 2t)` to the normal component
of the mesh displacement at side set 1. If side set 1 were a surface of constant x (viz.
normal in the x-direction) then this condition could be used to impose a piston motion
to the surface. Recall that *GD_LINEAR* cards are additive with each other and
*GD_TIME* is multiplicative with the previous cards. The first card is used to put a
constant of 1.0 in the equation, the second card (*GD_TIME* card) multiplies that
constant with the sinusoidal time function, and the third card is used to put the linear
term on mesh position. Note carefully the signs used.

Invoking with a maximum time is done using the optional parameter:

::

    BC = GD_LINEAR SS 1 R_MOMENTUM1 0 VELOCITY1 0 {web_speed/time_max} 0.
    BC = GD_TIME SS 1 R_MOMENTUM1 0 LINEAR 0 0 1. {time_max}
    BC = GD_LINEAR SS 1 R_MOMENTUM1 0 VELOCITY1 0 0. -1.

-------------------------
**Technical Discussion**
-------------------------

This boundary condition building block is very useful for imposing time-dependent
boundary conditions with some fairly standard functional forms without the
inconvenience of writing a user-defined boundary condition. Boundary conditions for
pulsating flow, piston motion, roll-eccentricity effects in coating, time-evolving
temperature transients, etc. can all be constructed using this card. The examples at the
end of this section on *GD_** options will help the user construct such functions.

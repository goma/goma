********************************
AC (Phase Velocity of Level Set)
********************************

::

    AC = LSV <integer_list> <VX|VY|VZ> <float1>

-----------------------
Description / Usage
-----------------------

This type of augmenting condition is used to connect a boundary condition to an 
augmenting constraint on the average velocity of one of the phases of the level set field. 
This is useful for performing a simulation in the reference frame of a moving material 
when the velocity is not known a priori, such as a bubble rising through a quiescent 
fluid.

The <integer_list> has four parameters; definitions of the input parameters are as 
follows:

LSV | LS_VELOCITY
    A mandatory string indicating the name of the augmenting 
    condition card.

<integer1>
    An integer parameter giving the material identification 
    number over which the phase volume constraint is applied. 
    That is to say, if more than one material is present, the phase
    velocity constraint can only be applied to the space 
    occupied by the level set phase in one material.

<integer2>
    An integer parameter giving the index of the boundary
    condition card whose parameter is being used as the new 
    degree of freedom. Numbering begins with zero, starting 
    with the first boundary condition in the input file. For this 
    augmenting condition, this is generally a velocity boundary 
    condition.

<integer3>
    An integer parameter that identifies the floating point 
    parameter of the boundary condition that will be varied. It is 
    an index that starts at zero with the left-most float value on 
    the <int2> boundary condition cast and increments upward 
    from left to right.

<integer4>
    An integer parameter which specifies the level set phase 
    whose velocity will be constrained by this augmenting 
    condition. 

    - > 0: The constraint will be placed on the phase where the level set function has positive values.
    - < 0: The constraint will be placed on the phase where the level set function has negative phase.

<VX|VY|VZ>
    This string parameter specifies which component of the 
    velocity will be constrained by the augmenting condition.

    - VX: The x component of the phase velocity will be constrained.
    - VY: The y component of the phase velocity will be constrained.
    - VZ: The z component of the phase velocity will be constrained.

<float1>
    A float parameter that fixes the volume averaged velocity of 
    the level set phase specified by <int4>.

------------
Examples
------------

The following is an example of using this augmenting condition to perform a 
calculation in the moving reference frame of the negative level set phase where the 
reference frame velocity is not known a priori:

::

    AC = LSV 1 0 0 -1 VX 0.0

Applied to boundary condition:

::

    BC = U NS 4 0.0 1.0

The x component of the velocity of the negative level set phase of material #1 is set to 
zero and the first floating point parameter of the first velocity boundary condition is 
used as the new degree of freedom. Note that the Dirichlet boundary condition has the 
additional float parameter 1.0 so that the boundary condition equation is retained as a 
residual equation. The initial guess for the boundary velocity in this case is 0.

-------------------------
Technical Discussion
-------------------------

The velocity of the specified phase of the level set field is calculated using a volume 
average over the domain where the level set fill function has positive or negative values 
as specified. The boundary between the two phases is defined as the zero contour line 
of that fill function. The 'thickness' of the region around the embedded interface where 
material properties transition between the properties of the two pure phases is specified 
by the "Level Set Length Scale" card in the input deck. Within the diffuse interface 
region, the contribution to the averaged velocity integral is weighted by a smoothed 
Dirac delta function. Additional information about the embedded interface and the 
form of the delta function are available in the technical discussion section of the "Level 
Set Length Scale" card.

See the technical discussion section of the AC (Boundary Condition) card for 
additional information regarding augmenting conditions.

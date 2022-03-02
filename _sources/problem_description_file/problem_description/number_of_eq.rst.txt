****************
**Number of EQ**
****************

::

	Number of EQ = <integer>

-----------------------
**Description / Usage**
-----------------------

This card is required for each material section in the *Problem Description File*. It
specifies how many equations (i.e., equation cards, *[EQ =]*) follow for this material
section, including the mesh motion equations if appropriate. This number of equations
is only for the current material, since each material has its own equation section.

The single input parameter is defined as

========= ====================================================================
<integer> The number of EQ cards following this card. Only the first
          Number of EQ equations are read; if there are more EQ
          cards than specified by <integer>, *Goma* ignores the extras.
          If <integer> is set to -1, *Goma* will automatically count the
          number of EQ cards between the *Number of EQ* card and
          the *END OF EQ* card.
========= ====================================================================

------------
**Examples**
------------

The following is a sample card that sets the number of equations to 5:
::

   Number of EQ = 5

-------------------------
**Technical Discussion**
-------------------------

For equation specification in *Goma*, it is important to remember that a scalar equation
has a single equation entry (e.g. fill, species, voltage, shear rate, etc.), while a vector equation (e.g. momentum, mesh, mom_solid, etc.) has an entry for each component of
the vector. Thus, if you were solving a two-dimension flow problem, you would need
to specify both U1 and U2 components of the momentum equation explicitly. The same
holds true for tensor equations (e.g. stress and velocity gradient); each term of the
tensor is specified explicitly. The one exception to this rule is for multimode
viscoelasticity where the first mode equations are specified through the equation card
and then the auxiliary modes are set by the *Number of viscoelastic modes* card. Please
see the viscoelastic tutorial memo (Rao, 2000) for a detailed discussion of multimode
viscoelasticity.



--------------
**References**
--------------

GT-014.1: Tutorial for Running Viscoelastic Flow Problems with GOMA, June 21,
2000, R. R. Rao

______________________________________________________________________

Equation Cards

Following the *Number of EQ card, the equation cards, or records, are racked as intended up to the *END OF EQ* card or to the number specified, with one equation record per line. Each card begins with the “*EQ* =” string, followed by the equation name, e.g., *energy*, some basis function and trial function information, and finally a series of term multipliers. These multipliers are intended to provide a means of activating or deactivating terms of an equation, and hence should be set to zero or one. However, one can use these multipliers as a way of adjusting the scaling of individual terms. Exercise caution in using these factors as expedients for transport coefficients; for instance
the equation term multiplier for the momentum diffusion term affects both the isotropic stress term (pressure) and the deviatoric stress. It is recommended that you consult the example tutorial menus and problems to get a feel for the structure of this section. A sample input file structure including the EQ section is shown in the figure at the beginning of this chapter.
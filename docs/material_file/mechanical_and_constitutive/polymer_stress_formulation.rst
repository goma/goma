**************************
Polymer Stress Formulation
**************************

::

   Polymer Constitutive Equation = {model_name}

-----------------------
**Description / Usage**
-----------------------

This card specifies which formulation of the polymer constitutive equation should be
used. Valid options are


**EVSS_G**       
   Uses the classic elastic-viscous stress splitting of Rajagopalan (1990) where the stress is the elastic     
   stress only without a Newtonian component. This option is the default if this *Polymer Stress Formulation*  
   card is not supplied. This formulation is almost never used. Prefer EVSS_F

**EVSS_F** 
   Uses the EVSS formulation of Guenette and Fortin (1995) that solves the standard stress equation with the   
   addition of a new term to the momentum equation. This formulation is used most often.                       

**EVSS_L**       
   Uses a research formulation for viscoelasticity that includes a level set discretization that switches the  
   equations from solid to fluid. This option is not currently in production usage. Partial level set support is
   included in EVSS_F formulation

**LOG_CONF**       
   Log-conformation tensor formulation from Fattal and Kupferman 2004, uses DEVSS-G

**LOG_CONF_GRADV**       
   Log-conformation tensor formulation from Fattal and Kupferman 2004, uses DEVSS-G but all gradient terms in constitutive
   equation are the field variable :math:`\nabla v` instead of the projection :math:`G`

**SQRT_CONF**       
   sqrt-conformation tensor formulation from Balci et al. 2011, uses DEVSS-G stabilization

------------
**Examples**
------------

The following is a sample card that sets the polymer stress formulation to EVSS_F:

::

   Polymer Stress Formulation = EVSS_F

-------------------------
**Technical Discussion**
-------------------------

If using *SQRT_CONF* with no guess for the square root of stress tensor, :math:`b`,
recommended initial guess is the identity tensor for all modes.

Use post processing card *Map Conf Stress* to output the stress values, otherwise the usual S values are the given conformation tensor 
base form such as the SQRT being :math:`b` in :math:`b^Tb = c` or LOG being :math:`s = log c` 


--------------
**References**
--------------

Guenette, R. and M. Fortin, “A New Mixed Finite Element Method for Computing
Viscoelastic Flow,” J. Non-Newtonian Fluid Mech., 60 (1995) 27-52.

Rajagopalan, D., R. C. Armstrong and R. A. Brown, “Finite Element Methods for
Calculation of Viscoelastic Fluids with a Newtonian Viscosity”, J. Non-Newtonian
Fluid Mech., 36 (1990) 159-192.
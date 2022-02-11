*****************************
Polymer Constitutive Equation
*****************************

::

   Polymer Constitutive Equation = {model_name}

-------------------
Description / Usage
-------------------

This required card is used to specify the polymer constitutive equation. A single input
parameter must be defined, that being the {model_name}.

{model_name}
    Name of the constitutive equation model, being one of the following values: **NOPOLYMER, OLDROYDB, GIESEKUS,
    PTT, SARAMITO_OLDROYDB, SARAMITO_GIESEKUS, SARAMITO_PTT, WHITE-METZNER**. Several of these polymer          
    constitutive models require additional parameters for the polymer properties that are entered via           
    additional cards, as described below. Please see the Example sectionand the tutorial referenced below.      

Thus,

NOPOLYMER        
    For Newtonian and generalized Newtonian models. No floating point values are required.                      
OLDROYDB         
    For the Oldroyd-B constitutive model. This option requires four floating point values, which are described  
    below.                                                                                                      
GIESEKUS
    For the Giesekus model. This option requires five floating point values, which are described below.         
PTT
    For the Phan-Thien Tanner model. This option requires six floating point values, which are described below. 
SARAMITO_OLDROYDB
    For the Oldroyd-B model used with the Saramito yield model. This option requires six floating point values, 
    which are described below.                                                                                  
SARAMITO_GIESEKUS
    For the Giesekus model used with the Saramito yield model. This option requires seven floating point values,
    which are described below.                                                                                  
SARAMITO_PTT
    For the Giesekus model used with the Saramito yield model. This option requires eight floating point values,
    which are described below.                                                                                  
WHITE_METZNER    
    For the White-Metzner model. This option is not currently working.                                          

--------
Examples
--------

The following is a sample card that sets the polymer constitutive equation to
**NOPOLYMER**. This option does not require any additional cards since it indicates
that there is no polymer constitutive equation present.

::

   Polymer Constitutive Equation = NOPOLYMER

The following is a sample card that sets the polymer constitutive equation to
**OLDROYDB**. This option requires four cards describing the polymer stress
formulation, weight function, viscosity and time constant.

::

   Polymer Constitutive Equation = OLDROYDB

::

   Polymer Stress Formulation = EVSS_F

::

   Polymer Weight Function = GALERKIN

::

   Polymer Viscosity = CONSTANT 1.

::

   Polymer Time Constant = CONSTANT 0.02

The following is a sample card that sets the polymer constitutive equation to
**GIESEKUS**. This option requires five cards describing the polymer stress formulation,
weight function, viscosity, time constant and mobility parameter.

::

   Polymer Constitutive Equation = GIESEKUS

::

   Polymer Stress Formulation = EVSS_F

::

   Polymer Weight Function = GALERKIN

::

   Polymer Viscosity = CONSTANT 1.

::

   Polymer Time Constant = CONSTANT 0.2

::

   Mobility Parameter = CONSTANT 0.1

The following is a sample card that sets the polymer constitutive equation to 
**PHANTHIEN TANNER** (or **PTT**). This option requires six additional cards that set the
polymer stress formulation, weight function for the stress equation, viscosity, time
constant and nonlinear PTT parameters.:

::

   Polymer Consitutive Equation = PTT

::

   Polymer Stress Formulation = EVSS_F

::

   Polymer Weight Function = GALERKIN

::

   Polymer Viscosity = CONSTANT 8000.

::

   Polymer Time Constant = CONSTANT 0.01

::

   PTT Xi parameter = CONSTANT 0.10

::

   PTT Epsilon parameter = CONSTANT 0.05

The following is a sample card that sets the polymer constitutive equation to
**SARAMITO_OLDROYDB**. This option requires six cards describing the polymer
stress formulation, weight function, viscosity, time constant, yield stress, 
and yield exponent.

::

   Polymer Constitutive Equation = SARAMITO_OLDROYDB

::

   Polymer Stress Formulation = EVSS_F

::

   Polymer Weight Function = GALERKIN

::

   Polymer Viscosity = CONSTANT 1.

::

   Polymer Time Constant = CONSTANT 0.02

::

  Polymer Yield Stress = CONSTANT 15.

::

  Yield Exponent = CONSTANT 0.

The following is a sample card that sets the polymer constitutive equation to
**SARAMITO_GIESEKUS**. This option requires seven cards describing the polymer stress
formulation, weight function, viscosity, time constant, mobility parameter, yield 
stress, and yield exponent.

::

   Polymer Constitutive Equation = SARAMITO_GIESEKUS

::

   Polymer Stress Formulation = EVSS_F

::

   Polymer Weight Function = GALERKIN

::

   Polymer Viscosity = CONSTANT 1.

::

   Polymer Time Constant = CONSTANT 0.2

::

  Polymer Yield Stress = CONSTANT 12.

::

  Yield Exponent = CONSTANT 1.0

::

   Mobility Parameter = CONSTANT 0.1

The following is a sample card that sets the polymer constitutive equation to 
**SARAMITO_PTT**. This option requires eight additional cards that set the
polymer stress formulation, weight function for the stress equation, viscosity, time
constant, nonlinear PTT parameters, yield stress, and yield exponent.

::

   Polymer Consitutive Equation = SARAMITO_PTT

::

   Polymer Stress Formulation = EVSS_F

::

   Polymer Weight Function = GALERKIN

::

   Polymer Viscosity = CONSTANT 8000.

::

   Polymer Time Constant = CONSTANT 0.01

::

  Polymer Yield Stress = CONSTANT 200.

::

  Yield Exponent = CONSTANT 0.5

::

   PTT Xi parameter = CONSTANT 0.10

::

   PTT Epsilon parameter = CONSTANT 0.05

The following is a sample card that sets the polymer constitutive equation to
**WHITE_METZNER**. This option is not currently functional for multimode
viscoelasticity. If needed it could be resurrected with only minimal changes to the input
parser.

::

   Polymer Consitutive Equation = WHITE_METZNER

--------------------
Technical Discussion
--------------------

The viscoelastic tutorial is helpful for usage issues such as extensions from single mode
to multimodes.



----------
References
----------

GT-014.1: Tutorial for Running Viscoelastic Flow Problems with GOMA, June 21,
2000, R. R. Rao


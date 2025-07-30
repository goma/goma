
Continuation Specifications
###############################

This section of input records is used to direct all automatic continuation procedures. The entire
section is completely optional. Basically, automatic continuation can be accomplished in steady
state simulations (see *Time Integration* card) through any one or combination of parameters.
These parameters can be any one or combination of the input floats required on the boundary
condition cards (see Section 4.10) or material property cards (see Chapter 5). The cards in this
section are used to specify the parameters that will be marched automatically, the method of
marching (e.g. zero-order, first-order, multiparameter first-order, etc.), the limits of parameter
values, and other sundry options. Much of this capability can now be managed from the LOCA
library package (Library of Continuation Algorithms - Salinger et al. 2002).

Input specifications for this section of input records is discussed in a separate, comprehensive
manual (Gates, et. al., 2000); an update to this manual has been completed during the summer of
2006 (Labreche, et. al., 2006).

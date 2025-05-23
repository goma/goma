**********
**LS_ADC**
**********

::

	BC = LS_ADC SS <bc_id> <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

**(Special/LEVEL SET)**

This boundary condition is used exclusively with level set interface tracking. It is used
to simulate contact and dewetting events. It employs a probabilistic model applied to
elements on a boundary that contain an interface to determine whether contact or
dewetting occurs there. It then uses a direct, brute force algortihm to manipulate the
level set field to enforce contact or dewetting.

A description of the input parameters follows:

============= ================================================================
**LS_ADC**    Name of the boundary condition.
**SS**        This string indicates that this boundary is applied to a
              sideset.
<bc_id>       This is a sideset id where contact or dewetting processes are
              anticipated. Only elements that border this sideset will be
              considered as possibilities for ADC events.
<float1>      :math:`\theta_c`, the capture angle in degrees.
<float2>      :math:`\alpha_c`, the capture distance (L)
<float3>      :math:`N_c`, the capture rate ( 1/ :math:`L^2` -T)
============= ================================================================

------------
**Examples**
------------

An example:
::

   BC = LS_ADC SS 10 15.0   0.2   100.0


-------------------------
**Technical Discussion**
-------------------------

It has been found that level set interface tracking problems that involve contact or
dewetting of the interfacial representation pose special problems for our numerical
method. To a certain extent, we can model this type of event by making special
modifications to the slipping properties of the boundary in question, however, this does
not always work, especially in the case of dewetting events.

What seems to be the trouble is that we are attempting to use continuum-based models
to simulate phenomena that essentially are due to molecular forces being expressed
over non-molecular length scales. These length scales, while big with respect to
molecules, are small with respect to our problem size. Hence, they are difficult to
include in the context of reasonable mesh spacing.

The approach this boundary condition takes to inclusion of contact and dewetting
phenomena is not attempt to model the finer details, but to simply note that they are due
to “molecular weirdness” and thus take place outside of ordinary continuum
mechanics. Therefore, there is some justification for, very briefly and in a localized
area, dispensing with continuum mechanics assumption and simply imposing a contact
or dewetting event. We refer to these as ADC events and will describe them in more
detail later.

The parameters supplied with the card are used to determine where and when such an
ADC event occurs. We have chosen to introduce a probabalistic model for this
purpose. The reasoning for this comes from reflecting on the dewetting problem. If
one imagines a thin sheet of fluid on a wetting substrate, it is clear that dewetting will
occur eventually at some point on that sheet. Where that event occurs is somewhat
random for a detached perspective. Introduction of a probability model for ADC
events attempts to capture this.

Whether an ADC event occurs at an element on the sideset is determined by the
following requirements:

* The interface surface passes through the element.

* There isn’t a contact line in the element.

* The angle between the interface normal and the sideset surface normal is less than
  or equal to the capture angle, :math:`\theta_c`.

* A random number in the range (0,1) determined by the standard C rand()
  function is less than a probability, *P*, given by

.. figure:: /figures/195_goma_physics.png
	:align: center
	:width: 90%

where *d* is the average distance of the interface to the sideset in that element, Dt is the
time step size, and *h* is the side length of the element (Note for 2D problems :math:`h^2` is
replaced by *h* where the other dimension is assumed unity in the *z* direction).

Interpretation of this probability relation might take the following course. Given that
the fluid interface lies within :math:`\alpha)c` of the surface, the length of time necessary before and
ADC event is certain to occur is given by 1/ :math:`N_ch^2` . Hence, the bigger the capture rate
parameter the faster this is likely to occur. The functional form for the case of *d* > 
:math:`\alpha_c` is
included merely to ensure that the probability drops smoothly to zero as quickly as
possible. One might point out that the probability at a specific element tends towards
zero as the element size decreases. Of course, in that context, the number of elements
should increase in number so that the overall probability of an ADC event should not
be a function of the degree of mesh refinement. A second point is that this boundary
condition can be made to function as means to initiate contact without delay by simply
choosing a capture rate that is large enought with respect to the current time step.

Application of an ADC event in a element that meets the preceding criteria is illustrated
in the cartoon below:

.. figure:: /figures/196_goma_physics.png
	:align: center
	:width: 90%

It is a simple manipulation of the level set values in that element so that the interface
will follow the path indicated by the dashed curve in the lower figure. No effort is
made in preservation of volume when this is done. The assumption is that these events
will occur infrequently enough that this is not a significant problem. However, the user
should be aware of this assumption and be careful that these events do not occur on a
regular basis as then the mass loss might be more significant.



--------------
**References**
--------------

No References. 

..
	TODO - Line 90 has a picture that needs to be exhanged with an equation.
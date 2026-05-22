***********************
Elliptic fxi/geta/hzeta
***********************

::

   Elliptic fxi = {model_name} {float_list}
   Elliptic geta = {model_name} {float_list} 
   Elliptic hzeta = {model_name} {float_list} 

-----------------------
**Description / Usage**
-----------------------

This is used to control mesh spacing in the elliptic mesh generation method.

Definitions of the input parameters are as follows:

{model_name}
  CONSTANT or SIMPLE_ABS
{float_list}
  See discussion below for the required number of floating point numbers and their meaning for each model.

CONSTANT
  takes one floating point number, :math:`val`, and returns that value for all inputs.

SIMPLE_ABS
  takes six floating point numbers, :math:`x_0`, :math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4` and :math:`a_5`

  :math:`val = a_5 + \frac{a_1}{a_2 + a_3 * |x - x_0|} + a_4 * |x - x_0|}`

DUAL_ABS
  takes five floating point numbers, :math:`x_0`, :math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`

  :math:`val = a_4 + \frac{a_1}{a_2 + a_3 * (|x| - |x_0|)}`

------------
**Examples**
------------

The following is a sample card:

::

   Elliptic geta = SIMPLE_ABS 0. 10.0 3.0 -1.0 0. 0.

-------------------------
**Technical Discussion**
-------------------------

--------
**FAQs**
--------
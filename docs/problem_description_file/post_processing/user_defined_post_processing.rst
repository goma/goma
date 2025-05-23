********************************
**User-Defined Post Processing**
********************************

::

   User-Defined Post Processing = {yes | no} <float_list>

-----------------------
**Description / Usage**
-----------------------

This option enables user-defined postprocessing options in *Goma*. An arbitrary number
of floating point constants can be loaded to use in the user-defined subroutine
user_post (user_post.c). This variable is called **USER** in the output EXODUS
II file and can be contoured or processed just like any other nodal variable in a postprocessing
visualization package.

============= ================================================================
**yes**       Calculate and write the user-defined postprocessing variable
              to the output EXODUSII file.
**no**        Do not calculate the user-defined postprocessing.
<float_list>  An arbitrary number (including zero) of floating point
              numbers, which can be accessed in file user_post
============= ================================================================

------------
**Examples**
------------

Consider the following sample input card:
::

   User-Defined Post Processing = yes 100.

Suppose you would like to contour the speed of a fluid in a two-dimensional problem
using this card, with your intent being to multiply the calculated value by a factor of
100.0 for some unit conversion or something. You must add

::

   post_value = param[0]*sqrt(fv->v[0]*fv->v[0] + fv->v[1]*fv->*v[1]) ;

to user_post.c. Note also that you have to comment out the error handler line just
above the location you enter the post_value code. The comments in the routine help
guide you through the process.

-------------------------
**Technical Discussion**
-------------------------

See the function user_post in user_post.c.




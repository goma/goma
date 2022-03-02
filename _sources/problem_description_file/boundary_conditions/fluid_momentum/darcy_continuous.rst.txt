********************
**DARCY_CONTINUOUS**
********************

::

	BC = DARCY_CONTINOUS SS <bc_id> <integer1> <integer2> [float1]

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MOMENTUM)**

This condition enforces continuity of mass flux at an interface between a continuous
medium and a saturated or partially saturated porous medium. In other words,
*DARCY_CONTINUOUS* is a boundary condition that equates the velocity component
in the liquid phase normal to the interface with the Darcy velocity in the porous phase,
normal to the same interface, with proper accounting for conservation of mass using
the liquid phase densities in the material files.

=================== ==========================================================
**DARCY_CONTINOUS** Name of the boundary condition (<bc_name>).
**SS**              Type of boundary condition (<bc_type>), where **SS**
                    denotes side set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (side set
                    in EXODUS II) in the problem domain.
<integer1>          Element block ID of porous phase from the EXODUS II  
                    database.
<integer2>          Element block ID of continuous fluid phase from the
                    EXODUS II database.
[float1]            An optional floating point that is used for level-set free
                    surface problems. This floating point represents a length
                    scale over which “contact” of a liquid free surface
                    represented by a level set field and a porous medium. It
                    should be set to some small integer multiple of the
                    smallest element size along the boundary. Note that this
                    length scale is only required in cases where “sharp”
                    interfaces using subelement integration are used. It is
                    not required for diffuse interface representations.
=================== ==========================================================

------------
**Examples**
------------

The boundary condition
::

    BC = DARCY_CONTINUOUS SS 5 2 1

applies to the interface defined by side set 5 which joins EXODUS II block 2 (porous
phase) and block 1 (continuous phase).

-------------------------
**Technical Discussion**
-------------------------

The *DARCY_CONTINUOUS* boundary condition imposes the following requirement at
the interface between a continuous medium and a saturated or partially saturated
porous medium:

.. figure:: /figures/115_goma_physics.png
	:align: center
	:width: 90%

where :math:`\underline{n}` is the outward-pointing normal to the surface, 
:math:`\underline{q}` is the Darcy flux, :math:`p_l` is the
liquid density, presumed to be the same in the adjacent phases, :math:`\underline{v}` 
is the fluid velocity and :math:`\underline{v}_s` is the mesh velocity.

Typically this boundary condition is applied between two blocks, one being of a
*LAGRANGIAN* mesh motion type (see *Mesh Motion* card) and the other being of an
*ARBITRARY* mesh motion type. Within the *LAGRANGIAN* material the *Media Type*
card is set to *POROUS_SATURATED, POROUS_UNSATURATED*, or
*POROUS_TWO_PHASE*. The other block is of type *CONTINOUS*.

Refer to the citations below where this boundary condition is discussed in more detail.


--------
**FAQs**
--------

**Important troubleshooting note**: Density, as specified in the material files for the
continuous and porous phase, MUST be the same for this boundary condition to make
sense.

--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

GT-028.0: Liquid Drop Impact on a Porous Substrate: a level-set tutorial, August 15,
2005.

.. TODO - Line 64 have photos that needs to be replaced with the real equation.
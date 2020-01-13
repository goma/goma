/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/
 

/*
 *$Id: el_elm_info.c,v 5.6 2010-03-03 22:33:57 prschun Exp $
 */

#include <string.h>
#include <stdio.h>

#include "std.h"
#include "el_elm.h"
#include "mm_as_const.h"
#include "mm_eh.h"
#include "el_elm_info.h"

#define GOMA_EL_ELM_INFO_C


/*************** R O U T I N E S   I N   T H I S   F I L E ********************
*
*  NAME                         TYPE            CALL_BY
* ---------------               -------         ------------------------
*  elem_info ()                 int             "mm_fill.c" matrix_fill 
*  find_stu  ()                 void            "mm_fill.c" matrix_fill
*  find_surf_st  ()             void            "mm_fill.c" matrix_fill
*  Gq_weight ()                 double          "mm_fill.c" matrix_fill
*  Gq_surf_weight ()            double          "mm_fill.c" matrix_fill
*  in_list   ()                 int              multiple routines
*  get_type  ()                 int
*
*
******************************************************************************/

int
elem_info(const int info,
          const int ielem_type )
     /*
      *        Function which returns the various parameters
      *       for the elements, e.g., polynomial order, number of
      *       Gaussian-quadrature points, etc, based upon what a code
      *       passed from the calling routine.
      *
      *        Author:          Scott Hutchinson (1421)
      *        Date:            15 May 1992
      *        Revised:         8 December 1992
      *        Revised:         1 December 2011 - SAR - TRISHELL
      *
      *   The routine currently handles the following requests for information
      *   about element type, ielem_type:
      *
      *       NNODES          Number of nodes in the element
      *       NQUAD           Number of volume quadrature points
      *       NDIM            Dimension of the element
      *       NQUAD_SURF      Number of surface quadrature points
      */
{
  int answer;
  
  /* return desired element information */
  
  switch( ielem_type ){           /* select type of element */

  case LINEAR_TRI:                /* linear triangle */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 3;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 3;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:              /* number of surface quad points */
      answer = 2;
      break;
    case NQUAD_EDGE:              /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

    case QUAD_TRI:                   /* quadratic triangle*/
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 6;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 6;
      /*
      answer = 4;
      */
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 3;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

    case QUAD6_TRI:               /* quadratic triangle with 6th order quadrature */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 6;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 12;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 5;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case P0_QUAD:                   /* constant on quadrilateral */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 1;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 1;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 1;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;  
    
  case P1_QUAD:                   /* linear on quadrilateral */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 3;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 4;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 1;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;  
    
    
  case BILINEAR_QUAD:             /* bilinear quadrilateral */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 4;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 4;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 2;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case C_BILINEAR_QUAD:              /* bilinear quadrilateral with additional centroid node */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 5;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 4;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 2;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;
    
  case S_BIQUAD_QUAD:          /* biquadratic serendipity quadrilateral */
    switch( info ){              /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 8;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 9;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 3;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;
    
  case BIQUAD_QUAD:               /* biquadratic quadrilateral */
    switch( info ){              /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 9;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 9;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 3;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case BIQUAD_QUAD_LS:            /* biquadratic quadrilateral for level set*/
    switch( info ){              /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 9;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 25;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 3;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

    
  case P0_HEX:                    /* bilinear quadrilateral */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 1;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 1;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 1;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;
    
  case P1_HEX:                    /* bilinear quadrilateral */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 4;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 8;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 1;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;
    
  case TRILINEAR_HEX:             /* trilinear hexahedron */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 8;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 8;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 4;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 2;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case C_TRILINEAR_HEX:           /* trilinear hexahedron with additional centroid node */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 9;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 8;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 4;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 2;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;
    
  case S_TRIQUAD_HEX:            /* serendipity triquadratic hexahedron */
    switch( info ){              /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 20;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 27;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 9;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 3;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;
    
  case TRIQUAD_HEX:               /* triquadratic hexahedron */
    switch( info ){               /* select type of information required*/
    case NNODES:                  /* number of nodes */
      answer = 27;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 27;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 9;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 3;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case LINEAR_TET:                /* linear tetrahedron */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 4;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 4;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 4;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 2;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;    
    /*
     * For 1D elements "surface" and "edge" are odd terms. Perhaps DIM-1 and DIM-2
     * is the way to think of them where, for BAR elements, DIM = 1.  Hence,
     * "edge" is a -1 dimensional thing and is meaningless...  This idea is consistent
     * with what existed in find_surf_st() and find_edge_s().
     */
    
  case LINEAR_BAR:                /* 1D, linear element for shell elements */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 2;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 2;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 1;
      break;
    case NQUAD_SURF:            /* number of "surface" quad points */
      answer = 1;
      break;
    case NQUAD_EDGE:            /* number of "edge" quad points */
      answer = 0;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;
    
  case QUAD_BAR:                /* 1D, quadratic element for shell elements */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 3;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 3;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 1;
      break;
    case NQUAD_SURF:            /* number of "surface" quad points */
      answer = 1;
      break;
    case NQUAD_EDGE:            /* number of "edge" quad points */
      answer = 0;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case BILINEAR_SHELL:                /* 1D, linear element for shell elements */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 4;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 4;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of "surface" quad points */
      answer = 2;
      break;
    case NQUAD_EDGE:            /* number of "edge" quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case BIQUAD_SHELL:                /* 1D, linear element for shell elements */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 9;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 9;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of "surface" quad points */
      answer = 2;
      break;
    case NQUAD_EDGE:            /* number of "edge" quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case BILINEAR_TRISHELL:         /* linear triangular shell */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 3;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 3;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:              /* number of surface quad points */
      answer = 2;
      break;
    case NQUAD_EDGE:              /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case P1_SHELL:                  /* linear discontinuous on shell */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 3;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 4;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 1;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;

  case P0_SHELL:                   /* constant discontinuous on shell */
    switch( info ){               /* select type of information required */
    case NNODES:                  /* number of nodes */
      answer = 1;
      break;
    case NQUAD:                   /* number of quadrature points */
      answer = 1;
      break;
    case NDIM:                    /* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:            /* number of surface quad points */
      answer = 1;
      break;
    case NQUAD_EDGE:            /* number of edge quad points */
      answer = 1;
      break;
    default:
      fprintf(stderr, "Unknown quantity\n");
      answer = -1;
      break;
    }
    break;


  default:
    fprintf(stderr, "Element itype = %d\n", ielem_type);
    fprintf(stderr, "Unknown or unimplemented element type.\n");
          answer = -1;
          break;
  }
  return answer;
} /* END of routine elem_info   */
/*****************************************************************************/

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int
dof_lnode_interp_type(const int n, const int Element_Type,
                      const int interp_type, const int edge)
    
    /**********************************************************************
     *
     * dof_lnode_interp_type():
     *
     * This routine returns the number of degrees of freedom for a variable
     * in an elemental interpolation given the local node number,
     * an element_type, and the interpolation type.
     *
     *
     * Element_Shape    LINE_SEGMENT    (1D)
     *
     *                  TRIANGLE        (2D)
     *                  QUADRILATERAL   (2D)
     *                  SHELL           (2D/3D)
     *
     *                  TETRAHEDRON     (3D)
     *                  PRISM           (3D)
     *                  HEXAHEDRON      (3D)
     *
     *
     * Note that in addition to the element shape, the node ordering convention 
     * is important so that, based on node number alone, we can identify corners, 
     * midside nodes, centroid nodes, etc. We take the node ordering convention
     * from EXODUS II which is merely the PATRAN node ordering convention...
     *
     * For example, we rely on the fact that local node [8] for quadrilaterals is
     * the centroid node. Thus, even if you only desired Q1 velocity interpolation
     * and P0 pressure at the centroid, you would need to generate 9 node quads
     * so that the centroid node is included.
     *
     * The argument, edge, is needed for the I_SP interpolation. In this
     * interpolation the basis function are bilinear in the interiors, but
     * biquadratic at the boundaries. Nodes on boundaries have edge equal to
     * one, but those in the interior have edge equal to zero.
     *
     * This is now changed so that pressures always appear at the 1st (0th) node
     * of any element (this has no effect on the results, I think)
     *
     * return value:
     *          0 -- if variable is not active at this node in this element
     *          1 -- if variable has 1 dof at this node in this element
     *          n -- if variable has n dof at this node in this element
     *         -1 -- if something is horribly wrong.
     *
     * NOTE: HKM: This function is currently unused. It's really a
     *       duplicate of node_info() without the issues involved with
     *       discontinuous variables thrown in. It will later replace
     *       node_info()'s functionality after discontinuous variables
     *       are reworked.
     ***********************************************************************/
{
  /*
   * All local variables made static to speed up calling overhead of this
   * function
   */
  static int Element_Shape;
  if (n < 0)  return (-1);

  /*
   * Now convert from SHM element type (cf. el_elm.h) into more fundamental
   * Element_Shape...
   */
  Element_Shape = type2shape(Element_Type);

  /*
   * If the variable has no interpolation then it means
   * that it is not active and hence no contribution to 
   * degrees of freedom at this node.
   */
  if (interp_type == I_NOTHING) return(0);

  /*
   * Let's categorize cases based on the dimension so that the number of
   * the local node can be used as a hint of whether this node has an
   * active variable or not.
   */
  switch (Element_Shape) {
 
    /* 
     * One dimensional line segments...
     */               
  case LINE_SEGMENT:
      switch (interp_type) {
      case I_Q1:                /* 2 node, 1 dof/node, Lagrangian linear */
          return( ( n < 2 ) ? 1 : 0 );
      case I_Q2:                /* 3 node, 1 dof/node, Lagrangian quadratic */
          return( ( n < 3 ) ? 1 : 0 );
      case I_Q3:                /* 4 node, 1 dof/node, Lagrangian cubic */
          return( ( n < 4 ) ? 1 : 0 );
      case I_Q4:                /* 5 node, 1 dof/node, Lagrangian quartic */
          return( ( n < 5 ) ? 1 : 0 );
      case I_P0:                /* 1 node, 1 dof/node, piecewise constant  */
          return( ( n == 0 ) ? 1 : 0 );
      case I_P1:                /* 1 node, 2 dof/node, piecewise linear */
          return( ( n == 0 ) ? 2 : 0 );
      case I_H3:                /* 2 node 2-dof Hermite cubic */
          return( ( n < 2 ) ? 2 : 0 );
      default:
          EH(-1, "Unrecognized line segment interpolation.");
          break;
      }
      break;

      /*
       * Two dimensional triangles...
       */
  case TRIANGLE:
  case TRISHELL:
      switch (interp_type) {
      case I_Q1:                /* 3 node, 1 dof/node, Lagrangian linear */
          return( ( n < 3 ) ? 1 : 0 );
      case I_Q2:                /* 6 node, 1 dof/node, Lagrangian quadratic */
        return( ( n < 6 ) ? 1 : 0 );
      case I_Q2_LSA:
          return( ( n < 6 ) ? 1 : 0 );
      case I_P0:                /* 1 node, 1 dof/node, piecewise constant */
          return( ( n == 0 ) ? 1 : 0 );
      case I_P1:                /* 1 node, 3 dof/node, piecewise linear */
          return( ( n == 0 ) ? 3 : 0 );
      default:
          EH(-1, "node_interp_info: Unrecognized triangle interpolation.");
          break;
      }
      break;

      /*
       * Two dimensional quadrilaterals...
       */
  case QUADRILATERAL:
  case SHELL:
      switch (interp_type) {
      case I_Q1:                /* 4 node, 1 dof/node, Lagrangian bilinear */
      case I_Q1_D:              /* 4 node, 1 dof/node, Lagrangian bilinear
                                   at interfaces that are discontinuous */
      case I_Q1_GP:
      case I_Q1_GN:
          return( ( n < 4 ) ? 1 : 0 );
      case I_Q1_G:              /* 4 node, 2 dof/node, Lagrangian bilinear */
      case I_Q1_XV:
      case I_Q1_XG:
          return( ( n < 4 ) ? 2 : 0 );
      case I_Q1_HV:             /* 4 node, 1 dof/node, Lagrangian bilinear + discontinuous enrichment */
      case I_Q1_HG:
          if ( n == 8 ) return 1;
          return( ( n < 4 ) ? 1 : 0 );
      case I_Q1_HVG:            /* 4 node, 1 dof/node, Lagrangian bilinear + 2 dof discontinuous enrichment */
          if ( n == 8 ) return 2;
          return( ( n < 4 ) ? 1 : 0 );
      case I_Q2:                /* 9 node, 1 dof/node, Lagrangian biquadratic */
      case I_Q2_D:              /* 9 node, 1 dof/node, Lagrangian biquadratic
                                   at interfaces that are discontinuous */
      case I_Q2_LSA:
      case I_Q2_D_LSA:
      case I_Q2_GP:
      case I_Q2_GN:
          return( ( n < 9 ) ? 1 : 0 );
      case I_Q2_G:              /* 9 node, 2 dof/node, Lagrangian biquadratic */
      case I_Q2_XV:
      case I_Q2_XG:
          return( ( n < 9 ) ? 2 : 0 );
      case I_Q2_HV:             /* 9 node, 1 dof/node, Lagrangian bilinear + discontinuous enrichment */
      case I_Q2_HG:
          if ( n == 8 ) return 2;
          return( ( n < 9 ) ? 1 : 0 );
      case I_Q2_HVG:            /* 9 node, 1 dof/node, Lagrangian bilinear + 2 dof discontinuous enrichment */
          if ( n == 8 ) return 3;
          return( ( n < 9 ) ? 1 : 0 );
      case I_SP:                /* 4 node, 1 dof/node, Lagrangian bilinear, but
                                 * biquadratic at all boundaries */
          if (edge) return( ( n < 8 ) ? 1 : 0 );
          else      return( ( n < 4 ) ? 1 : 0 );
      case I_S2:                /* 8 node, 1 dof/node, serendipity */
          return( ( n < 8 ) ? 1 : 0 );
      case I_P0:                /* 1 node, 1 dof/node, piecewise constant */
      case I_P0_GP:
      case I_P0_GN:
          switch (Element_Type) {
          case BILINEAR_QUAD:
              /* use first node for pressure in bilinear elements */
              return( ( n == 0 ) ? 1 : 0 );
          case C_BILINEAR_QUAD:
              /* use node 4 for pressure in bilinear elements */
              return( ( n == 4 ) ? 1 : 0 );
          case BIQUAD_QUAD:
          case BIQUAD_QUAD_LS:
              /* use centroid node for pressure in biquadratic elements */
              return( ( n == 8 ) ? 1 : 0 );
          case S_BIQUAD_QUAD:
              return( ( n == 7 ) ? 1 : 0 );
          default:
              EH( -1, "node_interp_info: unrecognized element type ");
              break;
          }
          break;
      case I_P0_G:              /* 1 node, 2 dof/node, extended piecewise constant */
      case I_P0_XV:
          switch (Element_Type) {
          case BILINEAR_QUAD:
              /* use first node for pressure in bilinear elements */
              return( ( n == 0 ) ? 2 : 0 );
          case C_BILINEAR_QUAD:
              /* use node 4 for pressure in bilinear elements */
              return( ( n == 4 ) ? 2 : 0 );
          case BIQUAD_QUAD:
          case BIQUAD_QUAD_LS:
              /* use centroid node for pressure in biquadratic elements */
              return( ( n == 8 ) ? 2 : 0 );
          case S_BIQUAD_QUAD:
              return( ( n == 7 ) ? 2 : 0 );
          default:
              EH( -1, "node_interp_info: unrecognized element type ");
              break;
          }
          break;
      case I_P1:                /* 1 node, 3 dof/node, piecewise linear */
      case I_P1_GP:
      case I_P1_GN:
          switch (Element_Type) {
          case BILINEAR_QUAD:
              /* use first node for pressure in bilinear elements */
              return( ( n == 0 ) ? 3 : 0 );
          case C_BILINEAR_QUAD:
              /* use centroid node for discontinuous dofs */
              return( ( n == 4 ) ? 3 : 0 );
          case BIQUAD_QUAD:
          case BIQUAD_QUAD_LS:
              /* use centroid node for pressure in biquadratic elements */
              return( ( n == 8 ) ? 3 : 0 );
          case S_BIQUAD_QUAD:
              return( ( n == 7 ) ? 3 : 0 );
          default:
              EH(-1, "node_intero_info: unrecognized element type ");
              break;
          }
          break;
      case I_P1_G:              /* 1 node, 6 dof/node, piecewise linear */
      case I_P1_XV:
          switch (Element_Type) {
          case BILINEAR_QUAD:
              /* use first node for pressure in bilinear elements */
              return( ( n == 0 ) ? 6 : 0 );
          case C_BILINEAR_QUAD:
              /* use centroid node for discontinuous dofs */
              return( ( n == 4 ) ? 6 : 0 );
          case BIQUAD_QUAD:
          case BIQUAD_QUAD_LS:
              /* use centroid node for pressure in biquadratic elements */
              return( ( n == 8 ) ? 6 : 0 );
          case S_BIQUAD_QUAD:
              return( ( n == 7 ) ? 6 : 0 );
          default:
              EH(-1, "node_intero_info: unrecognized element type ");
              break;
          }
          break;
      case I_H3:                /* 4 node, 4 dof/node, Hermite bicubic */
          return( ( n < 4 ) ? 4 : 0 );
      case I_Q3:                /* 16 node, 1 dof/node, Lagrangian bicubic */
          return( ( n < 16 ) ? 1 : 0 );
      case I_PQ1:
          switch (Element_Type) {
          case C_BILINEAR_QUAD:
              /* use centroid node for discontinuous dofs  */
              return( ( n == 4 ) ? 4 : 0 );
          case BIQUAD_QUAD:
          case BIQUAD_QUAD_LS:
              /* Use centroid node for all dofs */
              return( ( n == 8 ) ? 4 : 0 );
          default:
              EH(-1,
                 "PQ1 interpolation not implemented for this Element Type.");
              break;
          }
          break;
      case I_PQ2:
          switch(Element_Type) {
          case BIQUAD_QUAD:
          case BIQUAD_QUAD_LS:
              /* Use centroid node for all dofs */
              return( ( n == 8 ) ? 9 : 0 );
          default:
              EH(-1,
                 "PQ2 interpolation not implemented for this Element Type.");
              break;
          }
          break;
      default:
          EH(-1, "Unrecognized quadrilateral interpolation.");
          break;
      }
      break;

      /*
       * Three dimensional tetrahedrons...
       */
  case TETRAHEDRON:
      switch (interp_type) {
      case I_Q1:                /* 4 node, 1 dof/node, Lagrangian linear */
      case I_Q1_D:              /* 4 node, 1 dof/node, Lagrangian linear 
                                   at interfaces that are discontinuous */
          return( ( n < 4 ) ? 1 : 0 );
      case I_Q2:                /* 10 node, 1 dof/node, Lagrangian quadratic */
          return( ( n < 10 ) ? 1 : 0 );
      case I_P0:                /* 1 node, 1 dof/node, piecewise constant */
          return( ( n == 0 ) ? 1 : 0 );
      case I_P1:                /* 1 node, 4 dof/node, piecewise linear */
          return( ( n == 0 ) ? 4 : 0 );
      default:
          EH(-1, "Unrecognized tetrahedron interpolation.");
          break;
      }
      break;

      /*
       * Three dimensional prisms...
       */
  case PRISM:
      switch (interp_type) {
      case I_Q1:                /* 6 node, 1 dof/node, Lagrangian linear */
          return( ( n < 6 ) ? 1 : 0 );
      case I_Q2:                /* 15 node, 1 dof/node, Lagrangian quadratic */
          return( ( n < 15 ) ? 1 : 0 );
      case I_P0:                /* 1 node, 1 dof/node, piecewise constant */
          return( ( n == 0 ) ? 1 : 0 );
      case I_P1:                /* 1 node, 4 dof/node, piecewise linear */
          return( ( n == 15 ) ? 4 : 0 );
      default:
          EH(-1, "Unrecognized prism interpolation.");
          break;
      }
      break;

      /*
       * Three dimensional hexahedrons...
       */
  case HEXAHEDRON:
      switch (interp_type) {
      case I_Q1:                /* 8 node, 1 dof/node, Lagrangian linear */
      case I_Q1_D:              /* 8 node, 1 dof/node, Lagrangian linear
                                   at interfaces that are discontinuous */
      case I_Q1_GP:
      case I_Q1_GN:
          return( ( n < 8 ) ? 1 : 0 );
      case I_Q1_G:              /* 8 node, 2 dof/node, Lagrangian bilinear */
      case I_Q1_XV:
      case I_Q1_XG:
        return( ( n < 8 ) ? 2 : 0 );
      case I_S2:                /* 20 node, 1 dof/node, serendipity */
          return( ( n < 20 ) ? 1 : 0 );
      case I_Q2:                /* 27 node, 1 dof/node, Lagrangian quadratic */
      case I_Q2_D:              /* 27 node, 1 dof/node, Lagrangian quadratic
                                   at interfaces that are discontinuous */
      case I_Q2_GP:
      case I_Q2_GN:
          return( ( n < 27 ) ? 1 : 0 );
      case I_Q2_G:              /* 27 node, 2 dof/node, Lagrangian biquadratic */
      case I_Q2_XV:
      case I_Q2_XG:
          return( ( n < 27 ) ? 2 : 0 );
      case I_P0:                /* 1 node, 1 dof/node, piecewise constant */
      case I_P0_GP:
      case I_P0_GN:
          switch (Element_Type) {
          case TRILINEAR_HEX:
              return( ( n == 0 ) ? 1 : 0 );
          case C_TRILINEAR_HEX:
              return( ( n == 8 ) ? 1 : 0 );
          case TRIQUAD_HEX:
              return( ( n == 20 ) ? 1 : 0 );/* centroid node */
          default:
              EH(-1, "Unrecognized hexahedron interpolation.");
              break;
          }
          break;
      case I_P0_G:              
      case I_P0_XV:             /* 1 node, 2 dof/node, piecewise constant */
          switch (Element_Type) {
          case TRILINEAR_HEX:
              return( ( n == 0 ) ? 2 : 0 );
          case C_TRILINEAR_HEX:
              return( ( n == 8 ) ? 2 : 0 );
          case TRIQUAD_HEX:
              return( ( n == 20 ) ? 2 : 0 );/* centroid node */
          default:
              EH(-1, "Unrecognized hexahedron interpolation.");
              break;
          }
          break;
      case I_P1:                /* 1 node, 4 dof/node, piecewise linear */
      case I_P1_GP:
      case I_P1_GN:
          switch (Element_Type) {
          case TRILINEAR_HEX:
              return( ( n == 0 ) ? 4 : 0 );
          case C_TRILINEAR_HEX:
              return( ( n == 8 ) ? 4 : 0 );
          case TRIQUAD_HEX:
              return( ( n == 20 ) ? 4 : 0 ); /* centroid node */
          default:
              EH(-1, "Unrecognized hexahedron interpolation.");
              break;
          }
          break;
      case I_PQ1:
          switch (Element_Type) {
          case TRILINEAR_HEX:
              return( ( n == 0 ) ? 8 : 0 );
          case C_TRILINEAR_HEX:
              /* use centroid node for discontinuous dofs  */
              return( ( n == 4 ) ? 8 : 0 );
          case TRIQUAD_HEX:
              /* Use centroid node for all dofs */
              return( ( n == 20 ) ? 8 : 0 );
          default:
              EH(-1,"PQ1 interpolation not implemented for this Element Type.");
              break;
          }
          break;
      case I_PQ2:
          switch (Element_Type) {
          case TRIQUAD_HEX:
              /* Use centroid node for all dofs */
              return( ( n == 20 ) ? 27 : 0 );
          default:
              EH(-1,"PQ2 interpolation not implemented for this Element Type.");
              break;
          }
          break;
      case I_H3:                /* 8 node, 8 dof/node, Hermite tricubic */
          return( ( n < 8 ) ? 8 : 0 );
      default:
          EH(-1, "Unrecognized hexahedron interpolation.");
          break;
      }
      break;
  default:
      EH(-1, "Bad element shape.");
      break;
  }
  EH(-1, "node_interp_info: We should not be here.");
  return (-1);

}  /* END of routine node_info dof_lnode_interp_type  */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/*
 * type2shape() -- convert from general element types into basic shapes
 *
 * The element types tell both the interpolation and the basic shape. Sometimes
 * it's nice to be able to separate the two concepts. Hence, this code to
 * convert the legacy variables into basic shapes.
 */

int 
type2shape(const int element_type)
{
  int shape=-1;

  switch (element_type) {
  case LINEAR_TRI:
  case QUAD_TRI:
  case QUAD6_TRI:
      shape = TRIANGLE;
      break;
  case BILINEAR_QUAD:
  case C_BILINEAR_QUAD:              
  case S_BIQUAD_QUAD:
  case BIQUAD_QUAD_LS:
  case BIQUAD_QUAD:
  case P0_QUAD:
  case P1_QUAD:
    shape = QUADRILATERAL;
    break;
  case TRILINEAR_HEX:
  case C_TRILINEAR_HEX:
  case S_TRIQUAD_HEX:
  case TRIQUAD_HEX:
  case P0_HEX:
  case P1_HEX:
    shape = HEXAHEDRON;
    break;
  case LINEAR_BAR:
  case QUAD_BAR:
    shape = LINE_SEGMENT;
    break;
  case BILINEAR_SHELL:
  case BIQUAD_SHELL:
  case P1_SHELL:
  case P0_SHELL:
    shape = SHELL;
    break;
  case BILINEAR_TRISHELL:
    shape = TRISHELL;
    break;
  case LINEAR_TET:
    shape = TETRAHEDRON;
    break;
  default:
    fprintf(stderr,"type2shape ERROR: unknown element type: %d\b",
            element_type);
    EH(-1, "What basic shape are these new element types?");
    break;
  }
  return(shape);
}

/*
 * The basic shapes have various numbers of sides in 2D and 3D...
 */

int
shape2sides(const int element_shape)
{
  int num_sides=-1;

  switch (element_shape)
    {
    case QUADRILATERAL:
    case SHELL:
      num_sides = 4;
      break;

    case HEXAHEDRON:
      num_sides = 6;
      break;

    case TRIANGLE:
    case TRISHELL:
      num_sides = 3;
      break;

    case TETRAHEDRON:
      num_sides = 4;
      break;

    case LINE_SEGMENT:
      num_sides = 2;
      break;

    case PRISM:
      num_sides = 5;
      break;
      
    case PYRAMID:
      num_sides = 5;
      break;

    default:
      EH(-1, "Unrecognized element shape.");
      break;
    }

  return(num_sides);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int
getdofs(const int element_shape, const int interpolation)

    /**************************************************************************
     *
     * getdofs():
     *
     * For a basic element shape and given order of interpolation find out
     * the maximum number of degrees of freedom we have in the interpolation
     * for that basis function within the element. Do not include degrees of
     * freedom that are active in the element that may not be part of the
     * interpolation of the basis functions within the element, e.g., I_Q1_D
     * may have additional  degrees of freedom on the side of the element.
     * However, we won't count them below, so I_Q1_D returns 4 for a
     * quadralateral.
     *
     **************************************************************************/
{
  switch ( element_shape )
    {
      /* 
       * One dimensional line segments...
       */
    case LINE_SEGMENT:
      switch ( interpolation )
        {
        case I_Q1:              /* 2 node, 1 dof/node, Lagrangian linear */
        case I_Q1_D:
        case I_Q1_GP:
        case I_Q1_GN:
          return(2);
        case I_Q1_G:            /* 2 node, 2 dof/node, Lagrangian bilinear */
        case I_Q1_XV:
        case I_Q1_XG:
          return(4);

        case I_Q2:              /* 3 node, 1 dof/node, Lagrangian quadratic */
        case I_Q2_D:
        case I_Q2_GP:
        case I_Q2_GN:
          return(3);
        case I_Q2_G:            /* 3 node, 2 dof/node, Lagrangian biquadratic */
        case I_Q2_XV:
        case I_Q2_XG:
          return(6);

        case I_Q3:              /* 4 node, 1 dof/node, Lagrangian cubic */
          return(4);

        case I_Q4:              /* 5 node, 1 dof/node, Lagrangian quartic */
          return(5);

        case I_P0:              /* 1 node, 1 dof/node, piecewise constant  */
        case I_P0_GP:
        case I_P0_GN:
          return(1);

        case I_P1:              /* 1 node, 2 dof/node, piecewise linear */
        case I_P1_GP:
        case I_P1_GN:
          return(2);

        case I_H3:              /* 2 node 2-dof Hermite cubic */
          return(4);

        default:
          EH(-1, "Unrecognized line segment interpolation.");
          break;
        }
        break;
    case SHELL:
      switch ( interpolation )
        {
        case I_Q1:              /* 4 node, 1 dof/node, Lagrangian linear */
        case I_Q1_D:
        case I_Q1_GP:
        case I_Q1_GN:
          return(4);
        case I_Q1_G:            /* 4 node, 2 dof/node, Lagrangian bilinear */
        case I_Q1_XV:
        case I_Q1_XG:
          return(8);

        case I_Q2:              /* 9 node, 1 dof/node, Lagrangian quadratic */
        case I_Q2_D:
        case I_Q2_GP:
        case I_Q2_GN:
          return(9);
        case I_Q2_G:            /*  9 node, 2 dof/node, Lagrangian biquadratic */
        case I_Q2_XV:
        case I_Q2_XG:
          return(18);

        case I_Q3:              /* 4 node, 1 dof/node, Lagrangian cubic */
          return(4);

        case I_Q4:              /* 5 node, 1 dof/node, Lagrangian quartic */
          return(5);

        case I_P0:              /* 1 node, 1 dof/node, piecewise constant  */
        case I_P0_GP:
        case I_P0_GN:
          return(1);

        case I_P1:              /* 1 node, 2 dof/node, piecewise linear */
        case I_P1_GP:
        case I_P1_GN:
          return(4);

        case I_H3:              /* 4 node 2-dof Hermite cubic */
          return(8);

        default:
          EH(-1, "Unrecognized SHELL interpolation.");
          break;
        }
        break;

      /*
       * Two dimensional triangles...
       */
    case TRIANGLE:
    case TRISHELL:
      switch ( interpolation )
        {
        case I_Q1:              /* 3 node, 1 dof/node, Lagrangian linear */
        case I_Q1_D:
          return(3);

        case I_Q2:              /* 6 node, 1 dof/node, Lagrangian quadratic */
        case I_Q2_D:
        case I_Q2_LSA:
        case I_Q2_D_LSA:
          return(6);

        case I_P0:              /* 1 node, 1 dof/node, piecewise constant */
          return(1);

        case I_P1:              /* 1 node, 3 dof/node, piecewise linear */
          return(3);

        default:
          EH(-1, "Unrecognized triangle interpolation.");
          break;
        }
        break;

      /*
       * Two dimensional quadrilaterals...
       */
    case QUADRILATERAL:
      switch ( interpolation )
        {
        case I_Q1:              /* 4 node, 1 dof/node, Lagrangian bilinear */
        case I_Q1_D:
        case I_Q1_GP:
        case I_Q1_GN:
          return(4);
        case I_Q1_G:            /* 4 node, 2 dof/node, Lagrangian bilinear */
        case I_Q1_XV:
        case I_Q1_XG:
          return(8);
        case I_Q1_HV:           /* 4 node, 1 dof/node, Lagrangian bilinear + discontinuous enrichment */
        case I_Q1_HG:
          return(5);
        case I_Q1_HVG:          /* 4 node, 1 dof/node, Lagrangian bilinear + 2 dof discontinuous enrichment */
          return(6);
        case I_Q2:              /* 9 node, 1 dof/node, Lagrangian biquadratic */
        case I_Q2_D:
        case I_Q2_LSA:
        case I_Q2_D_LSA:
        case I_Q2_GP:
        case I_Q2_GN:
          return(9);
        case I_Q2_G:            /* 9 node, 2 dof/node, Lagrangian biquadratic */
        case I_Q2_XV:
        case I_Q2_XG:
          return(18);
        case I_Q2_HV:           /* 9 node, 1 dof/node, Lagrangian bilinear + discontinuous enrichment */
        case I_Q2_HG:
          return(10);
        case I_Q2_HVG:          /* 9 node, 1 dof/node, Lagrangian bilinear + 2 dof discontinuous enrichment */
          return(11);

        case I_SP:              /* 8 node, 0 or 1 dof/node, Lagrangian biquadratic subparametric*/
          return(8);

        case I_S2:              /* 8 node, 1 dof/node, serendipity */
          return(8);

        case I_P1:              /* 1 node, 3 dof/node, piecewise linear */
        case I_P1_GP:
        case I_P1_GN:
          return(3);

        case I_P0:              /* 1 node, 1 dof/node, piecewise constant */
        case I_P0_GP:
        case I_P0_GN:
          return(1);

        case I_P0_G:            /* 1 node, 2 dof/node, extended piecewise constant */
        case I_P0_XV:
          return(2);

        case I_P1_G:            /* 1 node, 6 dof/node, extended piecewise linear */
        case I_P1_XV:
          return(6);

        case I_H3:              /* 4 node, 4 dof/node, Hermite bicubic */
          return(16);

        case I_Q3:              /* 16 node, 1 dof/node, Lagrangian bicubic */
          return(16);

        case I_PQ1:            /* centroid node, 4 dof/node, bilinear discontinous */
          return(4);

        case I_PQ2:
          return(9);           /* centroid node, 9 dof/node, biquadratic discontinous */

        default:
          EH(-1, "Unrecognized quadrilateral interpolation.");
          break;
        }
        break;

      /*
       * Three dimensional tetrahedrons...
       */
    case TETRAHEDRON:
      switch ( interpolation )
        {
        case I_Q1:              /* 4 node, 1 dof/node, Lagrangian linear */
        case I_Q1_D:
          return(4);

        case I_Q2:              /* 10 node, 1 dof/node, Lagrangian quadratic */
        case I_Q2_D:
          return(10);

        case I_P0:              /* 1 node, 1 dof/node, piecewise constant */
          return(1);

        case I_P1:              /* 1 node, 4 dof/node, piecewise linear */
          return(4);

        default:
          EH(-1, "Unrecognized tetrahedron interpolation.");
          break;
        }
        break;

      /*
       * Three dimensional prisms...
       */
    case PRISM:
      switch ( interpolation )
        {
        case I_Q1:              /* 6 node, 1 dof/node, Lagrangian linear */
        case I_Q1_D:
          return(6);

        case I_Q2:              /* 15 node, 1 dof/node, Lagrangian quadratic */
        case I_Q2_D:
          return(15);

        case I_P0:              /* 1 node, 1 dof/node, piecewise constant */
          return(1);

        case I_P1:              /* 1 node, 4 dof/node, piecewise linear */
          return(4);

        default:
          EH(-1, "Unrecognized prism interpolation.");
          break;
        }
        break;


      /*
       * Three dimensional hexahedrons...
       */
    case HEXAHEDRON:
      switch ( interpolation )
        {
        case I_Q1:              /* 8 node, 1 dof/node, Lagrangian linear */
        case I_Q1_D:
        case I_Q1_GP:
        case I_Q1_GN:
          return(8);

        case I_Q1_G:            /* 8 node, 2 dof/node, Lagrangian bilinear */
        case I_Q1_XV:
        case I_Q1_XG:
          return(16);

        case I_S2:              /* 20 node, 1 dof/node, serendipity */
          return(20);

        case I_Q2:              /* 27 node, 1 dof/node, Lagrangian quadratic */
        case I_Q2_D:
        case I_Q2_GP:
        case I_Q2_GN:
          return(27);

        case I_Q2_G:            /* 27 node, 2 dof/node, Lagrangian biquadratic */
        case I_Q2_XV:
        case I_Q2_XG:
          return(54);

        case I_P0:              /* 1 node, 1 dof/node, piecewise constant */
        case I_P0_GP:
        case I_P0_GN:
          return(1);

        case I_P0_XV:           /* 1 node, 1 dof/node, piecewise constant */
        case I_P0_G:
          return(2);

        case I_P1:              /* 1 node, 4 dof/node, piecewise linear */
        case I_P1_GP:
        case I_P1_GN:
          return(4);
        case I_P1_XV:
        case I_P1_G:
          return(8);

        case I_H3:              /* 8 node, 8 dof/node, Hermite tricubic */
          return(64);

        default:
          EH(-1, "Unrecognized hexahedron interpolation.");
          break;
        }
        break;

    default:
      EH(-1, "Bad element shape.");
      return -1;
    }
  EH(-1, "We should not be here.");
  return (-1);
  
} /* END of routine getdofs  */
/*****************************************************************************/

void
find_stu(const int   iquad,     /* current GQ index  */
         const int   ielem_type, /* element type      */
         dbl   *s,              /* local          */
         dbl   *t,              /* GQ coordinates */
         dbl   *u  )            /* (returned)     */

/*
*       Function which determines the Gaussian-quadrature points s, t and u
*       from the given index.
*
*        Author:          Scott Hutchinson (1421)
*        Date:            19 May 1992
*        Revised:         24 May 1995 RAC
*
*        NOTE: this routine currently assumes element is a quad or a hex
*
*        Revision Histor:
*        PKN, 6 May 2002, extend for line segment (1D) elements.
*        PRS, 24 August 2009, extend for 2D shell elements.
*        SAR, 4 November 2011, fixed error for tets.
*/

{
/* LOCAL VARIABLES */
                           /*  1.0 / sqrt (3.0)  */
  static const double Ftemp1 =  0.57735026918962584208;
                           /*  sqrt (3.0/5.0)   */
  static const double Ftemp2 =  0.77459666924148340428;

  static const double Tri1 = 0.16666666666666666667;
  static const double Tri2 = 0.66666666666666666667;
  static const double Tri3 = 0.091576213509771;
  static const double Tri4 = 0.816847572980459;
  static const double Tri5 = 0.108103018168070;
  static const double Tri6 = 0.445948490915965;
  // static const double Tri7 = 0.33333333333333333333;
  static const double Tri8 = 0.873821971016996;
  static const double Tri9 = 0.063089014491502;
  static const double Tri10= 0.501426509658179;
  static const double Tri11= 0.249286745170910;
  static const double Tri12= 0.636502499121399;
  static const double Tri13= 0.310352451033785;
  static const double Tri14= 0.053145049844816;
  static const double quad5_1 =  0.9061798459;
  static const double quad5_2 =  0.5384693101;

  switch( ielem_type ){                 /* select element */

  case LINEAR_TRI:                   /* linear triangle */
  case BILINEAR_TRISHELL:            /* linear triangular shell */
    *s = (iquad == 0) ? Tri2 : Tri1;
    *t = (iquad == 1) ? Tri2 : Tri1;
    *u = 0.0;
    break;

  case QUAD_TRI:                     /* quadratic triangle */
    /* 6 point */
    if ( iquad < 3) {
      *s = (iquad == 0) ? Tri4 : Tri3;
      *t = (iquad == 1) ? Tri4 : Tri3;
      *u = 0.0;
    } else {
      *s = (iquad == 4) ? Tri5 : Tri6;
      *t = (iquad == 5) ? Tri5 : Tri6;
      *u = 0.0;
    }
    /* 4 point
    if ( iquad == 0) {
      *s = Tri7;
      *t = Tri7;
      *u = 0.0;
    } else {
      *s = (iquad == 1) ? 0.6 : 0.2;
      *t = (iquad == 2) ? 0.6 : 0.2;
      *u = 0.0;
    }
    */
    break;

  case QUAD6_TRI:                     /* quadratic triangle with 6th order quadrature */
    /* 12 point */
    if ( iquad < 3) {
      *s = (iquad == 0) ? Tri8 : Tri9;
      *t = (iquad == 1) ? Tri8 : Tri9;
      *u = 0.0;
    } else if ( iquad < 6) {
      *s = (iquad == 4) ? Tri10 : Tri11;
      *t = (iquad == 5) ? Tri10 : Tri11;
      *u = 0.0;
    } else if ( iquad < 8) {
      *s = (iquad == 6) ? Tri12 : Tri13;
      *t = (iquad == 7) ? Tri12 : Tri13;
      *u = 0.0;
    } else if ( iquad < 10) {
      *s = (iquad == 8) ? Tri14 : Tri12;
      *t = (iquad == 9) ? Tri14 : Tri12;
      *u = 0.0;
    } else {
      *s = (iquad == 10) ? Tri14 : Tri13;
      *t = (iquad == 11) ? Tri14 : Tri13;
      *u = 0.0;
    }
    /* 4 point
    if ( iquad == 0) {
      *s = Tri7;
      *t = Tri7;
      *u = 0.0;
    } else {
      *s = (iquad == 1) ? 0.6 : 0.2;
      *t = (iquad == 2) ? 0.6 : 0.2;
      *u = 0.0;
    }
    */
    break;

  case P0_QUAD:                         /* constant discontinuous on quadrilateral */
  case P0_SHELL:                        /* constant discontinuous on shell */
  case P0_HEX:                          /* constant discontinuous on hexahedron */
    *s = 0.0;
    *t = 0.0;
    *u = 0.0;
    break;

  case BILINEAR_QUAD:                   /* bilinear quadrilateral */
  case C_BILINEAR_QUAD:              /* bilinear quadrilateral with additional centroid node */
  case P1_QUAD:                         /* linear discontinuous on quadrilateral */
  case P1_SHELL:                        /* linear discontinuous on shell */
    *s = (iquad%2 == 0) ? Ftemp1 : -Ftemp1;
    *t = (iquad < 2)    ? Ftemp1 : -Ftemp1;
    *u = 0.0;
    break;

  case S_BIQUAD_QUAD:                /* biquadratic serendipity quadrilateral */
  case   BIQUAD_QUAD:                     /* biquadratic quadrilateral */
    if (iquad%3 == 0)
      *s = Ftemp2;
    else if ( (iquad-1)%3 == 0)
      *s = 0.0;
    else
      *s = -Ftemp2;

    if (iquad < 3)
      *t = Ftemp2;
    else if (iquad < 6)
      *t = 0.0;
    else
      *t = -Ftemp2;

    *u = 0.0;
    break;

 /* biquadratic quadrilateral for level set */
  case   BIQUAD_QUAD_LS:
    if (iquad%5 == 0)
      *s = quad5_1;
    else if ( (iquad-1)%5 == 0)
      *s = quad5_2;
    else if ( (iquad-2)%5 == 0)
      *s = 0.0;
    else if ( (iquad-3)%5 == 0)
      *s = -quad5_2;
    else
      *s = -quad5_1;

    if (iquad < 5)
      *t = quad5_1;
    else if (iquad < 10)
      *t = quad5_2;
    else if (iquad < 15)
      *t = 0.0;
    else if (iquad < 20)
      *t = -quad5_2;
    else
      *t = -quad5_1;

    *u = 0.0;
    break; 


  case TRILINEAR_HEX:                   /* trilinear hexahedron */
  case C_TRILINEAR_HEX:                 /* trilinear hexahedron with additional centroid node */
  case P1_HEX:                          /* linear discontinuous on hexahedron */
    *s = (iquad%2 == 0)            ? Ftemp1 : -Ftemp1;
    *t = (iquad - 4*(iquad/4) < 2) ? Ftemp1 : -Ftemp1;
    *u = (iquad < 4)               ? Ftemp1 : -Ftemp1;
    break;

  case S_TRIQUAD_HEX:                 /* serendipity triquadratic hexahedron */
  case   TRIQUAD_HEX:                     /* triquadratric hexahedron */
    if (iquad%3 == 0)
      *s = Ftemp2;
    else if ( (iquad-1)%3 == 0)
      *s = 0.0;
    else
      *s = -Ftemp2;

    if (iquad - 9*(iquad/9) < 3)
      *t = Ftemp2;
    else if (iquad - 9*(iquad/9) < 6)
      *t = 0.0;
    else
      *t = -Ftemp2;

    if (iquad < 9)
      *u = Ftemp2;
    else if (iquad < 18)
      *u = 0.0;
    else
      *u = -Ftemp2;
    break;

  case LINEAR_BAR:
    if ( iquad == 0)
      *s = Ftemp1;
    else
      *s = -Ftemp1;
    *t = *u = 0.0;
    break;
      
  case QUAD_BAR:
    if ( iquad == 0)
      *s = Ftemp2;
    else if ( iquad == 1)
      *s = 0.0;
    else
      *s = -Ftemp2;
    *t = *u = 0.0;
    break;

  case BILINEAR_SHELL:
    *s = (iquad%2 == 0) ? Ftemp1 : -Ftemp1;
    *t = (iquad < 2)    ? Ftemp1 : -Ftemp1;
    *u = 0.0;
    break;
      
  case BIQUAD_SHELL:
    if (iquad%3 == 0)
      *s = Ftemp2;
    else if ( (iquad-1)%3 == 0)
      *s = 0.0;
    else
      *s = -Ftemp2;

    if (iquad < 3)
      *t = Ftemp2;
    else if (iquad < 6)
      *t = 0.0;
    else
      *t = -Ftemp2;

    *u = 0.0;
    break;

  case LINEAR_TET:
    //*s = *t = *u = 0.25;
  { static const double alpha = 0.585410196624969;
    static const double beta = .138196601125011;
    switch (iquad) {
    case 0:
      *s = alpha;
      *t = *u = beta;
      break;
    case 1:
      *s = *u = beta;
      *t = alpha;
      break;
    case 2:
      *s = *t = beta;
      *u = alpha;
      break;
    case 3:
      *s = *t = *u = beta;
      break;
    }
  }

 /*   switch (iquad )
      {
      case 0: *s = *t = *u = 0.25; break;
      case 1: *s = one_sixth; *t = one_sixth; *u = one_sixth; break;
      case 2: *s = one_sixth; *t = one_sixth; *u = one_half ; break;
      case 3: *s = one_sixth; *t = one_half ; *u = one_sixth; break;
      case 4: *s = one_half ; *t = one_sixth; *u = one_sixth; break;
      } */
    break;
      
  default:
    EH(-1, "Unknown or unimplemented element type.\n");
    break;
  }

} /* END of routine find_stu  */
/*****************************************************************************/

void
find_surf_st(const int iquad,           /* current GQ index */
             const int ielem_type,      /* element type */
             const int iside,           /* current side of element */
             const int dim,             /* dimensions of element */
             double xi[DIM],      /* (returned) local GQ coordinates for surface integral */
             double *s,           /* Gaussian-quadrature points (s, t) */
             double *t,
             double *u)

/*
*       Function which determines the Gaussian-quadrature points (s, t)
*       for surface integration from the given index, iquad, and element
*       type, ielem_type.
*
*        Author:          Scott Hutchinson (1421)
*        Date:            19 May 1992
*        Revised:         24 May 1995 RAC
*        Revised:         04 Nov 2011 SAR tets
*
*        NOTE: this routine currently assumes element is a quad or a hex
*
*  TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98.
*/

{
  int i_s = -1;                 /* index of first local coord */
  int i_t = -1;                 /* index of second local coord */
                           /*  1.0 / sqrt (3.0)  */

  static double one_third = 3.333333333333333e-01;
  // static double one_sixth = 1.666666666666667e-01;

  static const double Ftemp1 =  0.57735026918962584208;

                           /*  sqrt (3.0/5.0)   */
  static const double Ftemp2 =  0.77459666924148340428;

  static const double quad5_1 =  0.90617984593866399280;
  static const double quad5_2 =  0.53846931010568309104;

  static const double Tri1 =  0.5*(1.-0.57735026918962584208);
  static const double Tri2 =  0.5*(1.+0.57735026918962584208);
  static const double Tri3 =  0.5*(1.-0.77459666924148340428);
  static const double Tri4 =  0.5*(1.+0.77459666924148340428);

  static const double Tri5 =  0.5*(1.-0.90617984593866399280);
  static const double Tri6 =  0.5*(1.+0.90617984593866399280);
  static const double Tri7 =  0.5*(1.-0.53846931010568309104);
  static const double Tri8 =  0.5*(1.+0.53846931010568309104);

  // static double Tri_temp1 = 0.2113248654051871;  /* 1/2 * (1 - 1/sqrt(3)) */
  // static double Tri_temp2 = 0.1127016653792583;  /* 1/2 * (1 - sqrt(3/5)) */

  /* the code further below is appropriate for line_segs, quads, and 
   * hexes but not for triangles
   */
  if ( ielem_type == LINEAR_TRI || ielem_type == BILINEAR_TRISHELL ) {

    xi[2] = 0.;
    *t = 0.;
    switch (iside) {
    case 1:
      switch( iquad ){
        case 0: xi[0] = 1. - Tri1; xi[1] = *s =  Tri1; break;
        case 1: xi[0] = 1. - Tri2; xi[1] = *s =  Tri2; break;

      }
      break;
    case 2:
      switch( iquad ){
        case 0: xi[0] = 0.; xi[1] = *s =  Tri1; break;
        case 1: xi[0] = 0.; xi[1] = *s =  Tri2; break;
      }
      break;
    case 3:
      switch( iquad ){
        case 0: xi[0] = *s = Tri1; xi[1] = 0.; break;
        case 1: xi[0] = *s = Tri2; xi[1] = 0.; break;
      }
      break;
    default:
      EH(-1,"Illegal side number for LINEAR_TRI element");
      break;
    }
    return;
  }

  if ( ielem_type == QUAD_TRI ) {
    xi[2] = 0.;
    *t = 0.;
    switch (iside) {
    case 1:
      switch( iquad ){
        case 0: xi[0] = 1. - Tri3; xi[1] = *s =  Tri3; break;
        case 1: xi[0] = 0.5; xi[1] = *s =  0.5; break;
        case 2: xi[0] = 1. - Tri4; xi[1] = *s =  Tri4; break;
      }
      break;
    case 2:
      switch( iquad ){
        case 0: xi[0] = 0.; xi[1] = *s =  Tri4; break;
        case 1: xi[0] = 0.; xi[1] = *s =  0.5; break;
        case 2: xi[0] = 0.; xi[1] = *s =  Tri3; break;
      }
      break;
    case 3:
      switch( iquad ){
        case 0: xi[0] = *s = Tri3; xi[1] = 0.; break;
        case 1: xi[0] = *s = 0.5; xi[1] = 0.; break;
        case 2: xi[0] = *s = Tri4; xi[1] = 0.; break;
      }
      break;
    default:
      EH(-1,"Illegal side number for QUAD_TRI element");
      break;
    }
    return;
  }
  if ( ielem_type == QUAD6_TRI ) {
    xi[2] = 0.;
    *t = 0.;
    switch (iside) {
    case 1:
      switch( iquad ){
        case 0: xi[0] = 1. - Tri5; xi[1] = *s =  Tri5; break;
        case 1: xi[0] = 1. - Tri7; xi[1] = *s =  Tri7; break;
        case 2: xi[0] = 0.5; xi[1] = *s =  0.5; break;
        case 3: xi[0] = 1. - Tri8; xi[1] = *s =  Tri8; break;
        case 4: xi[0] = 1. - Tri6; xi[1] = *s =  Tri6; break;
      }
      break;
    case 2:
      switch( iquad ){
        case 0: xi[0] = 0.; xi[1] = *s =  Tri5; break;
        case 1: xi[0] = 0.; xi[1] = *s =  Tri7; break;
        case 2: xi[0] = 0.; xi[1] = *s =  0.5; break;
        case 3: xi[0] = 0.; xi[1] = *s =  Tri8; break;
        case 4: xi[0] = 0.; xi[1] = *s =  Tri6; break;
      }
      break;
    case 3:
      switch( iquad ){
        case 0: xi[0] = *s = Tri5; xi[1] = 0.; break;
        case 1: xi[0] = *s = Tri7; xi[1] = 0.; break;
        case 2: xi[0] = *s = 0.5; xi[1] = 0.; break;
        case 3: xi[0] = *s = Tri8; xi[1] = 0.; break;
        case 4: xi[0] = *s = Tri6; xi[1] = 0.; break;
      }
      break;
    default:
      EH(-1,"Illegal side number for QUAD6_TRI element");
      break;
    }
    return;
  }  
  *t = 0.; /* to avoid uninitialized memory read later */

/* first calculate local coordinates on the n-dimensional surface of the element 
 * then translate them to the n+1 dimensional elemental coords */
  switch (dim){
  case 1:
    xi[1] = 0.;
    xi[2] = 0.;
    switch(iside) {
    case 1:
      xi[0] = -1.;
      break;
    case 2:
      xi[0] = 1.;
      break;
    default:
      EH(-1,"Illegal side number for 1-D element");
      break;
    }
    break;
  case 2:
    xi[2] = 0.;
    switch (iside) {
    case 1:
      xi[1] = -1.;
      i_s = 0;
      break;
    case 2:
      xi[0] = 1.;
      i_s = 1;
      break;
    case 3:
      xi[1] =  1.;
      i_s = 0;
      break;
    case 4:
      xi[0] = -1.;
      i_s = 1;
      break;
    default:
      EH(-1,"Illegal side number for 2-D element");
      break;
    }
    break;
  case 3:
    switch (iside) {
    case 1:
      xi[1] = -1.;
      i_s = 0;
      i_t = 2;
      break;
    case 2:
      xi[0] =  1.;
      i_s = 1;
      i_t = 2;
      break;
    case 3:
      xi[1] =  1.;
      i_s = 0;
      i_t = 2;
      break;
    case 4:
      xi[0] = -1.;
      i_s = 1;
      i_t = 2;
      break;
    case 5:
      xi[2] = -1.;
      i_s = 0;
      i_t = 1;
      break;
    case 6:
      xi[2] =  1.;
      i_s = 0;
      i_t = 1;
      break;
    default:
      EH(-1,"Illegal side number for 3-D element");
      break;
    }
    break;
  default:
    EH(-1,"Illegal element dimension");
    break;
  }

  switch( ielem_type ){                 /* select element */

  case P0_QUAD:                         /*constant discontinuous quadrilateral */
  case P0_SHELL:                        /*constant discontinuous shell */
  case P1_QUAD:                         /*linear discontinuous quadrilateral */
  case P1_SHELL:                        /*linear discontinuous shell */
    switch( iquad ){
        case 0: xi[i_s] = *s =  0.; break;
    }
    break;

  case BILINEAR_QUAD:                   /* bilinear quadrilateral */
  case BILINEAR_SHELL:               /* bilinear shell */
  case C_BILINEAR_QUAD:              /* bilinear quadrilateral with additional centroid node */
    switch( iquad ){
        case 0: xi[i_s] = *s =  Ftemp1; break;
        case 1: xi[i_s] = *s = -Ftemp1; break;
    }
    break;

  case S_BIQUAD_QUAD:                /* biquadratic serendipity quadrilateral */
  case   BIQUAD_QUAD:                /* biquadratic quadrilateral */
  case   BIQUAD_SHELL:               /* biquadratic shell */
   switch( iquad ){
        case 0: xi[i_s] = *s = Ftemp2; break;
        case 1: xi[i_s] = *s = 0.0   ; break;
        case 2: xi[i_s] = *s =-Ftemp2; break;
    }
    break;

  case BIQUAD_QUAD_LS:               /* biquadratic quadrilateral for level set */
   switch( iquad ){
        case 0: xi[i_s] = *s = quad5_1; break;
        case 1: xi[i_s] = *s = quad5_2; break;
        case 2: xi[i_s] = *s = 0.0   ; break;
        case 3: xi[i_s] = *s =-quad5_2; break;
        case 4: xi[i_s] = *s =-quad5_1; break;
    }
    break;

  case P0_HEX:                   /* Constant discontinuous on hexahedron */
  case P1_HEX:                   /* linear discontinuous on hexahedron */
    switch( iquad ){
        case 0: xi[i_s] = *s =  0.; xi[i_t] = *t = 0.; break;
    }
    break;

  case LINEAR_BAR:
  case QUAD_BAR:
    // empty case so we don't get EH, switch on dim is enough
    break;

  case TRILINEAR_HEX:                   /* trilinear hexahedron */
  case C_TRILINEAR_HEX:                 /* trilinear hexahedron with additional centroid node */
    switch( iquad ){
        case 0: xi[i_s] = *s =  Ftemp1; xi[i_t] = *t =  Ftemp1; break;
        case 1: xi[i_s] = *s = -Ftemp1; xi[i_t] = *t =  Ftemp1; break;
        case 2: xi[i_s] = *s =  Ftemp1; xi[i_t] = *t = -Ftemp1; break;
        case 3: xi[i_s] = *s = -Ftemp1; xi[i_t] = *t = -Ftemp1; break;
    }
    break;

  case S_TRIQUAD_HEX:                  /* serendipity triquadratic hexahedron */
  case   TRIQUAD_HEX:                  /* triquadratric hexahedron */
    switch( iquad ){
        case 0: xi[i_s] = *s =  Ftemp2; xi[i_t] = *t =  Ftemp2; break;
        case 1: xi[i_s] = *s =     0.0; xi[i_t] = *t =  Ftemp2; break;
        case 2: xi[i_s] = *s = -Ftemp2; xi[i_t] = *t =  Ftemp2; break;
        case 3: xi[i_s] = *s =  Ftemp2; xi[i_t] = *t =     0.0; break;
        case 4: xi[i_s] = *s =     0.0; xi[i_t] = *t =     0.0; break;
        case 5: xi[i_s] = *s = -Ftemp2; xi[i_t] = *t =     0.0; break;
        case 6: xi[i_s] = *s =  Ftemp2; xi[i_t] = *t = -Ftemp2; break;
        case 7: xi[i_s] = *s =     0.0; xi[i_t] = *t = -Ftemp2; break;
        case 8: xi[i_s] = *s = -Ftemp2; xi[i_t] = *t = -Ftemp2; break;
    }
    break;

  case LINEAR_TET:
    
    /* 
     * I'm not sure where the original code came from, especially with
     * the "A" values.  But it doesn't seem to work the way I expect.  But
     * the new code is much simpler and makes more sense to me.  Using
     * a four-point integration rule for triangles here.
     * Scott A Roberts, 1514 - 2011-11-04
     */
    switch ( iquad ) {
    case 0: *s = one_third; *t = one_third; break;
    case 1: *s = 0.2; *t = 0.2; break;
    case 2: *s = 0.2; *t = 0.6; break;
    case 3: *s = 0.6; *t = 0.2; break;
    }
    double r = 1.0 - *s - *t;
    
    switch ( iside ) {
    case 1: xi[0] = *s; xi[1] = 0.; xi[2] = *t; break;
    case 2: xi[0] =  r; xi[1] = *s; xi[2] = *t; break;
    case 3: xi[0] = 0.; xi[1] = *t; xi[2] = *s; break;
    case 4: xi[0] = *t; xi[1] = *s; xi[2] = 0.; break;
    }

    break;

  default:
    EH(-1,"Unknown or unimplemented element type.\n");
    break;
  }

} /* END of routine find_surf_st */

/*****************************************************************************/

int
find_edge_s (const int iquad,                /* current GQ index  */
             const int ielem_type,           /* element type   */
             const int iedge,                /* current edge of element  */
             const int dim,                     /* dimensions of element    */
             double xi[DIM],           /* local GQ coordinates for surface 
                                           integral (these are returned)  */
             double *s  )

/*
*       Function which determines the Gaussian-quadrature points (s, t)
*       for surface integration from the given index, iquad, and element
*       type, ielem_type.
*
*        Author:          Scott Hutchinson (1421)
*        Date:            19 May 1992
*        Revised:         31 October 1996 RAC (converted to edge)
*
*        NOTE: this routine currently assumes element is a quad or a hex
*
*/

{
/* LOCAL VARIABLES */
  int i_s;                      /* index of first local coord */
                           /*  1.0 / sqrt (3.0)  */
  static const double Ftemp1 =  0.57735026918962584208;

                           /*  sqrt (3.0/5.0)   */
  static const double Ftemp2 =  0.77459666924148340428;

  i_s = -1;
  *s = 0.; /* to avoid uninitialized memory read later */


/* first calculate local coordinates on the n-dimensional surface of the element 
 * then translate them to the n+1 dimensional elemental coords */
  switch (dim){
  case 1:
      EH(-1,"cannot have edge for 1-D element");
    break;

  case 2:
    xi[2] = 0.;
    switch (iedge) {
    case 1:
      xi[0] = -1.;
      xi[1] = -1.;
      break;
    case 2:
      xi[0] =  1.;
      xi[1] = -1.;
      break;
    case 3:
      xi[0] =  1.;
      xi[1] =  1.;
      break;
    case 4:
      xi[0] = -1.;
      xi[1] =  1.;
      break;
    default:
      EH(-1,"Illegal edge number for 2-D element");
      break;
    }
    break;
  case 3:
    switch (iedge) {
    case 1:
      xi[0] = -1.;
      xi[2] = -1.;
      i_s = 1;
      break;
    case 2:
      xi[0] =  1.;
      xi[2] = -1.;
      i_s = 1;
      break;
    case 3:
      xi[1] = -1.;
      xi[2] = -1.;
      i_s = 0;
      break;
    case 4:
      xi[1] =  1.;
      xi[2] = -1.;
      i_s = 0;
      break;
    case 5:
      xi[0] = -1.;
      xi[2] =  1.;
      i_s = 1;
      break;
    case 6:
      xi[0] =  1.;
      xi[2] =  1.;
      i_s = 1;
      break;
    case 7:
      xi[1] = -1.;
      xi[2] =  1.;
      i_s = 0;
      break;
    case 8:
      xi[1] =  1.;
      xi[2] =  1.;
      i_s = 0;
      break;

    case 9:
      xi[0] = -1.;
      xi[1] = -1.;
      i_s = 2;
      break;
    case 10:
      xi[0] = -1.;
      xi[1] =  1.;
      i_s = 2;
      break;
    case 11:
      xi[0] =  1.;
      xi[1] = -1.;
      i_s = 2;
      break;
    case 12:
      xi[0] =  1.;
      xi[1] =  1.;
      i_s = 2;
      break;
    default:
      EH(-1,"Illegal edge number for 3-D element");
      break;
    }
    break;
  default:
    EH(-1,"Illegal element dimension");
    break;
  }

  switch( ielem_type ){                 /* select element */

  case BILINEAR_QUAD:                   /* bilinear quadrilateral */
  case C_BILINEAR_QUAD:              /* bilinear quadrilateral with additional centroid node */
  case S_BIQUAD_QUAD:                /* biquadratic serendipity quadrilateral */
  case   BIQUAD_QUAD:                /* biquadratic quadrilateral */
  case   BIQUAD_QUAD_LS:                /* biquadratic quadrilateral */

    /* do nothing */

    break;

  case TRILINEAR_HEX:                   /* trilinear hexahedron */
  case C_TRILINEAR_HEX:                 /* trilinear hexahedron with additional centroid node */
    switch( iquad ){
        case 0: xi[i_s] = *s =  Ftemp1; break;
        case 1: xi[i_s] = *s = -Ftemp1; break;
    }
    break;

  case S_TRIQUAD_HEX:                  /* serendipity triquadratic hexahedron */
  case   TRIQUAD_HEX:                  /* triquadratric hexahedron */
    switch( iquad ){
        case 0: xi[i_s] = *s =  Ftemp2;  break;
        case 1: xi[i_s] = *s =     0.0;  break;
        case 2: xi[i_s] = *s = -Ftemp2;  break;
    }
    break;

  default:
    EH(-1,"Unknown or unimplemented element type.\n");
    break;
  }

  return i_s;

} /* END of routine find_edge_s */
/*****************************************************************************/

void
find_surf_center_st (
                     const int ielem_type,      /* element type */
                     const int iside,           /* current side of element */
                     const int dim,             /* dimensions of element */
                     double xi[DIM],      /* (returned) local GQ coordinates for surface integral */
                     double *s,           /* Gaussian-quadrature points (s, t) */
                     double *t )
/*
*       Function which determines surface centroid (s, t)
*       element type, ielem_type.
*
*        Author:          Scott Hutchinson (1421)
*        Date:            19 May 1992
*        Revised:         24 May 1995 RAC
*
*        NOTE: this routine currently assumes element is a quad or a hex
*
*  TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98.
*/

{
  int i_s = -1;                 /* index of first local coord */
  int i_t = -1;                 /* index of second local coord */
                           /*  1.0 / sqrt (3.0)  */
  *t = 0.; /* to avoid uninitialized memory read later */
  static double one_third = 0.333333333333333;
  switch ( ielem_type )
  {
	  
	  case  BILINEAR_QUAD:          /* bilinear quadrilateral */
          case  BILINEAR_SHELL:                /*Bilinear shell */
	  case C_BILINEAR_QUAD:              /* bilinear quadrilateral with additional centroid node */
	  case S_BIQUAD_QUAD:                /* biquadratic serendipity quadrilateral */
	  case   BIQUAD_QUAD:                /* biquadratic quadrilateral */
          case BIQUAD_SHELL:                 /* biquadratic shell */
	  case   BIQUAD_QUAD_LS:             /* biquadratic quadrilateral for level set */
	  case TRILINEAR_HEX:                /* trilinear hexahedron */
	  case C_TRILINEAR_HEX:              /* trilinear hexahedron with additional centroid node */
	  case S_TRIQUAD_HEX:                /* serendipity triquadratic hexahedron */
	  case   TRIQUAD_HEX:                /* triquadratric hexahedron */
		  
		  /* first calculate local coordinates on the n-dimensional surface of the element 
		  * then translate them to the n+1 dimensional elemental coords */
		  switch (dim){
			  case 1:
				  xi[1] = 0.;
				  xi[2] = 0.;
				  switch(iside) {
					  case 1:
						  xi[0] = -1.;
						  break;
					  case 2:
						  xi[0] = 1.;
						  break;
					  default:
						  EH(-1,"Illegal side number for 1-D element");
						  break;
				  }
					  break;
			  case 2:
				  xi[2] = 0.;
				  switch (iside) {
					  case 1:
						  xi[1] = -1.;
						  i_s = 0;
						  break;
					  case 2:
						  xi[0] = 1.;
						  i_s = 1;
						  break;
					  case 3:
						  xi[1] = 1.;
						  i_s = 0;
						  break;
					  case 4:
						  xi[0] = -1.;
						  i_s = 1;
						  break;
					  default:
						  EH(-1,"Illegal side number for 2-D element");
						  break;
				  }
					  break;
			  case 3:
				  switch (iside) {
					  case 1:
						  xi[1] = -1.;
						  i_s = 0;
						  i_t = 2;
						  break;
					  case 2:
						  xi[0] =  1.;
						  i_s = 1;
						  i_t = 2;
						  break;
					  case 3:
						  xi[1] =  1.;
						  i_s = 0;
						  i_t = 2;
						  break;
					  case 4:
						  xi[0] = -1.;
						  i_s = 1;
						  i_t = 2;
						  break;
					  case 5:
						  xi[2] = -1.;
						  i_s = 0;
						  i_t = 1;
						  break;
					  case 6:
						  xi[2] =  1.;
						  i_s = 0;
						  i_t = 1;
						  break;
					  default:
						  EH(-1,"Illegal side number for 3-D element");
						  break;
				  }
				  break;
			  default:
				  EH(-1,"Illegal element dimension");
				  break;
		  }
		  break;
		  
	case LINEAR_TRI:
	case QUAD_TRI:   
	case QUAD6_TRI:
        case LINEAR_TET:
        case BILINEAR_TRISHELL:
	switch (dim)
	  {
	  case 1:
	    xi[1] = 0.;
	    xi[2] = 0.;
	    switch(iside) {
	    case 1:
	      xi[0] = 0.;
	      break;
	    case 2:
	      xi[0] = 1.;
	      break;
	    default:
	      EH(-1,"Illegal side number for 1-D element");
	      break;
	    }
	    break;
	    
	  case 2:
	    xi[2] = 0.;
	    switch (iside) {
	    case 1: 
	      xi[0] = *s = 0.5; xi[1] = 0.5; break;
	    case 2:
	      xi[0] = 0.0; xi[1] = *s = 0.5; break;
	    case 3: 
	      xi[0] = *s = 0.5; xi[1] = 0.0; break;
	    default:
	      EH(-1,"Illegal side number for 2D triangle \n");
	      break;
	    }
	    break;
	  case 3:
	    *s = *t = 0.5;
	    switch (iside) {
	    case 1: xi[0] = one_third; xi[1] = 0.0; xi[2] = one_third; break;
	    case 2: xi[0] = one_third; xi[1] = one_third; xi[2] = one_third; break;
	    case 3: xi[0] = 0.0; xi[1] = one_third; xi[2] = one_third; break;
	    case 4: xi[0] = one_third; xi[1] = one_third; xi[2] = 0.0; break;
	    }
	    break;
	  }
	break;
	default:
	break;
  }
  switch( ielem_type ){                 /* select element */

  case BILINEAR_QUAD:                   /* bilinear quadrilateral */
  case BILINEAR_SHELL:               /* bilinear shell */
  case C_BILINEAR_QUAD:              /* bilinear quadrilateral with additional centroid node */
        xi[i_s] = *s = 0; break;

  case S_BIQUAD_QUAD:                /* biquadratic serendipity quadrilateral */
  case   BIQUAD_QUAD:                /* biquadratic quadrilateral */
  case  BIQUAD_SHELL:                /* biquadratic shell */
  case   BIQUAD_QUAD_LS:             /* biquadratic quadrilateral for level set */
        xi[i_s] = *s = 0; break;

  case TRILINEAR_HEX:                   /* trilinear hexahedron */
  case C_TRILINEAR_HEX:                 /* trilinear hexahedron with additional centroid node */
        xi[i_s] = *s =  0; xi[i_t] = *t =  0; break;

  case S_TRIQUAD_HEX:                  /* serendipity triquadratic hexahedron */
  case   TRIQUAD_HEX:                  /* triquadratric hexahedron */
        xi[i_s] = *s =  0; xi[i_t] = *t =  0; break;

  case LINEAR_TRI:
  case QUAD_TRI:   
  case QUAD6_TRI:
  case LINEAR_TET:
  case BILINEAR_TRISHELL:
	break;
  default:
    EH(-1,"Unknown or unimplemented element type.\n");
    break;
  }

} /* END of routine find_surf_center_st */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
find_nodal_stu (const int inode,           /* current node index */
                const int ielem_type,      /* element type */
                double  *s,          /* local GQ coordinates  */
                double  *t,          /*     (returned    */
                double  *u )         /*      values )    */

    /************************************************************************
     *
     * find_nodal_stu()
     *
     *       Function which determines the Gaussian-quadrature points
     *       s, t and u, at a given local node.
     *
     * Input
     * ------
     *   inode    = Local node number of the node at which the values for
     *              the local element coordinates will be evaluated.
     *   ielem_type = element type
     *
     * Output
     * --------
     *   s, t, u = values of the local element coordinates.
     *
     *        Author:          Scott Hutchinson (1421)
     *        Date:            19 May 1992
     *        Revised:         16 May 1995 Richard A. Cairncross
     *        Revised:          6 May 2002 Patrick K. Notz
     *        Revised:         24 August 2009 Peter R. Schunk
     *
     * Note: Triangles are not supported yet in this routine.
     ************************************************************************/
{
  /*
   * Branch according to the element types
   */
  switch (ielem_type) {

  case LINEAR_BAR:
    *t = 0.0;
    *u = 0.0;
    switch (inode){
    case 0:
      *s = -1;
      break;
    case 1:
      *s = 1;
      break;
    default:
      EH(-1, "Trying to get nodal local stu for BILINEAR at illegal node");
    }
    break;

  case QUAD_BAR:
    *t = 0.0;
    *u = 0.0;
    switch (inode){
    case 0:
      *s = -1;
      break;
    case 1:
      *s = 1;
      break;
    case 2:
      *s = 0;
      break;
    default:
      EH(-1, "Trying to get nodal local stu for BILINEAR at illegal node");
    }
    break;

  case LINEAR_TRI:   /* linear triangle */
  case BILINEAR_TRISHELL:  /* linear triangular shell */
    *u = 0.0;
    switch (inode) {
    case 0:
      *s = 1.0;
      *t = 0.0;
      break;
    case 1:
      *s =  0.0;
      *t =  1.0;
      break;
    case 2:
      *s =  0.0;
      *t =  0.0;
      break;
    default:
      EH(-1, "Trying to get nodal local stu for LINEAR_TRI at illegal node");
    }
    break;

   case QUAD_TRI:   /* quadratic triangle */
   case QUAD6_TRI:  /* quadratic triangle with 6th order quadrature */
    *u = 0.0;
    switch (inode) {
    case 0:
      *s = 1.0;
      *t = 0.0;
      break;
    case 1:
      *s =  0.0;
      *t =  1.0;
      break;
    case 2:
      *s =  0.0;
      *t =  0.0;
      break;
    case 3:
      *s = 0.5;
      *t = 0.5;
      break;
    case 4:
      *s =  0.0;
      *t =  0.5;
      break;
    case 5:
      *s =  0.5;
      *t =  0.0;
      break;
    default:
      EH(-1, "Trying to get nodal local stu for QUAD_TRI at illegal node");
    }
    break;
        
  case BILINEAR_QUAD:   /* bilinear quadrilateral */
  case BILINEAR_SHELL:
  case C_BILINEAR_QUAD: /* bilinear quadrilateral 
                         * with additional centroid node */
    *u = 0.0;
    switch (inode) {
    case 0:
      *s = -1.0;
      *t = -1.0;
      break;
    case 1:
      *s =  1.0;
      *t = -1.0;
      break;
    case 2:
      *s =  1.0;
      *t =  1.0;
      break;
    case 3:
      *s = -1.0;
      *t =  1.0;
      break;
    case 4:
      *s = 0.0;
      *t = 0.0;
      break;
    default:
      EH(-1, "Trying to get nodal local stu for BILINEAR at illegal node");
    }
    break;

  case   BIQUAD_QUAD:                     /* biquadratic quadrilateral */
  case   BIQUAD_SHELL:
  case   BIQUAD_QUAD_LS:                  /* biquadratic quadrilateral for level set */
      *u = 0.0;
    switch( inode ){
    case 0:
      *s = -1.0;
      *t = -1.0;
      break;
    case 1:
      *s =  1.0;
      *t = -1.0;
      break;
    case 2:
      *s =  1.0;
      *t =  1.0;
      break;
    case 3:
      *s = -1.0;
      *t =  1.0;
      break;
    case 4:
      *s =  0.0;
      *t = -1.0;
      break;
    case 5:
      *s =  1.0;
      *t =  0.0;
      break;
    case 6:
      *s =  0.0;
      *t =  1.0;
      break;
    case 7:
      *s = -1.0;
      *t =  0.0;
      break;
    case 8:
      *s =  0.0;
      *t =  0.0;
      break;
    default:
      EH(-1, "Trying to get nodal local stu for BIQUAD at illegal node");
    }
    break;
  case LINEAR_TET:  /* trilinear tetrahedron */
    switch ( inode ) {
    case 0:
      *s = *t = *u = 0.;
      break;
    case 1:
      *s = 1.;
      *t = *u = 0.;
      break;
    case 2:
      *t = 1.;
      *s = *u = 0.;
      break;
    case 3:
      *u = 1.;
      *s = *t = 0.;
      break;
    }
    break;

  case TRILINEAR_HEX:  /* trilinear hexahedron */
  case C_TRILINEAR_HEX:/* trilinear hexahedron with additional centroid node */
    switch( inode ){
    case 0:
      *s = -1.0;
      *t = -1.0;
      *u = -1.0;
      break;
    case 1:
      *s =  1.0;
      *t = -1.0;
      *u = -1.0;
      break;
    case 2:
      *s =  1.0;
      *t =  1.0;
      *u = -1.0;
      break;
    case 3:
      *s = -1.0;
      *t =  1.0;
      *u = -1.0;
      break;
    case 4:
      *s = -1.0;
      *t = -1.0;
      *u =  1.0;
      break;
    case 5:
      *s =  1.0;
      *t = -1.0;
      *u =  1.0;
      break;
    case 6:
      *s =  1.0;
      *t =  1.0;
      *u =  1.0;
      break;
    case 7:
      *s = -1.0;
      *t =  1.0;
      *u =  1.0;
      break;
    case 8:
      *s =  0.0;
      *t =  0.0;
      *u =  0.0;
      break;
    default:
      EH(-1, 
         "Trying to get nodal local stu for TRILINEAR at illegal node\n");
      break;
    }
    break;

  case TRIQUAD_HEX:    /* triquadratric hexahedron */    
    switch (inode) {
    case 0:
      *s = -1.0;
      *t = -1.0;
      *u = -1.0;
      break;
    case 1:
      *s =  1.0;
      *t = -1.0;
      *u = -1.0;
      break;
    case 2:
      *s =  1.0;
      *t =  1.0;
      *u = -1.0;
      break;
    case 3:
      *s = -1.0;
      *t =  1.0;
      *u = -1.0;
      break;
    case 4:
      *s = -1.0;
      *t = -1.0;
      *u =  1.0;
      break;
    case 5:
      *s =  1.0;
      *t = -1.0;
      *u =  1.0;
      break;
    case 6:
      *s =  1.0;
      *t =  1.0;
      *u =  1.0;
      break;
    case 7:
      *s = -1.0;
      *t =  1.0;
      *u =  1.0;
      break;
    case 8:
      *s =  0.0;
      *t = -1.0;
      *u = -1.0;
      break;
    case 9:
      *s =  1.0;
      *t =  0.0;
      *u = -1.0;
      break;
    case 10:
      *s =  0.0;
      *t =  1.0;
      *u = -1.0;
      break;
    case 11:
      *s = -1.0;
      *t =  0.0;
      *u = -1.0;
      break;
    case 12:
      *s = -1.0;
      *t = -1.0;
      *u =  0.0;
      break;
    case 13:
      *s =  1.0;
      *t = -1.0;
      *u =  0.0;
      break;
    case 14:
      *s =  1.0;
      *t =  1.0;
      *u =  0.0;
      break;
    case 15:
      *s = -1.0;
      *t =  1.0;
      *u =  0.0;
      break;
    case 16:
      *s =  0.0;
      *t = -1.0;
      *u =  1.0;
      break;
    case 17:
      *s =  1.0;
      *t =  0.0;
      *u =  1.0;
      break;
    case 18:
      *s =  0.0;
      *t =  1.0;
      *u =  1.0;
      break;
    case 19:
      *s = -1.0;
      *t =  0.0;
      *u =  1.0;
      break;
    case 20:
      *s =  0.0;
      *t =  0.0;
      *u =  0.0;
      break;
    case 21:
      *s =  0.0;
      *t =  0.0;
      *u = -1.0;
      break;
    case 22:
      *s =  0.0;
      *t =  0.0;
      *u =  1.0;
      break;
    case 23:
      *s = -1.0;
      *t =  0.0;
      *u =  0.0;
      break;
    case 24:
      *s =  1.0;
      *t =  0.0;
      *u =  0.0;
      break;
    case 25:
      *s =  0.0;
      *t = -1.0;
      *u =  0.0;
      break;
    case 26:
      *s =  0.0;
      *t =  1.0;
      *u =  0.0;
      break;
    default:
      EH(-1, 
         "Trying to get nodal local stu for TRIQUAD_HEX at illegal node\n");
      break;
    }
    break;

  default:
    EH(-1,"Unknown or unimplemented element type in find_nodal_stu.\n");
    break;
  }
} /* END of routine find_nodal_stu  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double
Gq_weight(const int iquad,               /* current GQ index */
          const int ielem_type )         /* element type     */

     /*
      *       Function which determines the Gaussian-quadrature weight
      *       for a given GQ index.
      *
      *        Author:          Scott Hutchinson (1421)
      *        Date:            19 May 1992
      *        Revised:         4 June 1992
      *
      */
{
  static const double ftemp1 = 0.55555555555555555556;
  static const double ftemp2 = 0.88888888888888888888;

  static const double tri1 =   0.33333333333333333333;
  static const double tri2 =   0.109951743655322;
  static const double tri3 =   0.223381589678011;
  //static const double tri4 =  -0.5625;
  //static const double tri5 =   0.52083333333333333333;

  static const double tri6 =   0.050844906370207;
  static const double tri7 =   0.116786275726379;
  static const double tri8 =   0.082851075618374;

  static const double quad5_1 = 0.23692688505618908751;
  static const double quad5_2 = 0.47862867049936646804;
  static const double quad5_3 = 0.56888888888888888889;

  double weight = 0.0, weight_s, weight_t, weight_u;

  switch( ielem_type ){                 /* select element */

  case LINEAR_BAR:
    weight = 1.0;
    break;
    
  case QUAD_BAR:
    if (iquad == 0)
      weight = ftemp1;
    else if (iquad == 1)
      weight = ftemp2;
    else
      weight = ftemp1;
    break;

  case LINEAR_TRI:                      /* linear triangle */
  case BILINEAR_TRISHELL:               /* linear triangular shell */
    weight = tri1/2.0;
    break;
  case QUAD_TRI:                        /* quadratic triangle */
    /* 6 point rule */
    if ( iquad < 3)
      weight = tri2/2.0;
    else
      weight = tri3/2.0;
    /* 4 point rule
    if ( iquad == 0)
      weight = tri4;
    else
      weight = tri5;
    */
    break;
  case QUAD6_TRI:                        /* quadratic triangle with 6th order quadrature */
    /* 12 point rule */
    if ( iquad < 3)
      weight = tri6/2.0;
    else if ( iquad < 6)
      weight = tri7/2.0;
    else
      weight = tri8/2.0;
    break;

  case P0_QUAD:                         /* constant discontinuous on quadrilateral */
  case P0_SHELL:                        /* constant discontinuous on shell */
  case P1_QUAD:                         /* linear discontinuous on quadrilateral */
  case P1_SHELL:                        /* linear discontinuous on shell */
  case P0_HEX:                          /* constant discontinuous on hexahedron */
  case BILINEAR_QUAD:                   /* bilinear quadrilateral */
  case BILINEAR_SHELL:                  /* Bilinear shell */
  case C_BILINEAR_QUAD:                 /* bilinear quadrilateral with additional centroid node */
    weight = 1.0;
    break;

  case S_BIQUAD_QUAD:                   /* biquadratic serendipity quadrilateral */
  case   BIQUAD_QUAD:                   /* biquadratic quadrilateral */
  case   BIQUAD_SHELL:                  /* biquadratic shell */

    if (iquad%3 == 0)
      weight_s = ftemp1;
    else if ( (iquad-1)%3 == 0)
      weight_s = ftemp2;
    else
      weight_s = ftemp1;

    if (iquad < 3)
      weight_t = ftemp1;
    else if (iquad < 6)
      weight_t = ftemp2;
    else
      weight_t = ftemp1;

    weight = weight_s*weight_t;
    break;

  case   BIQUAD_QUAD_LS:                     /* biquadratic quadrilateral for level set */
    if (iquad%5 == 0)
      weight_s = quad5_1;
    else if ( (iquad-1)%5 == 0)
      weight_s = quad5_2;
    else if ( (iquad-2)%5 == 0)
      weight_s = quad5_3;
    else if ( (iquad-3)%5 == 0)
      weight_s = quad5_2;
    else
      weight_s = quad5_1;

    if (iquad < 5)
      weight_t = quad5_1;
    else if (iquad < 10)
      weight_t = quad5_2;
    else if (iquad < 15)
      weight_t = quad5_3;
    else if (iquad < 20)
      weight_t = quad5_2;
    else
      weight_t = quad5_1;

    weight = weight_s*weight_t;
    break;

  case LINEAR_TET:
    //weight = 1.0/6.0;
    weight = 1.0/24.0;
    /*
    if (iquad == 0)
      weight = wltet_1;
    else 
      weight = wltet_2;
    */
    break;


  case TRILINEAR_HEX:                   /* trilinear hexahedron */
  case P1_HEX:                          /* linear discontinuous on hexahedron */
  case C_TRILINEAR_HEX:                 /* trilinear hexahedron with additional centroid node */
    weight = 1.0;
    break;

  case S_TRIQUAD_HEX:                 /* triquadratic serendipity hexahedron */
  case   TRIQUAD_HEX:                     /* triquadratic hexahedron */
    if (iquad%3 == 0)
      weight_s = ftemp1;
    else if ( (iquad-1)%3 == 0)
      weight_s = ftemp2;
    else
      weight_s = ftemp1;
    
    if (iquad - 9*(iquad/9) < 3)
      weight_t = ftemp1;
    else if (iquad -9*(iquad/9) < 6)
      weight_t = ftemp2;
    else
      weight_t = ftemp1;

    if (iquad < 9)
      weight_u = ftemp1;
    else if (iquad < 18)
      weight_u = ftemp2;
    else
      weight_u = ftemp1;

    weight = weight_s*weight_t*weight_u;
    break;

  default:
    EH(-1,"Unknown or unimplemented element type.\n");
    break;
  }

  return weight;

} /* END of routine Gq_weight */
/*****************************************************************************/

double
Gq_surf_weight(const int iquad,               /* current GQ index  */
               const int ielem_type )         /* element type      */

/*
*       Function which determines the Gaussian-quadrature weight
*       for surface integration, for a given element type
*
*        Author:          Scott Hutchinson (1421)
*        Date:            19 May 1992
*        Revised:         4 June 1992
*
*/

{
  static const double five_ninths  = 0.55555555555555555556;
  static const double eight_ninths = 0.88888888888888888888;

  static const double quad5_1 = 0.23692688505618908751;
  static const double quad5_2 = 0.47862867049936646804;
  static const double quad5_3 = 0.56888888888888888889;

  double weight = 1e12;

  switch( ielem_type ){                 /* select element */

  case LINEAR_BAR:
  case QUAD_BAR:
    weight = 1.0; /* Note that this is a vertex */
    break;

  case BILINEAR_SHELL:
    weight = 1.0; /* Note that this is a curve, probably needs more detail */
    break;
    
  case BIQUAD_SHELL:
    weight = 1.0; /* Note that this is a curve, probably needs more detail */
    break;

  case LINEAR_TRI:                      /* linear triangle */
  case BILINEAR_TRISHELL:               /* linear triangular shell */
    weight = 0.5;
    break;

  case QUAD_TRI:                        /* quadratic triangle */
   switch ( iquad ) {
        case 0: weight = 0.5*five_ninths;  break;
        case 1: weight = 0.5*eight_ninths;  break;
        case 2: weight = 0.5*five_ninths;  break;
    }
    break;

  case QUAD6_TRI:                       /* quadratic triangle with 6th order quadrature */
   switch ( iquad ) {
        case 0: weight = 0.5*quad5_1;  break;
        case 1: weight = 0.5*quad5_2;  break;
        case 2: weight = 0.5*quad5_3;  break;
        case 3: weight = 0.5*quad5_2;  break;
        case 4: weight = 0.5*quad5_1;  break;
    }
    break;

  case P0_QUAD:                         /* constant discontinuous on quadrilateral */
  case P0_SHELL:                        /* constant discontinuous on shell */
  case P0_HEX:                          /* constant discontinuous on hexahedron */
  case P1_QUAD:                         /* linear discontinuous on quadrilateral */
  case P1_SHELL:                        /* linear discontinuous on shell */
  case P1_HEX:                          /* linear discontinuous on hexahedron */
  case BILINEAR_QUAD:                   /* bilinear quadrilateral */
  case C_BILINEAR_QUAD:                 /* bilinear quadrilateral with additional centroid node */
    weight = 1.0;
    break;

  case S_BIQUAD_QUAD:                /* biquadratic serendipity quadrilateral */
  case   BIQUAD_QUAD:                /* biquadratic quadrilateral */
   switch ( iquad ) {
        case 0: weight = five_ninths;  break;
        case 1: weight = eight_ninths;  break;
        case 2: weight = five_ninths;  break;
    }
    break;

  case   BIQUAD_QUAD_LS:             /* biquadratic quadrilateral for level set */
   switch ( iquad ) {
        case 0: weight = quad5_1;  break;
        case 1: weight = quad5_2;  break;
        case 2: weight = quad5_3;  break;
        case 3: weight = quad5_2;  break;
        case 4: weight = quad5_1;  break;
    }
    break;

  case LINEAR_TET:
    switch ( iquad ) {
    case 0: weight = -0.281250000000000; break;
    case 1: weight =  0.260416666666667; break;
    case 2: weight =  0.260416666666667; break;
    case 3: weight =  0.260416666666667; break;
    }
    break;
    
  case TRILINEAR_HEX:                   /* trilinear hexahedron */
  case C_TRILINEAR_HEX:                 /* trilinear hexahedron with additional centroid node */
    weight = 1.0;
    break;

  case S_TRIQUAD_HEX:                  /* triquadratic serendipity hexahedron */
  case   TRIQUAD_HEX:                     /* triquadratic hexahedron */
    switch ( iquad ) {
        case 0: weight = five_ninths  * five_ninths;  break;
        case 1: weight = eight_ninths * five_ninths;  break;
        case 2: weight = five_ninths  * five_ninths;  break;
        case 3: weight = five_ninths  * eight_ninths;  break;
        case 4: weight = eight_ninths * eight_ninths;  break;
        case 5: weight = five_ninths  * eight_ninths;  break;
        case 6: weight = five_ninths  * five_ninths;  break;
        case 7: weight = eight_ninths * five_ninths;  break;
        case 8: weight = five_ninths  * five_ninths;  break;
    }
    break;

  default:
    EH(-1,"Unknown or unimplemented element type.\n");
    break;
  }

  return weight;

} /* END of Gq_surf_weight  */
/*****************************************************************************/

double
Gq_edge_weight ( int iquad,                  /* current GQ index  */
                 int ielem_type )            /* element type      */

/*
*       Function which determines the Gaussian-quadrature weight
*       for surface integration, for a given element type
*
*        Author:          Scott Hutchinson (1421)
*        Date:            19 May 1992
*        Revised:         4 June 1992
*        Revised:         31 October 1996  RAC -> edge
*
*/

{
                        /*  = 5.0 / 9.0       */
  static const double ftemp1 = 0.55555555555555555556;

                        /*  = 8.0 / 9.0       */
  static const double ftemp2 = 0.88888888888888888888;

  double weight = 1e12;

  switch( ielem_type ){                 /* select element */

  case LINEAR_BAR:
  case QUAD_BAR:
  case BILINEAR_SHELL:
  case BIQUAD_SHELL:
  case BILINEAR_TRISHELL:
    EH(-1,"Edges are undefined for BAR elements.\n");
    break;

  case BILINEAR_QUAD:                   /* bilinear quadrilateral */
  case C_BILINEAR_QUAD:              /* bilinear quadrilateral with additional centroid node */
  case S_BIQUAD_QUAD:                /* biquadratic serendipity quadrilateral */
  case   BIQUAD_QUAD:                /* biquadratic quadrilateral */
  case   BIQUAD_QUAD_LS:             /* biquadratic quadrilateral for level set */
    weight = 1.0;
    break;

  case TRILINEAR_HEX:                   /* trilinear hexahedron */
  case C_TRILINEAR_HEX:                 /* trilinear hexahedron with additional centroid node */
    weight = 1.0;
    break;

  case S_TRIQUAD_HEX:                  /* triquadratic serendipity hexahedron */
  case   TRIQUAD_HEX:                     /* triquadratic hexahedron */
    switch ( iquad ) {
        case 0: weight = ftemp1;  break;
        case 1: weight = ftemp2;  break;
        case 2: weight = ftemp1;  break;  
    }
    break;

  default:
    EH(-1,"Unknown or unimplemented element type.\n");
    break;
  }

  return weight;

} /* END of routine Gq_edge_weight */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int 
in_list(const int ivalue, const int i, const int iend, int *ivector)

/*
*        This function searches an integer vector, ivector[i:iend-1],
*       for the presence of a number, ivalue.  It returns the index of the 
*       value, or -1, if the number, ivalue, is not found in the list.
*   
*        The function is used (amongst other purposes) to see if a local node 
*       number is in the adjacency list of an element and to return 
*       its position if so.
*
*        Author:          Scott Hutchinson (1421)
*        Date:            15 May 1992
*        Revised:         26 May 1992
*
*        Revised          13 Feb 1998 , Thomas Baer (9112)
*/
{
  int i2 = i;
  if (iend <= i2) return -1;
  if (!ivector) return -1;
  while (ivalue != ivector[i2] && ++i2 < iend);
  if (i2 < iend) return i2;
  return -1;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
get_type(char string[],         /* EXODUS name of parent element  */
         const int nodes,       /* number of nodes in this element  */
         const int attrs)       /* number of atrributes in this element  */

     /*
      *       Function which returns the element type according to this analysis code
      *       based on the EXODUS element type string and the number of nodes in the
      *       element.
      *       The types are defined in the file el_elm.h.
      *
      *        Author:          Scott Hutchinson (1421)
      *        Date:            29 June 1992
      *        Revised:         29 June 1992
      *        Revised:         24 August 2009 PRS
      *        Revised:         30 November 2011 SAR - TRISHELL3
      */
{
  char err_msg[MAX_CHAR_ERR_MSG];
  int answer = -1;

  stringup(string);
 
  if (strncmp(string, "BAR", 3) == 0)
    {  /* select element shape */
      switch (nodes){              /* select number of nodes in this element */
      case 2:                      /* bilinear quadralateral */
	answer = LINEAR_BAR;
	break;
      case 3:                      /* bilinear quadrilateral with additional centroid node */
	answer = QUAD_BAR;
	break;
      default:
	sprintf(err_msg,"Bar element with %d nodes not implemented.\n", nodes);
	EH(-1,err_msg);
      }
    }

  else if (strncmp(string, "SHELL", 5) == 0 ) 
    {  /* select element shape */
      switch (nodes){              /* select number of nodes in this element */
      case 4:                      /* bilinear quadralateral */
	answer = BILINEAR_SHELL;
	break;
      case 9:                      
	answer = BIQUAD_SHELL;
	break;
      default:
	sprintf(err_msg,"Shell element with %d nodes not implemented.\n", nodes);
	EH(-1,err_msg);
      }
    }
  
  else if (strncmp(string, "QUAD", 4) == 0 )
    {  /* select element shape */
      switch (nodes){              /* select number of nodes in this element */
      case 4:                      /* bilinear quadralateral */
	answer = BILINEAR_QUAD;
	break;
      case 5:                      /* bilinear quadrilateral with additional centroid node */
	answer = C_BILINEAR_QUAD;
	break;
      case 8:                      /* serendipity biquadratic quadralateral */
	answer = S_BIQUAD_QUAD;
	break;
      case 9:                      /* biquadratic quadrilateral */
	answer = BIQUAD_QUAD; 
	break;
      default:
	sprintf(err_msg,"Quadrilateral element with %d nodes not implemented.\n", nodes);
	EH(-1,err_msg);
      }
    }

  else if (strncmp(string, "TRI", 3) == 0)
    {  /* select element shape */
      switch (nodes){              /* select number of nodes in this element */
      case 3:
	switch (attrs) {
	case 0:                    /* bilinear triangle */
	  answer = LINEAR_TRI;
	  break;
	case 1:                    /* bilinear triangular shell */
	  answer = BILINEAR_TRISHELL;
	  break;
	default:
	  sprintf(err_msg,"TRI/TRISHELL element with %d dimensions not implemented.\n", attrs);
	  EH(-1,err_msg);
	}
	break;
	case 6:
	  answer = QUAD_TRI;
	  break;
      default:
	sprintf(err_msg,"TRIANGLE element with %d nodes not implemented.\n", nodes);
	EH(-1,err_msg);
      }
    }

  else if (strncmp(string, "TETRA4", 6) == 0 || strncmp(string, "TETRA", 5) == 0)
    {  /* select element shape */
      switch (nodes){              /* select number of nodes in this element */
      case 4:                      /* bilinear quadralateral */
        answer = LINEAR_TET;
	break;
      default:
	sprintf(err_msg,"TET element with %d nodes not implemented.\n", nodes);
	EH(-1,err_msg);
      }
    }
  
  else if (strncmp(string, "HEX", 3) == 0)
    {  /* select element shape */
      switch (nodes){              /* select number of nodes in this element */
      case 8:                      /* trilinear hexahedron */
	answer = TRILINEAR_HEX;
	break;
      case 9:                      /* trilinear hexahedron with additional centroid node */
	answer = C_TRILINEAR_HEX;
	break;
      case 20:                      /* serendipity triquadratic hexahedron */
	answer = S_TRIQUAD_HEX;
	break;
      case 27:                      /* triquadratic hexahedron */
	answer = TRIQUAD_HEX;
	break;
      default:
	sprintf(err_msg,"Hexahedron element with %d nodes not implemented.\n", nodes);
	EH(-1,err_msg);
      }
    }
  else {
    sprintf(err_msg,"Element type %s not supported!\n", string);
    EH(-1,err_msg);
  }


  /* return desired element information */
  return answer;

} /* END of routine get_type */
/*****************************************************************************/
int
centroid_node( int elem_type )
{
  int elem_shape = type2shape(elem_type);

  switch ( elem_shape )
    {
    case LINE_SEGMENT:
      switch ( elem_type )
	{
	case QUAD_BAR:
	  return(2);
	default:
	  return(-1);
	}

    case SHELL:
      switch ( elem_type )
	{
	case BIQUAD_SHELL:
	  return(8);
	default:
	  return(-1);
	}

    case QUADRILATERAL:
      switch ( elem_type )
	{
	case C_BILINEAR_QUAD:
	  return(4);
        case BIQUAD_QUAD:
        case BIQUAD_QUAD_LS:
	  return(8);
	default:
	  return(-1);
	}
    case HEXAHEDRON:
      switch ( elem_type) 
	{
	case C_TRILINEAR_HEX:
	  return(8);
	case TRIQUAD_HEX:
	  return(20);
	default:
	  return(-1);
	}
    default:
      return(-1);
    }
}



int
load_surf_st( int ielem_type,
	      int id_side,
	      int dim,
	      double xi[DIM],
	      double s,
	      double t )

     /* load_surf_st
      * 
      *  Author: TAB
      *  Date  : August 2003
      *
      *  Description -
      *             This function performs the task of mapping a local surface coordinate to the corresponding
      *             point in the higher dimensional element to which that surface is attached.
      *             For example, for a two-dimesional element I can express a
      *             position on a side by a parameter that
      *             varies [-1,1], but if I want to evaluate variable values and shape function stuff at that point
      *             I must map that surface coordinate to the two-dimensional (s,t) 
      *             space of the element that the surfac
      *             is attached to.  This is a function of the elem_type and the value of id_side
      *  input - 
      *         elem_type - the type of element that the surface is a part of 
      *         id_side   - the surface of the element that the operations are being performed on
      *         dim       - the dimensional of the element (not the surface)
      *         s         - first surface coordinate ( the only one used for two-dimensional elementss
      *         t         - second surface coordinate ( used only with three dimensional elements )
      *
      *  output -
      *         xi        - array that contains the mapped element coordinates of the surface point
      */

{
  switch ( ielem_type )
    {
      /* all the quadrilateral elements */
    case P0_QUAD:
    case P0_SHELL:
    case P1_QUAD:
    case P1_SHELL:
    case BILINEAR_QUAD:
    case C_BILINEAR_QUAD:
    case S_BIQUAD_QUAD:
    case BIQUAD_QUAD:
    case BIQUAD_SHELL:
    case BILINEAR_SHELL:

      xi[2] = 0.0;

      switch ( id_side )
	{
	case 1:
	  xi[0] = s;
	  xi[1] = -1.0;
	  break;
	case 2:
	  xi[0] = 1.0;
	  xi[1] = s;
	  break;
	case 3:
	  xi[0] = s;
	  xi[1] = 1.0;
	  break;
	case 4:
	  xi[0] = -1.0;
	  xi[1] = s;
	  break;
	default:
	  EH(-1,"Illegal side number for 2D quadrilateral element.\n");
	  break;
	}
      break;

    case TRILINEAR_HEX:
    case C_TRILINEAR_HEX:
    case S_TRIQUAD_HEX:
    case TRIQUAD_HEX:

      switch( id_side)
	{
	case 1:
	  xi[0] = s;
	  xi[1] = -1.0;
	  xi[2] = t;
	  break;
	case 2:
	  xi[0] = 1.0;
	  xi[1] = s;
	  xi[2] = t;
	  break;
	case 3:
	  xi[0] = s;
	  xi[1] = 1.0;
	  xi[2] = t;
	  break;
	case 4:
	  xi[0] = -1.0;
	  xi[1] = s;
	  xi[2] = t;
	  break;
	case 5:
	  xi[0] = s;
	  xi[1] = t;
	  xi[2] = -1.0;
	  break;
	case 6:
	  xi[0] = s;
	  xi[1] = t;
	  xi[2] = 1.0;
	  break;
	}
      break;
    default:
      EH(-1,"Element type not yet implemented. \n");
      break;
    }
  return (0);
}

/* END of file el_elm_info.c  */
/*****************************************************************************/

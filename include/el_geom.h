/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

/*
 *$Id: el_geom.h,v 5.1 2007-09-18 18:53:41 prschun Exp $
 */

#ifndef GOMA_EL_GEOM_H
#define GOMA_EL_GEOM_H

/*------------------------------------------------------------------------------

  Include file containing element geometry variable declarations

        Author:          Scott Hutchinson (1421)
        Date:            24 November 1992
        Revised:	 8 January 1993
        Revised:         1997/08/01 09:38 MDT pasacki@sandia.gov

  Note: These EXODUS II variables are deprecated in favor of their counterparts
        in the master EXODUS II data structure 	Exo_DB *E in main.c,
        for example.
------------------------------------------------------------------------------*/

extern double **Coor; /* 2d dynamically allocated array containing
                         the physical space coordinates for each node.
                         The dimensions are Num_Dim by Num_Nodes */

extern int *Proc_Connect_Ptr; /* global connectivity ptrs for element nums */

extern int Num_Dim;  /* number of physical dimensions */
extern int Num_Node; /* total number of nodes in the mesh  */
extern int Num_Elem; /* total number of elements in the mesh */

extern int Max_NP_Elem; /* maximum number of nodes in any element   */

extern int Num_Internal_Elems;
/* Number of Elements on the local processor.
 This is equal to the number of "internal"
 plus the number of "border" elements.	      */

extern int *GNodes; /* data structure which contains the internal,
                       border and external nodes on each processor  */

extern int *GElems; /* Data structure which contains the internal
                       elements on each processor.  It is a map
                       from the local element number to the global
                       element number.
                       type: int vector of length Num_Internal_Elems*/

extern int Proc_Num_Elem_Blk;     /* number of element blocks in this processor   */
extern int *Proc_Num_Elem_In_Blk; /* number of elements in the processor's
                                   * element blocks */
extern int *Proc_Elem_Blk_Ids;    /* element block id's for the processor's element
                                    blocks                                       */
extern int *Proc_Elem_Blk_Types;  /* element block types for the processor's
                                     element blocks (integers, not charstrings)*/

extern int *Proc_Nodes_Per_Elem; /* # of nodes per element for each block on the
                                 current processor                            */

extern int *Proc_Num_Attr; /* # of attributes for each block on the current
                             processor                                    */

extern int *Proc_Elem_Connect; /* connectivity lists for the elements required
                                 by the current processor                     */

extern int Proc_Num_Node_Sets;       /* number node sets on current processor */
extern int Proc_NS_List_Length;      /* total length of all the node set lists*/
extern int *Proc_NS_Ids;             /* node sets ids for the node sets */
extern int *Proc_NS_Count;           /* node sets count record */
extern int *Proc_NS_Pointers;        /* node sets pointer record */
extern int *Proc_NS_List;            /* node sets list record */
extern double *Proc_NS_Dist_Fact;    /* node sets distribution factors */
extern int Proc_Num_Side_Sets;       /* number side sets on current processor */
extern int Proc_SS_Elem_List_Length; /* total length of all the side set elem
                                      * lists  */
extern int Proc_SS_Node_List_Length; /* total length of all the side set node
                                      * lists */
extern int *Proc_SS_Ids;             /* side sets ids */
extern int *Proc_SS_Elem_Count;      /* side sets count record for elements */
extern int *Proc_SS_Node_Count;      /* side sets count record for nodes */
extern int *Proc_SS_Elem_Pointers;   /* side sets pointer record for elements */
extern int *Proc_SS_Elem_List;       /* side sets element list record */

extern double *Proc_SS_Dist_Fact; /* side sets distribution factors */
extern int Num_Internal_Nodes;    /* that are owned exclusively by this proc */
extern int Num_Border_Nodes;      /* owned here but communicated elsewhere */
extern int Num_External_Nodes;    /* owned elsewhere, but needed here. */

#endif

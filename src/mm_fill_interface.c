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
 *$Id: mm_fill_interface.c,v 5.3 2010-03-18 23:47:45 hkmoffa Exp $
 */

#include <stdio.h>
#include <math.h>

#include "std.h"
#include "rf_vars_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "mm_species.h"
#include "rf_bc_const.h"
#include "mm_fill_jac.h"
#include "mm_interface.h"
#include "rd_mesh.h"
#include "mm_ns_bc.h"
#include "mm_eh.h"
#include "mm_qp_storage.h"
#include "el_elm.h"
#include "mm_elem_block_structs.h"
#include "mm_mp_const.h"
#include "mm_qtensor_model.h"

#define RGAS_CONST  8.314510E7  /* Gas Contant in g cm^2/(sec^2 g-mole K) */
#define RGAS_CALS   1.987093    /* Gas Constant in cal g-mole-1 K-1 */


/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
double
raoults_law_prxn(JACOBIAN_VAR_DESC_STRUCT *func_jac,
		 BOUNDARY_CONDITION_STRUCT *bc, int ip, 
		 ELEM_SIDE_BC_STRUCT *elem_side_bc,
		 double x_dot[MAX_PDIM], const double time, const double tt,
                 const double dt, const int intf_id)

    /**********************************************************************
     *
     * raoults_law_prxn():
     *
     *  Function which evaluates vapor-liquid equilibrium 
     * for adjacent phases using a discontinuous variable approach.
     * This function expects to be called from both sides of the interface
     * at the quadrature points on the interface.
     *
     * The expression returned is the value of - n dot j_k 
     *
     *  n dot (Y_k * rho * (v - v_s)) + W_k * S_k
     *
     * for both sides of the interface, depending upon which side of the
     * interface this function is being called from. 
     *		
     * Input
     * ----------
     *  wspec -> Species number on both sides of the interface whose
     *           concentrations are linked via a raoult's law-type 
     *           expression. ( bc->BC_Data_Int[0])
     *  eb_mat_liq = Element block id for the first phase -> this is 
     *               identified with the liquid phase  ( bc->BC_Data_Int[1])
     * 
     *  Output
     * --------
     * func[0] -> contains the contribution to the function from the
     *            current side of the interface.
     **********************************************************************/
{
  int i, num_terms, have_T, pos, speciesVT,
      index_is, index_lvdesc;
  int wspec = bc->BC_Data_Int[0];
  int ebID_mat_liq = bc->BC_Data_Int[1], ebIndex_liq, ebIndex_gas;
  int liquidSide = (Current_EB_ptr->Elem_Blk_Id == ebID_mat_liq);
  double func_value, massflux;
  double Y_wspec;
  MATRL_PROP_STRUCT *mp_liq, *mp_gas;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  INTERFACE_SOURCE_STRUCT *is, **is_hdl;

  if (af->Assemble_LSA_Mass_Matrix) return 0.0;

  /*
   * Figure out whether we have done this calculation before
   */
  ebIndex_liq = ebID_to_ebIndex(ebID_mat_liq);
  ebIndex_gas = ebID_to_ebIndex(bc->BC_Data_Int[2]);
  mp_liq = mp_glob[Matilda[ebIndex_liq]];
  mp_gas = mp_glob[Matilda[ebIndex_gas]];
  num_terms = mp_liq->Num_Species + mp_gas->Num_Species;
  have_T = (upd->vp[pg->imtrx][TEMPERATURE] != -1);
  if (have_T) num_terms++;

  /*
   * Find the storage for this bc from the linked list attached
   * to the element side. If it is not found, then malloc it.
   */
  is_hdl = (INTERFACE_SOURCE_STRUCT **)
      side_qp_storage_findalloc(VL_EQUIL_PRXN_BC, ip, elem_side_bc);
  if (*is_hdl == NULL) {
    *is_hdl = is_masstemp_create(mp_gas, mp_liq, TRUE);
  }
  is = *is_hdl;

  /*
   * Process the Source term 
   */
  if (! is[intf_id].Processed[wspec]) {
    is[intf_id].Processed[wspec] = TRUE;
    /*
     * Fill up the Var_Value[] list. Possibly change the type
     * of the species unknown vector at the same time.
     */
    is_masstemp_fillin(is, mp_gas, mp_liq, SPECIES_CONCENTRATION, time, intf_id);

    /*
     * Now, call the source routine and possibly do the
     * jacobian. The source term is defined as the
     * 
     */
    is[intf_id].Do_Jac = af->Assemble_Jacobian;
    source_vle_prxn(is, bc, mp_gas, mp_liq, have_T, intf_id);
    /*
     * We now need to possibly change the species variable type
     * for the source term and the dependent variable for the
     * Jacobian entries.
     */
    speciesVT = upd->Species_Var_Type;
    if (upd->Species_Var_Type == SPECIES_UNDEFINED_FORM) {
      speciesVT = SPECIES_MASS_FRACTION;
    }
    if (speciesVT != is[intf_id].SpeciesVT) {
      is_change1_speciesVT(is, wspec, 0, mp_gas, speciesVT, time, intf_id);
      is_change1_speciesVT(is, wspec, mp_gas->Num_Species, mp_liq,
			   speciesVT, time,intf_id);
      is_change1_speciesVT(is, mp_gas->Num_Species + wspec,
			   0, mp_gas, speciesVT, time, intf_id);
      is_change1_speciesVT(is, mp_gas->Num_Species + wspec,
			   mp_gas->Num_Species, mp_liq,
			   speciesVT, time, intf_id);
      convert_species_var(speciesVT, mp_gas, is[intf_id].SpeciesVT, 
			  is[intf_id].Var_Value, time);
      convert_species_var(speciesVT, mp_liq, is[intf_id].SpeciesVT, 
			  is[intf_id].Var_Value + mp_gas->Num_Species, time);
      is[intf_id].SpeciesVT = speciesVT;
      /*
       * Possibly Convert the format of what's held constant
       * during the partial derivatives wrt species variable
       * to one in which the sum of the mass fractions are
       * constant is a constraint
       */
      is_change1_lastspecies(is, wspec, 0, mp_gas, intf_id);
      is_change1_lastspecies(is, wspec, mp_gas->Num_Species, mp_liq, intf_id);
      is_change1_lastspecies(is, mp_gas->Num_Species + wspec,
			     0, mp_gas, intf_id);
      is_change1_lastspecies(is, mp_gas->Num_Species + wspec,
			     mp_gas->Num_Species, mp_liq, intf_id);
    }
  }

  /*
   *  First calculate the generic total mass flux to the surface
   *  term, and associated derivatives -> rho (v - vs) dot n
   *  for the current side of the interface that we are on.
   */
  massflux = mass_flux_surface(func_jac, x_dot, tt, dt);
  /*
   *  Next we need to adjust this expression by multiplication
   *  by Y_i
   */
  Y_wspec = mp->StateVector[SPECIES_UNK_0 + wspec];
  func_value = massflux * Y_wspec;

  /*
   * Residual Contribution
   *  -> is->SourceTerm[wspec] contains the source of 
   *     species k in the gas phase (mol cm-2 sec-1 if the
   *     species var type is SPECIES_CONCENTRATION and 
   *     gm cm-2 sec-1 if species var_type is SPECIES_MASS_FRACTION),
   *     which is the negative of the source term in the liquid
   *     phase.
   *     The source term in the liquid phase is located in 
   *           wspec + mp_gas->Num_Species
   */ 
  
  if (liquidSide) {
    pos = wspec + mp_gas->Num_Species;
    func_value += is[intf_id].SourceTerm[pos] * mp_liq->molecular_weight[wspec];
  } else {
    pos = wspec;
    func_value += is[intf_id].SourceTerm[pos] * mp_gas->molecular_weight[wspec];
  }

  /*
   * Jacobian terms
   */
  if (af->Assemble_Jacobian) {
  
    /*
     * Modify the existing entries in func_jac to account for the
     * extra Y_wspec multiplication
     */
    for (i = 0; i < func_jac->Num_lvdesc; i++) {
      func_jac->JacCol[i] *= Y_wspec;
    }
    for (i = 0; i < func_jac->Num_lvdof; i++) {
      func_jac->Jac_lvdof[i] *= Y_wspec;
    }
    /*
     * Add in the extra jacobian entry from the mass fraction
     * dependence
     */
    vd = is[intf_id].Var_List[pos];
    index_lvdesc = ei[pg->imtrx]->VDindex_to_Lvdesc[vd->List_Index];
    jacobianVD_addEntry(func_jac, index_lvdesc, massflux);

    /*
     * Make sure funcjac is big enough to accept all of the terms below
     */
    jacobianVD_realloc(&func_jac, is[intf_id].Num_Terms, 0);
    for (index_is = 0; index_is < is[intf_id].Num_Terms; index_is++) {
      /*
       * Look up what variable description this dependence is for
       */
      vd = is[intf_id].Var_List[index_is];
      /*
	* Use this vd structure to find the lvdesc index for
	* the corresponding variable description
	*
	* HKM -> HAVE TO GUARD AGAINST LAST SPECIES IN PHASE.
	*      -> solution is to create a vd entry for that
	*         species even if it is not part of the solution
	*         We may need it anyway for postprocessing.
	* 
	*/
      index_lvdesc = ei[pg->imtrx]->VDindex_to_Lvdesc[vd->List_Index];
      if (index_lvdesc >= 0 &&
	  DOUBLE_NONZERO(is[intf_id].JacMatrix[pos][index_is])) {
	/*
	* Transfer the pertinent information to the
	*  Jacobian_Var_Desc structure
	*/
	if (liquidSide) {
	  jacobianVD_addEntry(func_jac, index_lvdesc, 
			      is[intf_id].JacMatrix[pos][index_is] * 
			      mp_liq->molecular_weight[wspec]);
	} else {
	  jacobianVD_addEntry(func_jac, index_lvdesc,
			      is[intf_id].JacMatrix[pos][index_is] * 
			      mp_gas->molecular_weight[wspec]);
	}
      }
    }
  }
  return func_value;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void 
source_vle_prxn(INTERFACE_SOURCE_STRUCT *is, BOUNDARY_CONDITION_STRUCT *bc,
		MATRL_PROP_STRUCT *mp_gas, MATRL_PROP_STRUCT *mp_liq, 
		int have_T, const int intf_id)

    /**********************************************************************
    *
    * source_vle_prxn():
    *
    *  Implements Sl_k = k_f * ( Cg_k - K_k * Cl_k)
    *
    *  where K_k = P_v/ (RT C_l_mix) 
    *
    *  and stores the Jacobian as well. The dependent variable is assumed
    *  to be concentration. Sk_k is returned with units of mol cm-2 sec-1.
    *  Right now a source term for heat is not returned. However, it would
    *  be easy to implement.
    *  The source term vector contains only one entry, that for the 
    *  particular species at equilibrium with its vapor. However, the 
    *  column of the Jacobian corresponding to that species is dense,
    *  including the temperature dependence.
    *
    *  The following format of the InterfaceSource_Structure is assumed:
    *
    *     0 to mp_gas->Num_Species  -> gas phase concentrations
    *     mp_gas->Num_Species to plus mp_liq->Num_Secies
    *                               -> liquid phase concentrations
    *     T, iff have_T is true. 
    *
    *  Boundary condition structure is assumed to hold the following
    *
    *   BC_Data_Int[0] -> The species undergoing the equilibrium condition
    *   BC_Data_Float[0] -> The pseudo reaction rate constant.
    *   BC_Data_Int[1] -> The element block curresponding to the liquid
    *                     Phase, -> mp_liq is input from the parameter
    *                     list. The two should correspond to the same block
    *   BC_Data_Int[2] -> The element block corresponding to the gas
    *                     Phase, -> mp_gas is input from the parameter
    *                     list. The two should correspond to the same block.
    *  
    *  Fields Used in the Materials Property Structure
    *     mp_liq->VaporPressureModel[wspec] -> Get the vapor pressure model
    *     mp_liq->u_vapor_pressure[wspec] -> parameters for the model
    *     mp_liq->vapor_pressure[wspec] -> used if vp is constant.
    *
    ***********************************************************************/
{
  int k, wspec, wspec_liq, pos, index_T = 0;
  double C_l_mix, k_f, T;
  double psat_w, dpsatdT_w, *vv = is[intf_id].Var_Value, K_equil;
  double **Jac;

  /*
   * Zero the calculable fields in the interface structure
   */
  interface_source_zero(&is[intf_id]);

  /*
   * Species number that the boundary condition is being applied on
   */
  wspec = bc->BC_Data_Int[0];
  wspec_liq = wspec + mp_gas->Num_Species;

  /*
   * Find the vapor pressure
   */
  if        (mp_liq->VaporPressureModel[wspec] == ANTOINE) {
    antoine_psat(mp_liq->u_vapor_pressure[wspec], &psat_w, &dpsatdT_w);
  } else if (mp_liq->VaporPressureModel[wspec] == RIEDEL) {
    riedel_psat(mp_liq->u_vapor_pressure[wspec], &psat_w, &dpsatdT_w);
  } else {
    psat_w = mp_liq->vapor_pressure[wspec];
  }

  /*
   * Find the concentration in the liquid phase, C_l_mix
   */
  if (is[intf_id].SpeciesVT != SPECIES_CONCENTRATION) {
    EH(GOMA_ERROR, "unimplemented complication");
  }
  C_l_mix = 0.0;
  for (k = 0, pos = mp_gas->Num_Species; k < mp_liq->Num_Species; 
       k++, pos++) {
    C_l_mix += vv[pos];
  }

  /*
   * Find a value of the temperature, even if it is not part of the
   * input.
   */
  if (have_T) {
    index_T = mp_gas->Num_Species + mp_liq->Num_Species;
    T = vv[index_T];
  } else {
    T = mp_liq->reference[TEMPERATURE];
  } 

  /*
   *  Make sure that the interface solution structure understands
   *  that we have carried out the calculation assuming that the
   *  species comprising the unknowns in that structure are
   *  species concentrations, C_k, and that they have units
   *  of mol cm-3.
   */
  is[intf_id].SpeciesVT = SPECIES_CONCENTRATION;

  /*
   * Pseudo reaction rate constant (cm sec-1)
   */
  k_f = bc->BC_Data_Float[0];
  K_equil = psat_w / (RGAS_CONST * T * C_l_mix);

  is[intf_id].SourceTerm[wspec] = k_f * (K_equil * vv[wspec_liq] - vv[wspec]);
  is[intf_id].SourceTerm[wspec_liq] = -  is[intf_id].SourceTerm[wspec];

  if (is[intf_id].Do_Jac) {
    Jac = is[intf_id].JacMatrix;
    Jac[wspec][wspec] = - k_f;
    for (k = 0, pos = mp_gas->Num_Species; k < mp_liq->Num_Species;
	 k++, pos++) {
       Jac[wspec][pos] = - k_f * K_equil * vv[wspec_liq] / C_l_mix;
    }
    Jac[wspec][wspec_liq] +=  k_f * K_equil;
    if (have_T) {
      Jac[wspec][index_T] = 
	  (k_f * vv[wspec_liq] * 
	   (dpsatdT_w / (RGAS_CONST * T * C_l_mix) - K_equil/T));
    }
    for (k = 0; k < is[intf_id].Num_Terms; k++) {
      Jac[wspec_liq][k] = - Jac[wspec][k];
    }
  }
  is[intf_id].Processed[wspec] = TRUE;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double
is_equil_prxn(JACOBIAN_VAR_DESC_STRUCT *func_jac,
	      BOUNDARY_CONDITION_STRUCT *bc, int ip, 
	      ELEM_SIDE_BC_STRUCT *elem_side_bc,
	      double x_dot[MAX_PDIM], const double time,
              const double tt, const double dt, const int intf_id)

    /**********************************************************************
     *
     * is_equil_prxn():
     *
     * Function which evaluates  equilibrium 
     * for adjacent phases using a discontinuous variable approach, 
     * assuming that each phase is an ideal solution.
     * This function expects to be called from both sides of the interface
     * at the quadrature points on the interface.
     *
     * The expression returned is the value of - n dot j_k = 
     *
     *  n dot (Y_k * rho * (v - v_s)) + W_k * S_k
     *
     * for both sides of the interface, depending upon which side of the
     * interface this function is being called from. v is the mass
     * average velocity, while v_s is the velocity of the interface.
     *
     * for the SPECIES_CONCENTRATION and SPECIES_MOLE_FRACTION variable
     * types the equation returned is :
     *
     *  n dot j_k / M_k =   n dot ( C_k * (v - v_s)) +  S_k
     *		
     * Input
     * ----------
     *  wspec -> Species number on both sides of the interface whose
     *           concentrations are linked via a one-to-one 
     *           expression. ( bc->BC_Data_Int[0])
     *  eb_mat_a = Element block id for the first phase -> this is 
     *             identified as the "reactants"  ( bc->BC_Data_Int[1])
     *  eb_mat_b = Element block id for the second phase -> this is 
     *             identified as the "products"  ( bc->BC_Data_Int[2])
     * 
     *  Output
     * --------
     *  func[0] -> contains the contribution to the function from the
     *             current side of the interface.
     **********************************************************************/
{
  int i, num_terms, have_T, pos, speciesVT,
      index_is, index_lvdesc;
  int wspec_a = bc->BC_Data_Int[0], wspec_b;
  int ebID_mat_a = bc->BC_Data_Int[1], ebIndex_a, ebIndex_b;
  int a_Side = (Current_EB_ptr->Elem_Blk_Id == ebID_mat_a);
  double func_value = 0.0, volflux = 0.0, C_wspec = 0.0;
  MATRL_PROP_STRUCT *mp_a = 0, *mp_b = 0;
  INTERFACE_SOURCE_STRUCT *is, **is_hdl;

  VARIABLE_DESCRIPTION_STRUCT *vd;

  if (af->Assemble_LSA_Mass_Matrix) return 0.0;

  /*
   * Figure out whether we have done this calculation before
   */
  ebIndex_a = ebID_to_ebIndex(ebID_mat_a);
  ebIndex_b = ebID_to_ebIndex(bc->BC_Data_Int[2]);
  mp_a = mp_glob[Matilda[ebIndex_a]];
  mp_b = mp_glob[Matilda[ebIndex_b]];
  num_terms = mp_a->Num_Species + mp_b->Num_Species;
  have_T = (upd->vp[pg->imtrx][TEMPERATURE] != -1);
  if (have_T) num_terms++;

  /*
   * Currently we limit this bc to one type of species
   * var type. This will be relaxed in the future.
   */
  if (upd->Species_Var_Type != SPECIES_CONCENTRATION) {
    EH(GOMA_ERROR,"unimplemented");
  }

  /*
   * Find the storage for this bc from the linked list attached
   * to the element side. If it is not found, then malloc it.
   */
  is_hdl = (INTERFACE_SOURCE_STRUCT **)
      side_qp_storage_findalloc(IS_EQUIL_PRXN_BC, ip, elem_side_bc);
  if (*is_hdl == NULL) {
    *is_hdl = is_masstemp_create(mp_a, mp_b, TRUE);
  }
  is = *is_hdl;

  /*
   * Process the Source term if it hasn't been already
   */
  if (! is[intf_id].Processed[wspec_a]) {
    is[intf_id].Processed[wspec_a] = TRUE;
    /*
     * Fill up the Var_Value[] list. Possibly change the type
     * of the species unknown vector at the same time.
     */
    is_masstemp_fillin(is, mp_a, mp_b, SPECIES_CONCENTRATION, time, intf_id);

    /*
     * Now, call the source routine and possibly do the
     * jacobian. The source term is defined as having
     * units of mol cm-2 sec-1. The dependent variable is
     * the species concentration.
     */
    is[intf_id].Do_Jac = af->Assemble_Jacobian;
    source_is_equil_prxn(is, bc, mp_a, mp_b, have_T,intf_id);

    /*
     * We now need to possibly change the species variable type
     * for the source term and the dependent variable for the
     * Jacobian entries.
     */
    speciesVT = upd->Species_Var_Type;
    if (upd->Species_Var_Type == SPECIES_UNDEFINED_FORM) {
      speciesVT = SPECIES_CONCENTRATION;
    }

    wspec_b = mp_a->Num_Species + wspec_a;

    if (speciesVT != is[intf_id].SpeciesVT) {
      /*
       * Change the gas and liquid source term derivative wrt Ck to a
       * derivitive wrt the dependent variable form of the species
       * equations for a and b phases (4 in all).
       */
      is_change1_speciesVT(is, wspec_a, 0, mp_a, speciesVT, time, intf_id);
      is_change1_speciesVT(is, wspec_a, mp_a->Num_Species, mp_b,
			   speciesVT, time, intf_id);
      is_change1_speciesVT(is, wspec_b, 0, mp_a, speciesVT, time, intf_id);
      is_change1_speciesVT(is, wspec_b, mp_a->Num_Species, mp_b,
			   speciesVT, time, intf_id);
      /*
       * Convert the storred values of the dependent variables
       * in the is structure to their proper forms.
       */
      convert_species_var(speciesVT, mp_a, is[intf_id].SpeciesVT, 
			  is[intf_id].Var_Value, time);
      convert_species_var(speciesVT, mp_b, is[intf_id].SpeciesVT, 
			  is[intf_id].Var_Value + mp_a->Num_Species, time);
      /*
       * Ok update SpeciesVT -> we have finished changing the
       * species variable type in the is structure
       */
      is[intf_id].SpeciesVT = speciesVT;
      /*
       * Convert the format of what's held constant
       * during the partial derivatives wrt species variable
       * to one in which the (sum of the mass fractions are
       * constant) constraint replaces the (last species is
       * constant) constraint.
       */
      is_change1_lastspecies(is, wspec_a, 0, mp_a,intf_id);
      is_change1_lastspecies(is, wspec_a, mp_a->Num_Species, mp_b, intf_id);
      is_change1_lastspecies(is, wspec_b, 0, mp_a, intf_id);
      is_change1_lastspecies(is, wspec_b, mp_a->Num_Species, mp_b, intf_id);
    }
  }

  /*
   *  First calculate the generic total vol flux to the surface
   *  term, and associated derivatives ->(v - vs) dot n
   *  for the current side of the interface that we are on.
   */
  volflux = vol_flux_surface(func_jac, x_dot, tt, dt);

  /*
   *  Next we need to adjust this expression by multiplication
   *  by wspec concentration. Note, mp refers to the current
   *  side of the interface.
   */
  if (upd->Species_Var_Type == SPECIES_CONCENTRATION) {
    C_wspec = mp->StateVector[SPECIES_UNK_0 + wspec_a];
    func_value = volflux * C_wspec;
  } else {
    EH(GOMA_ERROR, "unimplemented");
  }

  /*
   * Residual Contribution
   *  -> is->SourceTerm[wspec] contains the source of 
   *     species k in the a phase (mol cm-2 sec-1 if the
   *     species var type is SPECIES_CONCENTRATION
   *
   *     The source term in the b phase is located in 
   *           wspec_a + mp_b->Num_Species
   */
  
  if (a_Side) {
    pos = wspec_a;
    func_value += is[intf_id].SourceTerm[pos];
  } else {
    pos = wspec_a + mp_a->Num_Species;
    func_value += is[intf_id].SourceTerm[pos];
  }

  /*
   * Jacobian terms
   */
  if (af->Assemble_Jacobian) {
  
    /*
     * Modify the existing entries in func_jac to account for the
     * extra C_wspec multiplication
     */
    for (i = 0; i < func_jac->Num_lvdesc; i++) {
      func_jac->JacCol[i] *= C_wspec;
    }
    for (i = 0; i < func_jac->Num_lvdof; i++) {
      func_jac->Jac_lvdof[i] *= C_wspec;
    }
    /*
     * Add in the extra jacobian entry from the concentration
     * dependence
     */
    vd = is[intf_id].Var_List[pos];
    index_lvdesc = ei[pg->imtrx]->VDindex_to_Lvdesc[vd->List_Index];
    jacobianVD_addEntry(func_jac, index_lvdesc, volflux);

    /*
     * Make sure funcjac is big enough to accept all of the terms below
     */
    jacobianVD_realloc(&func_jac, is[intf_id].Num_Terms, 0);
    for (index_is = 0; index_is < is[intf_id].Num_Terms; index_is++) {
      /*
       * Look up what variable description this dependence is for
       */
      vd = is[intf_id].Var_List[index_is];
      /*
	* Use this vd structure to find the lvdesc index for
	* the corresponding variable description
	*
	* HKM -> HAVE TO GUARD AGAINST LAST SPECIES IN PHASE.
	*      -> solution is to create a vd entry for that
	*         species even if it is not part of the solution
	*         We may need it anyway for postprocessing.
	* 
	*/
      index_lvdesc = ei[pg->imtrx]->VDindex_to_Lvdesc[vd->List_Index];
      if (index_lvdesc >= 0 &&
	  DOUBLE_NONZERO(is[intf_id].JacMatrix[pos][index_is])) {
	/*
	 * Transfer the pertinent information to the
	 *  Jacobian_Var_Desc structure
	 */
	jacobianVD_addEntry(func_jac, index_lvdesc, 
			    is[intf_id].JacMatrix[pos][index_is]);
      }
    }
  }
  return func_value;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void 
source_is_equil_prxn(INTERFACE_SOURCE_STRUCT *is, 
		     BOUNDARY_CONDITION_STRUCT *bc, 
		     MATRL_PROP_STRUCT *mp_a, MATRL_PROP_STRUCT *mp_b, 
		     int have_T, const int intf_id)

    /**********************************************************************
    *
    * source_is_equil_prxn():
    *
    *  Implements Sb_k = k_f * ( Ca_k - Cb_k / Kc_k)
    *
    *  where      Kc_k = Cb_mix / Ca_mix exp [ - Delta_Gk / RT ]
    *
    *             Delta_Gk = mu_b_k - mu_a_k 
    *
    *  and stores the Jacobian as well. The dependent variable is assumed
    *  to be concentration. Sb_k is returned with units of mol cm-2 sec-1.
    *  Right now a source term for heat is not returned. However, it would
    *  be easy to implement.
    *  The source term vector contains only one entry, that for the 
    *  particular species at equilibrium with its partner in the
    *  other phase. However, the 
    *  column of the Jacobian corresponding to that species is dense,
    *  including the temperature dependence.
    *
    *  The following format of the InterfaceSource_Structure is assumed:
    *
    *     0 to mp_a->Num_Species  -> material a concentrations
    *     mp_a->Num_Species to plus mp_b->Num_Secies
    *                               -> material b concentrations
    *     T, iff have_T is true.
    *
    *  Boundary condition structure is assumed to hold the following
    *
    *   BC_Data_Int[0] -> The species undergoing the equilibrium condition
    *   BC_Data_Float[0] -> The pseudo reaction rate constant, k_f
    *                            (cm/sec)
    *   BC_Data_Int[1] -> The element block curresponding to material a
    *   BC_Data_Int[2] -> The element block corresponding to material b
    *  
    *  Fields Used in the Materials Property Structure:
    *     mp->SSChemPot
    *     mp->PSChemPot
    *     mp->ChemPot
    ***********************************************************************/
{
  int k, wspec_a, wspec_b, kspec, pos, index_T = -1;
  double C_a_mix, C_b_mix, k_f, T, deltaG, Kc_equil;
  double  *vv = is[intf_id].Var_Value, tmp;
  double **Jac;

  /*
   * Zero the calculable fields in the interface structure
   */
  interface_source_zero(&is[intf_id]);

  /*
   * Species number that the boundary condition is being applied on
   */
  kspec =  bc->BC_Data_Int[0];
  wspec_a = kspec;
  wspec_b = kspec + mp_a->Num_Species;

  /*
   * Find the concentration in the a phase
   */
  C_a_mix = 0.0;
  for (k = 0; k < mp_a->Num_Species; k++) {
    C_a_mix += vv[k];
  }

  /*
   * Find the concentration in the b phase, C_b_mix
   */
  if (is[intf_id].SpeciesVT != SPECIES_CONCENTRATION) {
    EH(GOMA_ERROR, "unimplemented complication");
  }
  C_b_mix = 0.0;
  for (k = 0, pos = mp_a->Num_Species; k < mp_b->Num_Species; 
       k++, pos++) {
    C_b_mix += vv[pos];
  }

  /*
   * Find a value of the temperature, even if it is not part of the
   * input.
   */
  if (have_T) {
    index_T = mp_a->Num_Species + mp_b->Num_Species;
    T = vv[index_T];
  } else {
    T = mp_a->reference[TEMPERATURE];
  } 

  /*
   *  Make sure that the interface solution structure understands
   *  that we have carried out the calculation assuming that the
   *  species comprising the unknowns in that structure are
   *  species concentrations, C_k, and that they have units
   *  of mol cm-3.
   */
  is[intf_id].SpeciesVT = SPECIES_CONCENTRATION;

  /*
   * Calculate the Gibbs free energy change for the reaction
   * (cal mol-1 K-1)
   *  HKM -> (ignoring PSChemPot for the time being)
   */
  deltaG = mp_b->SSChemPotData[kspec] - mp_a->SSChemPotData[kspec];

  /*
   * Pseudo reaction rate constant (cm sec-1)
   */
  k_f = bc->BC_Data_Float[0];

  /*
   * Calculate the concentration equilibrium constant for this
   * reaction
   */
  Kc_equil = C_b_mix / C_a_mix  * exp( - deltaG / (RGAS_CALS * T));

  /*
   * Finally, calculate the source term
   */
  tmp = k_f * vv[wspec_b] / Kc_equil;
  is[intf_id].SourceTerm[wspec_b] = k_f * (vv[wspec_a] - vv[wspec_b] / Kc_equil);
  is[intf_id].SourceTerm[wspec_a] = - is[intf_id].SourceTerm[wspec_b];

  if (is[intf_id].Do_Jac) {
    Jac = is[intf_id].JacMatrix;

    Jac[wspec_b][wspec_a] = k_f;
    Jac[wspec_b][wspec_b] = - k_f / Kc_equil;
    for (k = 0; k < mp_a->Num_Species; k++) {
       Jac[wspec_b][k] += - tmp / C_a_mix;
    }
    for (k = mp_a->Num_Species; 
	 k < mp_a->Num_Species + mp_b->Num_Species; k++) {
       Jac[wspec_b][k] +=   tmp / C_b_mix;
    }
    if (have_T) {
      Jac[wspec_b][index_T] = 
	  tmp / Kc_equil * deltaG / (RGAS_CALS * T * T);
    }
    for (k = 0; k < is[intf_id].Num_Terms; k++) {
      Jac[wspec_a][k] = - Jac[wspec_b][k];
    }
  }
  	is[intf_id].Processed[kspec] = TRUE;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

INTERFACE_SOURCE_STRUCT *
is_masstemp_create(MATRL_PROP_STRUCT *mp_1, MATRL_PROP_STRUCT *mp_2,
		   int do_Jac)

    /**********************************************************************
    *
    * is_masstemp_create():
    *
    *  Sets up the storage for an interfacial source structure that
    *  includes the species concentrations for two materials plus
    *  a temperature unknown, tacked onto the end of the structure.
    *
    **********************************************************************/
{
  int num_terms, have_T, k, pos, i;
  INTERFACE_SOURCE_STRUCT *is;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  num_terms = mp_1->Num_Species + mp_2->Num_Species;
  have_T = (upd->vp[pg->imtrx][TEMPERATURE] != -1);
  if (have_T) num_terms++;
  is = interface_source_alloc(num_terms, 0, do_Jac);
  for(i=0 ; i<Num_Interface_Srcs ; i++)    {
  for (k = 0; k < mp_1->Num_Species; k++) {
    vd = find_or_create_vd(MASS_FRACTION, 1, mp_1->MatID, k, pg->imtrx);
    is[i].Var_List[k] = vd;
  }
  for (k = 0, pos = mp_1->Num_Species; k < mp_2->Num_Species;
       k++, pos++) {
    vd = find_or_create_vd(MASS_FRACTION, 1, mp_2->MatID, k, pg->imtrx);
    is[i].Var_List[pos] = vd;
  }
  if (have_T) {
    vd = find_or_create_vd(TEMPERATURE, 1, -1, 0, pg->imtrx);
    is[i].Var_List[num_terms-1] = vd;
  }
  }
  return is;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void 
is_masstemp_fillin(INTERFACE_SOURCE_STRUCT *is,
		   MATRL_PROP_STRUCT *mp_1, MATRL_PROP_STRUCT *mp_2,
		   int speciesVT_desired, const double time, const int intf_id)
    /**********************************************************************
    *
    * is_masstemp_fillin():
    *
    *  Sets up the storage for an interfacial source structure that
    *  includes the species concentrations for two materials plus
    *  a temperature unknown, tacked onto the end of the structure.
    *  It gets the information from the state variable vector within
    *  the material property structure.
    *  
    **********************************************************************/
{
  int num_terms, k;
  double *d_ptr, *is_ptr = is[intf_id].Var_Value;

  /*
   * Fill in the species concentrations for the first problem
   */
  d_ptr = mp_1->StateVector + SPECIES_UNK_0;
  for (k = 0; k < mp_1->Num_Species; k++) {
    *(is_ptr++) = d_ptr[k];
  }
  /*
   * Convert the species vector if necessary
   */
  if (speciesVT_desired !=  mp_1->StateVector_speciesVT) {
    (void) convert_species_var(speciesVT_desired, mp_1,
			       mp_1->StateVector_speciesVT,
			       is[intf_id].Var_Value, time);
  }

  /* 
   * Fill in the species concentrations for the second material
   */
  d_ptr = mp_2->StateVector + SPECIES_UNK_0;
  is_ptr = is[intf_id].Var_Value +  mp_1->Num_Species;
  for (k = 0; k < mp_2->Num_Species; k++) {
    *(is_ptr++) = d_ptr[k];
  }
  /*
   * Convert the species vector if necessary
   */
  if (speciesVT_desired !=  mp_2->StateVector_speciesVT) {
    (void) convert_species_var(speciesVT_desired, mp_2,
			       mp_2->StateVector_speciesVT,
			       is[intf_id].Var_Value + mp_1->Num_Species, time);
  }
  /*
   *  Fill in the temperature if it is part of the problem
   */
  num_terms =  mp_1->Num_Species +  mp_2->Num_Species;
  if (is[intf_id].Num_Terms > num_terms) {
    *is_ptr = mp_1->StateVector[TEMPERATURE];
  }

  /*
   * Store the type of species vector in the interface source 
   * structure
   */
  is[intf_id].SpeciesVT = speciesVT_desired;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

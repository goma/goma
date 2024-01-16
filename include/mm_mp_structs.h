/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 *$Id: mm_mp_structs.h,v 5.20 2010-07-30 21:14:52 prschun Exp $
 */

/* mm_mp_structs.h -- material property structure definitions
 *
 * Revised: Mon Mar 14 15:34:10 MST 1994 pasacki@sandia.gov
 *
 * These are initially just constant floating point values. Later, we
 * may want to use pointers to functions that return complicated values
 * depending on all kinds of other things.
 *
 * Nothing is said about the units -- it is up to the user to input a
 * consistent set of units to solve the transport equations that they
 * are interested in solving. I would recommend MKSA (SI) units throughout,
 * or, even better, have the user set up the constants ahead of time along
 * with an appropriate nondimensionalization.
 *
 * There may be a better way of setting this up, based on the idea that some
 * material properties are associated with bulk materials, some with surfaces or
 * interfaces between two materials or one material and a boundary, and some
 * properties associated with a point (like contact angle at a three phase
 * junction.)
 *
 * The possibility of functional dependence of material properties on not just
 * other variables, but on the spatial gradients of those variables (eg.,
 * viscosity = function of rate of strain tensor) requires some thought to
 * implementation so that mesh derivatives can be performed as easily as
 * possible.
 *
 * Most of these are currently set up with material properties depending
 * only on each variable type - not including multiplicity of concentration
 * which is easily added by adding a new dimension [MAX_CONC]
 */

#ifndef GOMA_MM_MP_STRUCTS_H
#define GOMA_MM_MP_STRUCTS_H

#include "mm_as_structs.h"
#include "mm_elem_block_structs.h"
#include "rf_fem_const.h"

/*
 * New size variables provide hints to interprocessor communications as to
 * how big each of the "u_...." variables really is.
 */

#ifndef GOMA_CK_NAME_DEF
#define GOMA_CK_NAME_DEF
typedef char CK_NAME[64]; /* Typedefs for common names used for naming
                             domains and species */
typedef char CK_NAME_STR[64];
#endif

struct Material_Properties {
  int MatID;                        /* Material ID number for this material. This is
                                     * >= 0, and consistent across all processors,
                                     * unique, and equal to the index into the global
                                     * array of Material property structures */
  char Material_Name[MAX_MATLNAME]; /* Character string name for the material */
  int Num_Matrl_Elem_Blk;           /* Number of element blocks comprising this
                                     * material                                     */

  int *Matrl_Elem_Blk_Ids; /* Malloced list of element block ids comprising
                            * this material
                            * Length = Num_Matrl_Elem_Blk */

  int DefaultDatabase; /* Default place to look for physical property data
                        *      GOMA_MAT = 0 = Default
                        *      CHEMKIN_MAT = 1: Look up all property data
                        *                       in the chemkin data bases
                        */
  int Num_Species;     /* Number of species defined for this material */
  int Num_Species_Eqn; /* Number of species equations solved for in this
                        * material: NOTE: this is usually one less than
                        * the total number of species in the material due
                        * to the implicit imposition of the sum MF = 1
                        * constraint                                    */
  int Dropped_Last_Species_Eqn;
  /* True if the last species equation and only
   * the last species equn has been dropped from
   * equation system */
  int NonDiluteFormulation; /* If this flag is set, then there are additional
                             * species, whose concentrations are determined
                             * by the sum MF = 1 constraint, and/or by the
                             * equation of state, or by a volumetric
                             * Dirichlet condition.
                             */
  char **Species_Names;     /* Pointer to a vector of species names for the
                             * current material */
  int PhaseID[MAX_CONC];    /* Phase Identification for all of the species in
                             * the current material domain. If there is just
                             * one phase in the domain, then this vector will
                             * be identically zero. */
  int Species_Var_Type;     /* Species Variable type for the species defined
                             * in this material.
                             * Possible values are SPECIES_MOLE_FRACTION,
                             * SPECIES_MASS_FRACTION, SPECIES_CONCENTRATION,
                             * SPECIES_CAP_PRESSURE, etc.
                             * The acceptable values are listed
                             * in rf_fem_const.h. This variable influences all
                             * aspects of the species conservation equation,
                             * as well as the names that are put into
                             * the output files.
                             * The value of this entry will also determine
                             * the prefix associated with the name in the
                             * exodus file. */
  double StateVector[MAX_VARIABLE_TYPES];
  /* This is the state vector for calculation of
   * physical properties of the material.
   * Note: a state vector usually includes
   *   specification of the Temperature, Pressure,
   *   and species concentrations, and anything
   *   else that would nail down the current
   *   thermodynamic state of the system at a
   *   point.
   * -> Note, the species unknowns are put into
   *    the SPECIES_UNK_0 slots in this vector.
   */
  int StateVector_speciesVT; /* Current Species variable types for the
                              * species unknown in the state vector. Note
                              * these may be different than the Species_Var_Type
                              * field also in this structure.
                              */
  int Volumetric_Dirichlet_Cond[MAX_CONC];
  /* Flag for each species in the material to
   * indicate whether the degree of freedom is really
   * an algebraic constraint dictated by the sum MF=1
   * constraint or by the equation of state,
   * (a value of  TRUE), or whether it is a normal
   * degree of freedom  whose equation involves
   * advection, diffusion,  and reaction, etc
   * (FALSE).
   */
  int Num_Porous_Eqn; /* Number of porous media eqns solved for in this
                       * material
                       */

  int Num_Porous_Shell_Eqn; /* Number of porous media shell eqns solved for in this
                             * material
                             */

  int Porous_Eqn[MAX_PMV];             /* array containing the porous media equations
                                        * active in this material
                                        */
  int Porous_Shell_Eqn[MAX_POR_SHELL]; /* array containing the porous shell equations
                                        * active in this material
                                        */

  char **Porous_Names; /* Pointer to a vector of porous phase names for the
                        * current material */

  dbl thermal_conductivity; /* Yeah, you could make this a tensor... */
  dbl d_thermal_conductivity[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_thermal_conductivity;
  dbl *u_thermal_conductivity;
  int ConductivityModel;
  int thermal_conductivity_tableid;
  int Ewt_funcModel;
  int Energy_Div_Term;
  dbl Ewt_func;
  int Rst_funcModel;
  dbl Rst_func;
  dbl Rst_diffusion;
  dbl Rst_func_supg;
  int thermal_cond_external_field;

  dbl electrical_conductivity; /* Yeah, you could make this a tensor... */
  dbl d_electrical_conductivity[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_electrical_conductivity;
  dbl *u_electrical_conductivity;
  int Elec_ConductivityModel;
  int elec_cond_external_field;

  dbl permittivity;
  dbl permittivity_imag;
  dbl d_permittivity[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_permittivity;
  dbl *u_permittivity;
  int PermittivityModel;

  dbl elect_surf_diffusivity;
  dbl d_elect_surf_diffusivity[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_elect_surf_diffusivity;
  dbl *u_elect_surf_diffusivity;
  int Elect_Surf_DiffusivityModel;

  dbl magnetic_permeability;
  dbl d_magnetic_permeability[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_magnetic_permeability;
  dbl *u_magnetic_permeability;
  int MagneticPermeabilityModel;

  dbl shell_user_par;
  dbl d_shell_user_par[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_shell_user_par;
  dbl *u_shell_user_par;
  int Shell_User_ParModel;

  dbl acoustic_impedance;
  dbl d_acoustic_impedance[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_acoustic_impedance;
  dbl *u_acoustic_impedance;
  int Acoustic_ImpedanceModel;
  int acoustic_impedance_tableid;

  dbl wave_number; /* Acoustic Wave Number */
  dbl d_wave_number[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_wave_number;
  dbl *u_wave_number;
  int wave_numberModel; /* CONSTANT */
  int wave_number_tableid;

  dbl acoustic_ksquared_sign; /* Sign of wavenumber squared -- captures imaginary wavenumbers */
  int Ksquared_SignModel;     /* CONSTANT */

  dbl acoustic_absorption;
  dbl d_acoustic_absorption[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_acoustic_absorption;
  dbl *u_acoustic_absorption;
  int Acoustic_AbsorptionModel;
  int acoustic_absorption_tableid;

  dbl refractive_index;
  dbl d_refractive_index[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_refractive_index;
  dbl *u_refractive_index;
  int Refractive_IndexModel;
  int refractive_index_tableid;

  dbl light_absorption;
  dbl d_light_absorption[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_light_absorption;
  dbl *u_light_absorption;
  int Light_AbsorptionModel;
  int light_absorption_tableid;

  dbl extinction_index;
  dbl d_extinction_index[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_extinction_index;
  dbl *u_extinction_index;
  int Extinction_IndexModel;
  int extinction_index_tableid;

  int VoltageFormulation; /* Used to select k for potential equation
                           * Permittivity or conductivity (default) */

  //! Viscosity of the phase, usually evaluated at the gauss point
  dbl viscosity;
  int len_u_viscosity;
  dbl *u_viscosity;
  //! Derivative of the viscosity wrt the independent unknowns in the problem
  dbl d_viscosity[MAX_VARIABLE_TYPES + MAX_CONC];

  //! Second derivatives of the viscosity wrt the independent unknowns in the
  //! problem
  dbl d2_viscosity[MAX_VARIABLE_TYPES + MAX_CONC + 1];

  int ViscosityModel;
  int viscosity_tableid;

  dbl dilationalViscosity;
  int len_u_dilationalViscosity;
  dbl *u_dilationalViscosity;
  dbl dilationalViscosityRatio;
  dbl d_dilationalViscosityRatio[MAX_VARIABLE_TYPES + MAX_CONC];
  dbl dilationalViscosityMultiplier;
  //! Volume fraction of Gas for multiphase material models
  //!  (Currently there is one model, foam)
  dbl volumeFractionGas;
  //! Derivative of the volume fraction of gas wrt the independent
  //! unknowns in the problem
  dbl d_volumeFractionGas[MAX_VARIABLE_TYPES + MAX_CONC];

  //! Integer representing the model for the dilationa viscosity
  int DilationalViscosityModel;
  int dilationalViscosity_tableid;

  /* Integer for wall distance external field index */
  int dist_wall_ext_field_index;

  int Mwt_funcModel;
  dbl Mwt_func;

  int SAwt_funcModel;
  dbl SAwt_func;

  dbl surface_tension;
  int len_u_surface_tension;
  dbl *u_surface_tension;
  dbl d_surface_tension[MAX_VARIABLE_TYPES + MAX_CONC];
  int SurfaceTensionModel;

  dbl SurfaceDiffusionCoeffProjectionEqn;

  dbl heat_capacity;
  int len_u_heat_capacity;
  dbl *u_heat_capacity;
  dbl d_heat_capacity[MAX_VARIABLE_TYPES + MAX_CONC];
  int HeatCapacityModel;
  int heat_capacity_tableid;

  dbl Volume_Expansion;
  int len_u_Volume_Expansion;
  dbl *u_Volume_Expansion;
  dbl d_Volume_Expansion[MAX_VARIABLE_TYPES + MAX_CONC];
  int VolumeExpansionModel;

  dbl density;                    /* Value of the density given the current
                                   * value of the state variable, StateVariable,
                                   * also in this structure. */
  PROPERTYJAC_STRUCT *DensityJac; /* When necessary, this includes the derivatives
                                   * of the density wrt the state variables */
  int SBM_Length_enabled;

  int DensityModel;  /* Model type: for types, see mm_mp_const.h  */
  int len_u_density; /* Constants for user-defined density model*/
  dbl *u_density;    /* Constants for user-defined density model*/
  /*
   * HKM Note: d_density is now superfluous and will be eliminated eventually
   */
  dbl d_density[MAX_VARIABLE_TYPES + MAX_CONC]; /* Note this includes the thermal */
                                                /* expansion coefficient, beta */
                                                /* = d(density)/d(temperature) */
  /*
   * HKM Note: d_rho_dT and d_rho_dC contain duplicate information to
   *           DensityJac and will be either eliminated or changed to
   *           a form without the MDE dimension.
   */
  dbl d_rho_dT[MDE];
  dbl d_rho_dC[MAX_CONC][MDE];

  int PBE_BA_Type;

  int Spwt_funcModel;
  dbl Spwt_func;

  int SpSSPG_funcModel;
  dbl SpSSPG_func;

  int SpYZbeta_funcModel;
  dbl SpYZbeta_func;
  dbl SpYZbeta_value;

  int Momentwt_funcModel;
  dbl Momentwt_func;

  int MomentSSPG_funcModel;
  dbl MomentSSPG_func;

  int MomentShock_funcModel;
  dbl MomentShock_func;
  dbl MomentShock_Ref[MAX_MOMENTS];

  int MomentDiffusivityModel;
  dbl MomentDiffusivity;

  int MomentSecondLevelSetDiffusivityModel;
  dbl MomentSecondLevelSetDiffusivity;

  int MomentLevelSetDiffusionOnly;

  dbl diffusivity[MAX_CONC];
  int DiffusivityModel[MAX_CONC];
  int diffusivity_tableid[MAX_CONC];
  double SpeciesSecondLevelSetDiffusivity[MAX_CONC];
  int SpeciesOnlyDiffusion[MAX_CONC];   /* STANDARD or TAYLOR_GALERKIN */
  int SpeciesTimeIntegration[MAX_CONC]; /* STANDARD or TAYLOR_GALERKIN */
  int FreeVolSolvent[MAX_CONC];         /* Solvent identity array for multi-C FV*/
  int len_u_diffusivity[MAX_CONC];
  dbl *u_diffusivity[MAX_CONC];
  dbl d_diffusivity[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  dbl diffusivity_gen_fick[MAX_CONC][MAX_CONC]; /* generalized fickian diffusion ACS 4/00 */
  dbl loadfv[MAX_CONC][15];
  dbl d_diffusivity_gf[MAX_CONC][MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  /* Property parameters needed for the Stefan-Maxwell flux model; KSC: 7/98, 2/99, 2/02;
     Molecular weight is used in VL_POLY; ACS 8/98 */
  dbl diffusivity_Stefan_Maxwell[MAX_CONC][MAX_CONC];
  dbl u_diffusivity_Stefan_Maxwell[MAX_CONC][MAX_CONC][3]; /* KSC, 9/04 */
  dbl reaction_rate;
  int ReactionRateModel;
  int len_u_reaction_rate; /* Constants for user-defined reaction-rate model*/
  dbl *u_reaction_rate;    /* Constants for user-defined reaction-rate model*/
  dbl d_reaction_rate[MAX_VARIABLE_TYPES + MAX_CONC];
  dbl molecular_weight[MAX_CONC];
  int MolecularWeightModel[MAX_CONC];
  dbl d_molecular_weight[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
  dbl charge_number[MAX_CONC];
  int ChargeNumberModel[MAX_CONC];
  dbl d_charge_number[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
  dbl solution_temperature;
  int SolutionTemperatureModel;
  int len_u_solution_temperature; /* Constants for user-defined solution temperature model*/
  dbl *u_solution_temperature;    /* Constants for user-defined solution temperature model*/
  dbl d_solution_temperature[MAX_VARIABLE_TYPES + MAX_CONC];
  dbl electrolyte_temperature;
  dbl electrolyte_conductivity;
  dbl electrode_potential;
  dbl thermodynamic_potential;
  int ThermodynamicPotentialModel;
  int len_u_thermodynamic_potential;
  dbl *u_thermodynamic_potential;
  dbl d_thermodynamic_potential[MAX_VARIABLE_TYPES + MAX_CONC];
  dbl interfacial_area;
  int InterfacialAreaModel;
  int len_u_interfacial_area;
  dbl *u_interfacial_area;
  dbl d_interfacial_area[MAX_VARIABLE_TYPES + MAX_CONC];

  /* diffusivitives for each term  of HYDRO Model */

  int GamDiffType[MAX_CONC];
  int MuDiffType[MAX_CONC];
  int GravDiffType[MAX_CONC];
  int SBM_Type[MAX_CONC];
  int FickDiffType[MAX_CONC];
  int CurvDiffType[MAX_CONC];
  int QTensorDiffType[MAX_CONC];
  int NSCoeffType[MAX_CONC];
  int len_u_gadiffusivity[MAX_CONC]; /*this is currently defined for MPI*/
  int len_u_mdiffusivity[MAX_CONC];  /*this is currently defined for MPI*/
  int len_u_fdiffusivity[MAX_CONC];
  int len_u_gdiffusivity[MAX_CONC]; /*this is currently defined for MPI*/
  int len_SBM_Lengths2[MAX_CONC];   /*this is currently defined for MPI*/
  int len_u_cdiffusivity[MAX_CONC]; /*this is currently defined for MPI*/
  int len_u_qdiffusivity[MAX_CONC]; /*this is currently defined for MPI*/
  int len_u_nscoeff[MAX_CONC];

  dbl gam_diffusivity[MAX_CONC]; /* Kc from shear-gradient term */
  dbl *u_gadiffusivity[MAX_CONC];
  dbl mu_diffusivity[MAX_CONC]; /* Kmu from viscosity-gradient term */
  dbl *u_mdiffusivity[MAX_CONC];
  dbl f_diffusivity[MAX_CONC]; /* normal Fickian diffusion term */
  dbl *u_fdiffusivity[MAX_CONC];
  dbl g_diffusivity[MAX_CONC]; /* hindered settling function */
  dbl SBM_Lengths[MAX_CONC];   /* hindered settling function */
  dbl NSCoeff[MAX_CONC];
  dbl *u_nscoeff[MAX_CONC];
  dbl *u_gdiffusivity[MAX_CONC]; /*this is currently defined for MPI*/
  dbl *SBM_Lengths2[MAX_CONC];   /*this is currently defined for MPI*/
  dbl cur_diffusivity[MAX_CONC]; /* curvature induced migration term */
  dbl *u_cdiffusivity[MAX_CONC];
  dbl q_diffusivity[MAX_CONC][DIM]; /* Q tensor diffusion components. */
  dbl *u_qdiffusivity[MAX_CONC];

  /* Parameters for Ryan's Qtensor model */
  int QtensorExtensionPModel;
  int QtensorNctModel;
  dbl Qtensor_Extension_P;
  dbl Qtensor_Nct;

  /* Reference Concentration Model is defined by itself */

  int RefConcnModel[MAX_CONC];
  dbl reference_concn[MAX_CONC];
  dbl *u_reference_concn[MAX_CONC];
  int len_u_reference_concn[MAX_CONC];

  int AdvectiveScalingModel[MAX_CONC];
  dbl AdvectiveScaling[MAX_CONC];
  /*
   *   ExtrinsicIndependentSpeciesVar[i] is true if
   *   the species variable i is extrinsic. This means
   *   that there is an extra c (Del dot V) term in the
   *   advection operator when there is dilational flow.
   */
  int ExtrinsicIndependentSpeciesVar[MAX_CONC];

  /*  molar volume */
  /*
   *  HKM -> The MolarVolumeModel and SpecificVolumeModel are
   *         in direct competition with the DensityModel field.
   *         They are the same thing. And, should not have a
   *         [MAX_CONC] vector associated with the quantity.
   *         Moreover, these two fields are not used within the
   *         code, to my knowledge. Thus, they may disappear in
   *         the future.
   */
  int MolarVolumeModel[MAX_CONC];
  int SpecificVolumeModel[MAX_CONC];
  dbl molar_volume[MAX_CONC];
  dbl specific_volume[MAX_CONC];
  dbl flory_param[MAX_CONC][MAX_CONC];

  /*
   * HKM - Specification of the chemical potential
   *       -> Really, better to do it via calls to chemkin.
   *          But, I'll include it here for BC specification
   *          and code development
   *     ( mu[i] = RTln(act_coeff[i]) + mu_start_[T, P] )
   *     ( mu_star_i = RTln(P/1atm) + mu_ss_i[T])
   *        -> See Denbigh for equations
   */
  int SSChemPotModel[MAX_CONC];
  dbl SSChemPotData[MAX_CONC];
  int PSChemPotModel[MAX_CONC];
  dbl PSChemPotData[MAX_CONC];
  int ChemPotModel[MAX_CONC];
  dbl ChemPotData[MAX_CONC];

  /* some properties for porous materials */

  int PorousMediaType;
  int CapStress;

  int i_ys; /* species number of species eqn */

  dbl porosity;
  int PorosityModel;
  int len_u_porosity; /* Constants for user-defined porosity model*/
  dbl *u_porosity;    /* Constants for user-defined porosity model*/
  dbl d_porosity[MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV];
  int porosity_external_field_index;

  dbl porous_compressibility;     /* Porous (rock) compressibility */
  dbl initial_porosity;           /* Initializes the external field */
  int PorousCompressibilityModel; /* when this is set to "CONST_INIT" */
  int len_u_porous_compressibility;
  dbl *u_porous_compressibility;

  dbl matrix_density; /* density of the solid material of a porous medium       */
  int PorousMatrixDensityModel;
  dbl d_matrix_density[MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV];

  dbl specific_heat; /* specific heat of the solid material of a porous medium */
  int PorousSpecificHeatModel;
  dbl d_specific_heat[MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV];

  dbl permeability;
  dbl permeability_imag;
  dbl perm_tensor[DIM][DIM];
  int PermeabilityModel;
  dbl d_permeability[MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV];
  int len_u_permeability;
  dbl *u_permeability;
  dbl d_perm_tensor_dx[DIM][DIM][DIM][MDE];
  dbl d_perm_tensor[DIM][DIM][MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV];
  int perm_external_field_index;

  dbl PorousLiqCompress;
  int PorousLiquidCompressModel;
  dbl d_PorousLiquidCompres[MAX_VARIABLE_TYPES + MAX_CONC];

  dbl PorousLiqRefPress;
  int PorousLiqRefPressModel;
  dbl d_PorousLiqRefPress[MAX_VARIABLE_TYPES + MAX_CONC];

  dbl rel_gas_perm;
  int RelGasPermModel;
  int len_u_rel_gas_perm;
  dbl *u_rel_gas_perm;
  dbl d_rel_gas_perm[MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV];

  dbl rel_liq_perm;
  int RelLiqPermModel;
  int len_u_rel_liq_perm;
  dbl *u_rel_liq_perm;
  dbl d_rel_liq_perm[MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV];
  int rel_liq_perm_external_field_index;

  dbl saturation;
  int SaturationModel;
  int len_u_saturation;
  dbl *u_saturation;
  dbl d_saturation[MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV];
  dbl d_d_saturation[MAX_VARIABLE_TYPES + MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
  int saturation_tableid;
  int SAT_external_field_index;

  dbl cap_pres;
  int CapPresModel;
  int len_u_cap_pres;
  dbl *u_cap_pres;
  dbl d_cap_pres[MAX_VARIABLE_TYPES + MAX_CONC];
  dbl d_d_cap_pres[MAX_VARIABLE_TYPES + MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
  int cap_pres_tableid;
  int cap_pres_external_field_index;

  /*
   *  Porous_wt_funcModel: This is where you specify standard weighting
   *                       or SUPG. The wt_func should normally always be
   *                       equal to 1, because 1 corresponds to the theoretical
   *                       amount of compensating artificial diffusivity
   *                       that makes up for the lack of actual diffusivity
   *                       in Galerkin's weighting as you get to high
   *                       peclet numbers.
   */
  int Porous_wt_funcModel;
  dbl Porous_wt_func;

  /*
   *  Porous_Mass_Lump:  If true mass lumping of the porous media
   *                     inventory terms is carried out.
   */
  int Porous_Mass_Lump;

  dbl porous_diffusivity[MAX_PMV];
  int PorousDiffusivityModel[MAX_PMV];
  int PorousTimeIntegration[MAX_PMV]; /* STANDARD or TAYLOR_GALERKIN */

  int len_u_porous_diffusivity[MAX_PMV];
  dbl *u_porous_diffusivity[MAX_PMV];
  dbl d_porous_diffusivity[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC];

  dbl porous_latent_heat_vap[MAX_PMV];
  int PorousLatentHeatVapModel[MAX_PMV];
  dbl d_porous_latent_heat_vap[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC];

  dbl porous_latent_heat_fusion[MAX_PMV];
  int PorousLatentHeatFusionModel[MAX_PMV];
  dbl d_porous_latent_heat_fusion[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC];

  dbl porous_vapor_pressure[MAX_PMV];
  dbl d_porous_vapor_pressure[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_porous_vapor_pressure[MAX_PMV];
  dbl *u_porous_vapor_pressure[MAX_PMV];
  int PorousVaporPressureModel[MAX_PMV];

  dbl porous_gas_constants;
  dbl d_porous_gas_constants[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_porous_gas_constants;
  dbl *u_porous_gas_constants;
  int PorousGasConstantsModel;

  dbl porous_sink_constants;
  dbl d_porous_sink_constants[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_porous_sink_constants;
  dbl *u_porous_sink_constants;
  int PorousSinkConstantsModel;

  dbl porous_vol_expansion[MAX_PMV];
  dbl d_porous_vol_expansion[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC];
  int PorVolExpModel[MAX_PMV];
  int len_u_porous_vol_expansion[MAX_PMV];
  dbl *u_porous_vol_expansion[MAX_PMV];

  dbl porous_molecular_weight[MAX_PMV];
  int PorousMolecularWeightModel[MAX_PMV];
  dbl d_porous_molecular_weight[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC];

  int moment_growth_model;
  dbl moment_growth_scale;
  dbl moment_growth_reference_pressure;

  int moment_coalescence_model;
  dbl moment_coalescence_scale;

  /*
   * Source terms...
   */

  dbl momentum_source[DIM];
  dbl d_momentum_source[DIM][MAX_VARIABLE_TYPES + MAX_CONC];
  int MomentumSourceModel;
  int len_u_momentum_source;
  dbl *u_momentum_source;

  dbl heat_source;
  dbl d_heat_source[MAX_VARIABLE_TYPES + MAX_CONC];
  int HeatSourceModel;
  int len_u_heat_source;
  dbl *u_heat_source;

  dbl species_source[MAX_CONC];
  dbl d_species_source[MAX_VARIABLE_TYPES + MAX_CONC];
  int SpeciesSourceModel[MAX_CONC];
  int len_u_species_source[MAX_CONC];
  dbl *u_species_source[MAX_CONC];
  dbl Jac_Species_Source[MAX_CONC * MAX_CONC];
  int species_source_external_field_index;

  dbl species_vol_expansion[MAX_CONC];
  dbl d_species_vol_expansion[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
  int SpecVolExpModel[MAX_CONC];
  int len_u_species_vol_expansion[MAX_CONC];
  dbl *u_species_vol_expansion[MAX_CONC];

  dbl mass_source;
  dbl d_mass_source[MAX_VARIABLE_TYPES + MAX_CONC];
  int MassSourceModel;
  int len_u_mass_source;
  dbl *u_mass_source;

  dbl mesh_source[DIM];
  dbl real_solid_source[DIM];
  dbl d_mesh_source[DIM][MAX_VARIABLE_TYPES + MAX_CONC];
  int MeshSourceModel;
  int len_u_mesh_source;
  dbl *u_mesh_source;
  int RealSolidSourceModel;

  dbl current_source;
  dbl d_current_source[MAX_VARIABLE_TYPES + MAX_CONC];
  int CurrentSourceModel;
  int len_u_current_source;
  dbl *u_current_source;

  dbl moment_source;
  dbl d_moment_source[MAX_VARIABLE_TYPES + MAX_CONC];
  int MomentSourceModel;
  int len_u_moment_source;
  dbl *u_moment_source;

  dbl heightU;
  dbl d_heightU[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_heightU_function_constants;
  dbl *u_heightU_function_constants;
  int HeightUFunctionModel;
  int heightU_ext_field_index;
  int heightU_function_constants_tableid;

  dbl heightL;
  dbl d_heightL[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_heightL_function_constants;
  dbl *u_heightL_function_constants;
  int HeightLFunctionModel;
  int heightL_ext_field_index;
  int heightL_function_constants_tableid;

  dbl veloU[DIM];
  dbl d_veloU[DIM][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_veloU_function_constants;
  dbl *u_veloU_function_constants;
  int VeloUFunctionModel;

  dbl veloL[DIM];
  dbl d_veloL[DIM][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_veloL_function_constants;
  dbl *u_veloL_function_constants;
  int VeloLFunctionModel;

  dbl dcaU;
  dbl d_dcaU[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_dcaU_function_constants;
  dbl *u_dcaU_function_constants;
  int DcaUFunctionModel;

  dbl dcaL;
  dbl d_dcaL[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_dcaL_function_constants;
  dbl *u_dcaL_function_constants;
  int DcaLFunctionModel;

  int FSIModel;
  int LubIntegrationModel;
  int LubInt_NGP;
  dbl Lub_gpts[MAX_LUB_NGP];
  dbl Lub_wts[MAX_LUB_NGP];
  dbl LubInt_PL;

  int Lub_Curv_DiffModel;
  double Lub_Curv_Diff;

  int TurbulentLubricationModel;

  dbl lubsource;
  dbl d_lubsource[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_lubsource;
  dbl *u_lubsource_function_constants;
  int LubSourceModel;

  dbl lubmomsource[DIM];
  dbl d_lubmomsource[DIM][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_lubmomsource;
  int LubMomSourceModel;

  dbl FilmEvap;
  dbl d_FilmEvap[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_FilmEvap_function_constants;
  dbl *u_FilmEvap_function_constants;
  int FilmEvapModel;

  dbl DisjPress;
  dbl d_DisjPress[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_DisjPress_function_constants;
  dbl *u_DisjPress_function_constants;
  int DisjPressModel;

  dbl SlipCoeff;
  dbl d_SlipCoeff[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_SlipCoeff_function_constants;
  dbl *u_SlipCoeff_function_constants;
  int SlipCoeffModel;

  dbl DiffCoeff;
  dbl d_DiffCoeff[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_DiffCoeff_function_constants;
  dbl *u_DiffCoeff_function_constants;
  int DiffCoeffModel;

  dbl PorousShellClosedPorosity;
  dbl d_PorousShellClosedPorosity[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellClosedPorosity_function_constants;
  dbl *u_PorousShellClosedPorosity_function_constants;
  int PorousShellClosedPorosityModel;
  int por_shell_closed_porosity_ext_field_index;

  dbl PorousShellClosedHeight;
  dbl d_PorousShellClosedHeight[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellClosedHeight_function_constants;
  dbl *u_PorousShellClosedHeight_function_constants;
  int PorousShellClosedHeightModel;
  int por_shell_closed_height_ext_field_index;

  dbl PorousShellClosedRadius;
  dbl d_PorousShellClosedRadius[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellClosedRadius_function_constants;
  dbl *u_PorousShellClosedRadius_function_constants;
  int PorousShellClosedRadiusModel;
  int por_shell_closed_radius_ext_field_index;

  dbl PorousShellClosedP0;
  dbl d_PorousShellClosedP0[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellClosedP0_function_constants;
  dbl *u_PorousShellClosedP0_function_constants;
  int PorousShellClosedP0Model;

  dbl PorousShellPatm;
  dbl d_PorousShellPatm[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellPatm_function_constants;
  dbl *u_PorousShellPatm_function_constants;
  int PorousShellPatmModel;

  dbl PorousShellPref;
  dbl d_PorousShellPref[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellPref_function_constants;
  dbl *u_PorousShellPref_function_constants;
  int PorousShellPrefModel;

  dbl PorousShellCrossKappa;
  dbl d_PorousShellCrossKappa[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellCrossKappa_function_constants;
  dbl *u_PorousShellCrossKappa_function_constants;
  int PorousShellCrossKappaModel;
  int Xperm_external_field_index;

  dbl PorousShellInitPorePres;
  dbl d_PorousShellInitPorePres[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellInitPorePres_function_constants;
  dbl *u_PorousShellInitPorePres_function_constants;
  int PorousShellInitPorePresModel;

  dbl PorousShellDiffusivity;
  dbl d_PorousShellDiffusivity[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellDiffusivity_function_constants;
  dbl *u_PorousShellDiffusivity_function_constants;
  int PorousShellDiffusivityModel;

  dbl PorousShellRT;
  dbl d_PorousShellRT[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellRT_function_constants;
  dbl *u_PorousShellRT_function_constants;
  int PorousShellRTModel;

  dbl PorousShellHenry;
  dbl d_PorousShellHenry[MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellHenry_function_constants;
  dbl *u_PorousShellHenry_function_constants;
  int PorousShellHenryModel;

  dbl PorousShellPorosity[MAX_POR_SHELL];
  dbl d_PorousShellPorosity[MAX_POR_SHELL][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellPorosity[MAX_POR_SHELL];
  dbl *u_PorousShellPorosity[MAX_POR_SHELL];
  int PorousShellPorosityModel[MAX_POR_SHELL];
  int por_shell_porosity_ext_field_index[MAX_POR_SHELL];

  dbl PorousShellHeight[MAX_POR_SHELL];
  dbl d_PorousShellHeight[MAX_POR_SHELL][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellHeight[MAX_POR_SHELL];
  dbl *u_PorousShellHeight[MAX_POR_SHELL];
  int PorousShellHeightModel[MAX_POR_SHELL];
  int por_shell_height_ext_field_index[MAX_POR_SHELL];

  dbl PorousShellPermeability[MAX_POR_SHELL];
  dbl PorousShellPermTensor[MAX_POR_SHELL][DIM][DIM];
  dbl d_PorousShellPermeability[MAX_POR_SHELL][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellPermeability[MAX_POR_SHELL];
  dbl *u_PorousShellPermeability[MAX_POR_SHELL];
  int PorousShellPermeabilityModel[MAX_POR_SHELL];
  int por_shell_permeability_ext_field_index[MAX_POR_SHELL];

  dbl PorousShellCrossPermeability[MAX_POR_SHELL];
  dbl d_PorousShellCrossPermeability[MAX_POR_SHELL][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellCrossPermeability[MAX_POR_SHELL];
  dbl *u_PorousShellCrossPermeability[MAX_POR_SHELL];
  int PorousShellCrossPermeabilityModel[MAX_POR_SHELL];
  int por_shell_cross_permeability_ext_field_index[MAX_POR_SHELL];

  dbl PorousShellRelPerm[MAX_POR_SHELL];
  dbl d_PorousShellRelPerm[MAX_POR_SHELL][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellRelPerm[MAX_POR_SHELL];
  dbl *u_PorousShellRelPerm[MAX_POR_SHELL];
  int PorousShellRelPermModel[MAX_POR_SHELL];
  int por_shell_rel_perm_ext_field_index[MAX_POR_SHELL];

  dbl PorousShellCapPres[MAX_POR_SHELL];
  dbl d_PorousShellCapPres[MAX_POR_SHELL][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_PorousShellCapPres[MAX_POR_SHELL];
  dbl *u_PorousShellCapPres[MAX_POR_SHELL];
  int PorousShellCapPresModel[MAX_POR_SHELL];
  int por_shell_cap_pres_ext_field_index[MAX_POR_SHELL];
  int por_shell_cap_pres_hyst_num_switch_ext_field_index[MAX_POR_SHELL];
  int por_shell_cap_pres_hyst_curve_type_ext_field_index[MAX_POR_SHELL];

  // EM incident wave, mix between boundary condition / source / problem property
  int IncidentWaveModel;
  dbl incident_wave;
  int len_u_incident_wave;
  dbl *u_incident_wave;

  /*
   * Boundary conditions...(these quantities and the geometric surface
   * parameters may better be place in their own side_set/node_set sort
   * of pointer...
   */

  dbl reference[MAX_VARIABLE_TYPES]; /* Mostly for reference temperature for */
                                     /* Boussinesq term, but conceivably also */
                                     /* extended to concentration, and other */
                                     /* variables, too. */
  int ReferenceModel[MAX_VARIABLE_TYPES];

  dbl melting_point_liquidus;
  dbl d_melting_point_liquidus[MAX_VARIABLE_TYPES + MAX_CONC];
  int LiquidusModel;

  dbl melting_point_solidus;
  dbl d_melting_point_solidus[MAX_VARIABLE_TYPES + MAX_CONC];
  int SolidusModel;

  dbl latent_heat_vap[MAX_CONC];
  int LatentHeatVapModel[MAX_CONC];
  dbl d_latent_heat_vap[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  dbl latent_heat_fusion[MAX_CONC];
  int LatentHeatFusionModel[MAX_CONC];
  dbl d_latent_heat_fusion[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  dbl vapor_pressure[MAX_CONC];
  dbl d_vapor_pressure[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
  int len_u_vapor_pressure[MAX_CONC];
  dbl *u_vapor_pressure[MAX_CONC];
  int VaporPressureModel[MAX_CONC];

  dbl heat_transfer_coefficient;
  dbl d_heat_transfer_coefficient[MAX_VARIABLE_TYPES + MAX_CONC];

  dbl mass_transfer_coefficient[MAX_CONC];
  dbl d_mass_transfer_coefficient[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  dbl mass_flux[MAX_CONC];
  dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  dbl species_activity[MAX_CONC];
  dbl d_species_activity[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  dbl geometry_parameters[MAX_GEOMETRY_PARMS]; /* defined in rf_fem_const.h */

  //! FlowingLiquid Viscosity -> Used in the Brinkman Equation
  //  -> also reused in some situations for the FOAM_EPOXY model to hold
  //     the pure liquid phase viscosity and its derivatives.
  dbl FlowingLiquid_viscosity;
  dbl d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES + MAX_CONC];
  int FlowingLiquidViscosityModel;
  dbl *u_FlowingLiquid_viscosity;
  int len_u_FlowingLiquid_viscosity;

  //! Inertia_coefficient -> used in the Brinkman Equation
  dbl Inertia_coefficient;
  dbl d_Inertia_coefficient[MAX_VARIABLE_TYPES + MAX_CONC];
  int InertiaCoefficientModel;

  // TFMP structure for material properties function constants
  int tfmp_density_model;
  int len_tfmp_density_const;
  dbl *tfmp_density_const;

  int tfmp_viscosity_model;
  int len_tfmp_viscosity_const;
  dbl *tfmp_viscosity_const;

  // TFMP indicators for the stabilizing diffusivity corrector
  int tfmp_diff_model;
  int len_tfmp_diff_const;
  dbl *tfmp_diff_const;

  // TFMP variables for wt function application
  int tfmp_wt_model;
  int tfmp_wt_len;
  dbl tfmp_wt_const;

  int tfmp_mass_lump;
  int tfmp_clipping;
  dbl tfmp_clip_strength;

  int tfmp_rel_perm_model;
  int len_tfmp_rel_perm_const;
  dbl *tfmp_rel_perm_const;

  int tfmp_dissolution_model;
  int len_tfmp_dissolution_const;
  dbl *tfmp_dissolution_const;

  int tfmp_drop_lattice_model;
  int len_tfmp_drop_lattice_const;
  dbl *tfmp_drop_lattice_const;

  /* Properties for Fixed Deformable Roller
  CTPM "coupled two-phase membrane"
  */

  int shell_tangent_model;
  int len_shell_tangent_seed_vec_const;
  dbl *shell_tangent_seed_vec_const;

  int shell_moment_tensor_model;

  // Elastohydrodynamic Lubrication
  int ehl_gap_model;
  int ehl_normal_method;
  int ehl_integration_kind;

  int table_index;

  struct Data_Table *table;

  struct Second_LS_Phase_Properties *mp2nd;
};

//! Structure containing all of the physical constants related to viscosity
/*!
 *   Generalized Newtonian Models:
 *   Newtonian, Power Law, Carreau, Bingham,  Carreau_wlf
 *   or Carreau_Suspension etc
 *
 *  This structure is allocated for every material. For each material
 *  there are MAX_MODES numbers of these structures.
 */
struct Generalized_Newtonian {
  //! Integer describing the viscosity model
  /*!
   *   Generalized Newtonian Models:
   *   Newtonian, Power Law, Carreau, Bingham,  Carreau_wlf
   *   or Carreau_Suspension etc
   */
  int ConstitutiveEquation;

  dbl mu0;
  int mu0Model;
  dbl pos_ls_mup;
  int len_u_mu0;
  dbl *u_mu0;
  dbl nexp;
  int nexpModel;
  int len_u_nexp;
  dbl *u_nexp;
  dbl muinf;
  int muinfModel;
  int len_u_muinf;
  dbl *u_muinf;
  dbl lam;
  int lamModel;
  int len_u_lam;
  dbl *u_lam;
  dbl aexp;
  int aexpModel;
  int len_u_aexp;
  dbl *u_aexp;
  dbl atexp;
  int atexpModel;
  int len_u_atexp;
  dbl *u_atexp;
  /* CARREAU_WLF viscosity model  */
  dbl wlfc2;
  int wlfc2Model;
  int len_u_wlfc2;
  dbl *u_wlfc2;
  /* these are for the BINGHAM and HERSCHEL_BULKLEY
   *  yielding material model */
  dbl tau_y;
  int len_u_tau_y;
  dbl *u_tau_y;
  int tau_yModel;
  dbl fexp;
  int fexpModel;
  dbl epsilon;
  int epsilonModel;
  /* these are for SUSPENSION/FILLED_EPOXY models */
  dbl maxpack;
  int maxpackModel;
  int sus_species_no;
  /* these are for CURE/EPOXY/FILLED_EPOXY models */
  dbl gelpoint;
  int gelpointModel;
  dbl cureaexp;
  int cureaexpModel;
  dbl curebexp;
  int curebexpModel;
  dbl tgel0;
  int tgel0Model;
  int cure_species_no;
  dbl k1; /* rate coefficients etc. for Bond evolution viscosity model */
  dbl k2;
  dbl n0;
  dbl pexp;
  dbl qexp;
  dbl diff;
  //! Model for the dilational viscosity
  int DilViscModel;
  //! This is the constant used for the dilational viscosity
  dbl DilVisc0;
  dbl thixo_factor;
  int thixoModel;
  int len_u_thixo;
  dbl *u_thixo_factor;
};
typedef struct Generalized_Newtonian GEN_NEWT_STRUCT;
typedef struct PolymerTimeConstants {
  //! Integer describing the polymer time constant model
  /*!
   *   Constant, Power Law, Carreau, Bingham,  Carreau_wlf
   *   or Carreau_Suspension etc
   */
  int ConstitutiveEquation;

  dbl lambda0;
  int lambda0Model;
  dbl pos_ls_lambda;
  dbl nexp;
  int nexpModel;
  dbl lambdainf;
  int lambdainfModel;
  dbl carreau_lambda;
  int carreau_lambdaModel;
  dbl aexp;
  int aexpModel;
  dbl atexp;
  int atexpModel;
} POLYMER_TIME_CONST_STRUCT;

struct Positive_LS_Viscoelastic_Properties {

  double alpha; /* This is the Geisekus mobility parameter */

  double xi; /* This is the PTT upper convected / lower convected weight parameter */

  double eps; /* This is the PTT elongational parameter */
};

struct Viscoelastic_Constitutive {
  /* this struct contains the polymer viscosity
   * if it is shearthinning etc or NEWTONIAN
   */
  GEN_NEWT_STRUCT *gn;
  POLYMER_TIME_CONST_STRUCT *time_const_st;

  dbl alpha; /* This is the Geisekus mobility parameter */
  int alphaModel;

  dbl xi; /* This is the PTT upper convected / lower convected weight parameter */
  int xiModel;

  dbl eps; /* This is the PTT elongational parameter */
  int epsModel;

  // Rolie Poly
  dbl stretch_time;
  int stretchModel;

  dbl CCR_coefficient;
  int CCR_coefficientModel;

  dbl polymer_exponent;
  int polymer_exponentModel;

  dbl maximum_stretch_ratio;
  int maximum_stretch_ratioModel;

  // FENE
  dbl extensibility;
  int extensibilityModel;

  // level set
  struct Positive_LS_Viscoelastic_Properties pos_ls;
  dbl muJeffreys; /* 2nd viscosity used in modified Jeffreys model */
  int muJeffreysModel;
};
typedef struct Viscoelastic_Constitutive VISC_CONST_STRUCT;

struct Viscoelastic_Nonmodal {
  /*
   * Viscoelastic Constitutive Equation:
   * with proper coefficient choices it can become:
   * Giesekus Model
   * Maxwell Model
   * Oldroyd-B Model
   * White-Metzner Model
   * Leonov Model
   */
  int ConstitutiveEquation;
  dbl wt_func;
  int wt_funcModel;
  dbl shockcapture;
  int shockcaptureModel;

  int ptt_type;
  /* This is the adaptive viscosity scaling. If if it zero
   * we get the standard formulation, if nonzero we get
   * numerical viscosity that may stabilize the stress equations
   */
  dbl eps;

  int evssModel; /* this is to choose the EVSS model - either
                    Fortin's or Rajagopalans */
  int modes;

  int dg_J_model;

  dbl *dg_J_model_wt;
  int len_dg_J_model_wt; /* Sigh...everyone forgets to do this... */

  int shiftModel;
  dbl *shift;
  int len_shift; /*  time constant temperature shift */
};

struct Elastic_Constitutive {
  /*
   * Constants used in the Elasticity consitutive Equations
   */
  int ConstitutiveEquation;

  dbl lame_mu;
  int lame_mu_model;
  int len_u_mu;
  dbl *u_mu;
  dbl d_lame_mu[MAX_VARIABLE_TYPES + MAX_CONC];
  int lame_mu_tableid;

  dbl lame_lambda;
  int lame_lambda_model;
  int len_u_lambda;
  dbl *u_lambda;
  dbl d_lame_lambda[MAX_VARIABLE_TYPES + MAX_CONC];

  dbl lame_TempShift;
  int lameTempShiftModel;
  int len_u_lame_TempShift;
  dbl *u_lame_TempShift;
  dbl d_lame_TempShift[MAX_VARIABLE_TYPES + MAX_CONC];
  int lame_TempShift_tableid;

  dbl bend_stiffness;
  int bend_stiffness_model;
  int len_u_bend_stiffness;
  dbl *u_bend_stiffness;
  dbl d_bend_stiffness[MAX_VARIABLE_TYPES + MAX_CONC];

  dbl exten_stiffness;
  int exten_stiffness_model;
  int len_u_exten_stiffness;
  dbl *u_exten_stiffness;
  dbl d_exten_stiffness[MAX_VARIABLE_TYPES + MAX_CONC];

  dbl poisson;
  int poisson_model;
  int len_u_poisson;
  dbl *u_poisson;
  dbl d_poisson[MAX_VARIABLE_TYPES + MAX_CONC];

  dbl Strss_fr_sol_vol_frac;

  dbl v_mesh_sfs[DIM];
  int v_mesh_sfs_model; /* Looks as if this ought to be an int - pas */
  int len_u_v_mesh_sfs;
  dbl *u_v_mesh_sfs;

  dbl thermal_expansion; /*  thermo-elasticity properties   */
  int thermal_expansion_model;
  int len_u_thermal_expansion;
  dbl *u_thermal_expansion;
  int thermal_expansion_tableid;
  dbl solid_reference_temp;
  int solid_reference_temp_model;

  dbl solid_viscosity; /*  viscoelastic solid viscosity   */
  int solid_viscosity_model;
  int len_u_solid_viscosity;
  dbl *u_solid_viscosity;
  dbl solid_dil_viscosity; /*  viscoelastic solid dilational viscosity   */
  int solid_dil_viscosity_model;
  int len_u_solid_dil_viscosity;
  dbl *u_solid_dil_viscosity;
  dbl solid_retardation; /*  viscoelastic solid retardation time   */
  int solid_retard_model;
  int len_u_solid_retardation;
  dbl *u_solid_retardation;
};

typedef struct Elastic_Constitutive ELASTIC_CONST_STRUCT;

struct Viscoplastic_Constitutive {
  /*
   * Constants used in the Viscoplasticity consitutive Equations
   */
  int ConstitutiveEquation;

  int update_flag; /*This is a flag for updating and advancing the
                    *time integrals in the EVP hyperelastic formulation.
                    *We have no other place to put this flag unless we
                    *want to pass it into the depths of assemble_mesh. Because
                    *it is only used there, I saw it not necessary to cludder
                    *matrix fill arg list, and I made it global.
                    */

  dbl plastic_mu;
  int plastic_mu_model;
  int len_u_plastic_mu;
  dbl *u_plastic_mu;
  dbl d_plastic_mu[MAX_VARIABLE_TYPES + MAX_CONC];

  dbl yield;
  int yield_model;
  int len_u_yield;
  dbl *u_yield;
  dbl d_yield[MAX_VARIABLE_TYPES + MAX_CONC];

  /*
   * This part of this struct is used to allocate, transport, and deallocate
   * the global arrays
   * required for elastoviscoplastic models.  Those models, because of their Lagrangian
   * hyperelastic nature require a time integral along material paths.  The only
   * way to track those time integrals is to advanced these kinematic tensors globally,
   * i.e., they can not be in an element-by-element structure.  These babies can get
   * pretty big so they won't be allocated unless you are doing an EVP model
   */

  dbl ****F_vp_glob; /*viscoplastic deformation gradient = I + int(D_vp)dt */
  dbl ****F_vp_old_glob;
  dbl ****TT_glob; /*Plastic potential stress */
  dbl ****TT_old_glob;
  dbl ******dTT_dx_glob; /*Stress sensitivities */
  dbl ******dTT_dx_old_glob;
  dbl ******dTT_dc_glob; /*Stress sensitivities */
  dbl ******dTT_dc_old_glob;
};

/*
 * This is where this belongs.
 */

struct Variable_Initialization {
  int var;
  int ktype;
  double init_val;
  double init_val_min;
  double init_val_max;
  double init_val_minus; /* Value in negative LS phase */
  double init_val_plus;  /* Value in positive LS phase */
  int slave_block;       /* this is set TRUE in order that this initialization value
                           is not applied to nodes that are shared with adjacent blocks */
  int len_u_pars;
  double *u_pars;
};

struct Second_LS_Phase_Properties {
  int ViscosityModel;
  dbl viscosity;
  int viscositymask[2];
  dbl viscosity_phase[MAX_PHASE_FUNC];

  int DensityModel;
  dbl density;
  int densitymask[2];
  dbl density_phase[MAX_PHASE_FUNC];

  int HeatCapacityModel;
  dbl heatcapacity;
  int heatcapacitymask[2];
  dbl heatcapacity_phase[MAX_PHASE_FUNC];

  int ThermalConductivityModel;
  dbl thermalconductivity;
  int thermalconductivitymask[2];
  dbl thermalconductivity_phase[MAX_PHASE_FUNC];

  int MomentumSourceModel;
  dbl momentumsource[DIM];
  int momentumsourcemask[2];
  dbl momentumsource_phase[MAX_PHASE_FUNC][DIM];

  int HeatSourceModel;
  dbl heatsource;
  int heatsourcemask[2];
  dbl heatsource_phase[MAX_PHASE_FUNC];

  int AcousticImpedanceModel;
  dbl acousticimpedance;
  int acousticimpedancemask[2];
  dbl acousticimpedance_phase[MAX_PHASE_FUNC];

  int wavenumberModel;
  dbl wavenumber;
  int wavenumbermask[2];
  dbl wavenumber_phase[MAX_PHASE_FUNC];

  int AcousticAbsorptionModel;
  dbl acousticabsorption;
  int acousticabsorptionmask[2];
  dbl acousticabsorption_phase[MAX_PHASE_FUNC];

  int RefractiveIndexModel;
  dbl refractiveindex;
  int refractiveindexmask[2];
  dbl refractiveindex_phase[MAX_PHASE_FUNC];

  int LightAbsorptionModel;
  dbl lightabsorption;
  int lightabsorptionmask[2];
  dbl lightabsorption_phase[MAX_PHASE_FUNC];

  int ExtinctionIndexModel;
  dbl extinctionindex;
  int extinctionindexmask[2];
  dbl extinctionindex_phase[MAX_PHASE_FUNC];

  int SpeciesSourceModel[MAX_CONC];
  dbl speciessource[MAX_CONC];
  int speciessourcemask[2][MAX_CONC];
  dbl speciessource_phase[MAX_PHASE_FUNC][MAX_CONC];
  int use_species_source_width[MAX_CONC];
  dbl species_source_width[MAX_CONC];

  int FlowingLiquidViscosityModel;
  dbl FlowingLiquid_viscosity;
  int FlowingLiquid_viscositymask[2];
  dbl FlowingLiquid_viscosity_phase[MAX_PHASE_FUNC];
};

typedef struct Second_LS_Phase_Properties SECOND_LS_PHASE_PROP_STRUCT;

#endif

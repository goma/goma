/************************************************************************\
 * Copyright (c) 2002 Sandia Corporation.                               *
 *					         			*
 * Under the terms of Contract DE-AC04-94AL85000, there is a            *
 * non-exclusive license for use of this work by or on behalf of the    *
 * U.S. Government. Export of this program may require a license from   *
 * the United States Government.                                        *
 *					         			*
 * This software is the property of Sandia Corporation and discloses    *
 * material protectable under copyright laws of the United States.      *
 * Use, Duplication, or Disclosure is prohibited, unless allowed        *
 * subject to the terms of a separate license agreement.                *
 *					         			*
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES          *
 * DEPARTMENT OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR       *
 * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY    *
 * LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,    *
 * OR USEFULNESS OF ANY INFORMATION, APPARATUS OR PROCESS DISCLOSED,    *
 * OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.*
 *					         			*
\************************************************************************/

%{


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

#include "mm_parser.tab.h"
#include "mm_parser.h"
#include "std.h"
#include "mm_input.h"

#undef YY_INPUT
#define YY_INPUT(buf,result,max) (result = myinput(buf,max))


/*#undef yyinput*/
/*#undef yyunput*/

/*#define input() (*mystring++)*/
/*#define unput(c) (*--mystring = c)*/

/*#define yyinput(void) (return input(); )*/
/*#define yyunput(int c) (unput(c);)*/

extern char *mystring;
extern parser_file_mode;
extern int line_number /*, first_tab_flag*/ ;
extern FILE *parser_log, *error_log;
extern char debug_clue[100];
extern int error_found_in_last_line;
extern int error_number;
extern int ProcID;
YY_BUFFER_STATE hold_buffer;


char	msg[100];
int count_of_characters_in_this_line = 0;

void print_to_log();

%}
%e 200000
%p 300000
%a 200000
%n 50000
%k 4000

%%
%{
if (parser_file_mode == input_file) 		{parser_file_mode = reset; return PARSER_IN_INPUT_FILE_MODE;}
if (parser_file_mode == mat_file) 		{parser_file_mode = reset; return PARSER_IN_MAT_FILE_MODE;}

%}

"---Mechanical Properties and Constitutive Equations"	{print_to_log(); /*return(MECHANICAL_PROPERTIES_DELIMITER);*/}
"---Physical Properties"		{print_to_log(); /*return(PHYSICAL_PROPERTIES_DELIMITER);*/}
"="					{print_to_log(); strcpy(yylval.string, yytext);   return(EQUALS_);}
"3d"|"3D"				{print_to_log(); strcpy(yylval.string, yytext);   return(THREE_D_);}
"3dfile"|"3DFILE"			{print_to_log(); strcpy(yylval.string, yytext);   return(THREE_D_FILE_);}
"STEADY"|"steady"|"Steady"		{print_to_log(); strcpy(yylval.string, yytext);   return(STEADY_);}
"TRANSIENT"|"Transient"|"transient"	{print_to_log(); strcpy(yylval.string, yytext);   return(TRANSIENT_);}
"A" 					{print_to_log(); strcpy(yylval.string, yytext);   return(A_);}
"MAX."|"Max."|"max."			{print_to_log(); strcpy(yylval.string, yytext);   return(MAX_);}
"DESIRED"|"Desired"|"desired"		{print_to_log(); strcpy(yylval.string, yytext);   return(DESIRED_);}
"ALC"|"Alc"|"alc"			{print_to_log(); strcpy(yylval.string, yytext);   return(ALC_);}
"SENSITIVITY"|"Sensitivity"|"sensitivity"	{print_to_log(); strcpy(yylval.string, yytext);   return(SENSITIVITY_);}
"AGGRESSIVENESS"|"aggressiveness"|"Aggressiveness"	{print_to_log(); strcpy(yylval.string, yytext);   return(AGGRESSIVENESS_);}
"Anneal"|"ANNEAL"			{print_to_log(); strcpy(yylval.string, yytext);   return(ANNEAL_);}
"on"|"On"|"ON"				{print_to_log(); strcpy(yylval.string, yytext);   return(ON_);}
"off"|"Off"|"OFF"			{print_to_log(); strcpy(yylval.string, yytext);   return(OFF_);}
"AC"|"ac"				{print_to_log(); strcpy(yylval.string, yytext);   return(AC_);}
"MT"|"mt"				{print_to_log(); strcpy(yylval.string, yytext);   return(MT_);}
"ARBITRARY" 				{print_to_log(); strcpy(yylval.string, yytext);   return(ARBITRARY_);}
"Adaptive"				{print_to_log(); strcpy(yylval.string, yytext);   return(ADAPTIVE_);}
"Aexp" 					{print_to_log(); strcpy(yylval.string, yytext);   return(AEXP_);}
"Algorithm"|"algorithm"			{print_to_log(); strcpy(yylval.string, yytext);   return(ALGORITHM_);}  
"And"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(AND_);} 
"Augmenting"|"augmenting" 		{print_to_log(); strcpy(yylval.string, yytext);   return(AUGMENTING_);} 
"B"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(B_);}
"USERBC"|"userbc" 			{print_to_log(); strcpy(yylval.string, yytext);   return(USERBC_);}
"USERMAT"|"usermat" 			{print_to_log(); strcpy(yylval.string, yytext);   return(USERMAT_);}
"BASIS_ONCE"				{print_to_log(); strcpy(yylval.string, yytext);   return(BASIS_ONCE_);}
"BASIS_RESEED"				{print_to_log(); strcpy(yylval.string, yytext);   return(BASIS_RESEED_);}
"BC"|"bc"|"Bc"				{print_to_log(); strcpy(yylval.string, yytext);   return(BC_);} 
"MATERIAL"|"Material"|"material"	{print_to_log(); strcpy(yylval.string, yytext);   return(MATERIAL_);}
"TAG"|"Tag"|"tag"			{print_to_log(); strcpy(yylval.string, yytext);   return(TAG_);} 
"BINGHAM"				{print_to_log(); strcpy(yylval.string, yytext);   return(BINGHAM_);}
"BOUSS"					{print_to_log(); strcpy(yylval.string, yytext);   return(BOUSS_);}
"Body"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(BODY_);}
"Boundary"|"BOUNDARY"|"boundary"	{print_to_log(); strcpy(yylval.string, yytext);   return(BOUNDARY_);}
"Bulk"|"bulk"				{print_to_log(); strcpy(yylval.string, yytext);   return(BULK_);}
"BUTLER_VOLMER"				{print_to_log(); strcpy(yylval.string, yytext);   return(BUTLER_VOLMER_);}
"CC" 					{print_to_log(); strcpy(yylval.string, yytext);   return(CC_);}
"TP"|"tp"|"Tp"				{print_to_log(); strcpy(yylval.string, yytext);   return(TP_);}
"TC"|"tc"|"Tc"				{print_to_log(); strcpy(yylval.string, yytext);   return(TC_);}
"CA"|"CA_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(CA_);}
"HC"|"hc"				{print_to_log(); strcpy(yylval.string, yytext);   return(HC_);}
"CAPILLARY"|"CAPILLARY_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(CAPILLARY_);}
"CAP_REPULSE"|"CAP_PRESSURE_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(CAP_REPULSE_);}
"CARREAU"				{print_to_log(); strcpy(yylval.string, yytext);   return(CARREAU_);}
"CARTESIAN"				{print_to_log(); strcpy(yylval.string, yytext);   return(CARTESIAN_);}
"CA_EDGE"|"CA_EDGE_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(CA_EDGE_);}
"CA_EDGE_CURVE"|"CA_EDGE_CURVE_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(CA_EDGE_CURVE_);}
"CA_EDGE_CURVE_INT"|"CA_EDGE_CURVE_INT_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(CA_EDGE_CURVE_INT_);}
"CA_EDGE_INT"|"CA_EDGE_INT_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(CA_EDGE_INT_);}
"CA_EDGE_OR_FIX"|"CA_EDGE_OR_FIX_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(CA_EDGE_OR_FIX_);}
"CA_OR_FIX"|"CA_OR_FIX_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(CA_OR_FIX_);}
"CIRCLE" 				{print_to_log(); strcpy(yylval.string, yytext);   return(CIRCLE_);}
"COLLOCATE_EDGE" 			{print_to_log(); strcpy(yylval.string, yytext);   return(COLLOCATE_EDGE_);}
"COMPRESSIBLE"				{print_to_log(); strcpy(yylval.string, yytext);   return(COMPRESSIBLE_);}
"CONDITIONS"|"Conditions"|"conditions"	{print_to_log(); strcpy(yylval.string, yytext);   return(CONDITIONS_);}
"CONSTANT_LAST_CONC"|"Constant_Last_Conc"|"constant_last_conc"			{print_to_log(); strcpy(yylval.string, yytext);   return(CONSTANT_LAST_CONC_);}
"CONSTANT"|"Constant"|"constant"	{print_to_log(); strcpy(yylval.string, yytext);   return(CONSTANT_);}
"CONTACT_LINE"				{print_to_log(); strcpy(yylval.string, yytext);   return(CONTACT_LINE_);}
"CONT_NORM_VEL"|"CONT_NORM_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(CONT_NORM_VEL_);}
"CONT_TANG_VEL"|"CONT_TANG_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(CONT_TANG_VEL_);}
"CURE"					{print_to_log(); strcpy(yylval.string, yytext);   return(CURE_);}
"HZERO"|"Hzero"|"hzero"			{print_to_log(); strcpy(yylval.string, yytext);   return(HZERO_);}
"HFIRST"|"Hfirst"|"first"		{print_to_log(); strcpy(yylval.string, yytext);   return(HFIRST_);}
"CURRENT"				{print_to_log(); strcpy(yylval.string, yytext);   return(CURRENT_);} 
"Capacity"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(CAPACITY_);}
"Capilary" 				{print_to_log(); strcpy(yylval.string, yytext);   return(CAPILLARY_);}
"Coefficient" 				{print_to_log(); strcpy(yylval.string, yytext);   return(COEFFICIENT_);}
"Condition"|"CONDITION"|"condition"	{print_to_log(); strcpy(yylval.string, yytext);   return(CONDITION_);} 
"Conduction"			 	{print_to_log(); strcpy(yylval.string, yytext);   return(CONDUCTION_);} 
"Conductivity" 				{print_to_log(); strcpy(yylval.string, yytext);   return(CONDUCTIVITY_);}
"Constant2" 				{print_to_log(); strcpy(yylval.string, yytext);   return(CONSTANT_2_);}
"Constitutive"			 	{print_to_log(); strcpy(yylval.string, yytext);   return(CONSTITUTIVE_);}
"Contours"|"CONTOURS"|"contours"	{print_to_log(); strcpy(yylval.string, yytext);   return(CONTOURS_);}
"Convection" 				{print_to_log(); strcpy(yylval.string, yytext);   return(CONVECTION_);}  
"Convective"				{print_to_log(); strcpy(yylval.string, yytext);   return(CONVECTIVE_);}
"Coordinate"	  			{print_to_log(); strcpy(yylval.string, yytext);   return(COORDINATE_);}
"Cure" 					{print_to_log(); strcpy(yylval.string, yytext);   return(CURE_);} 
"Current"|"CURRENT"|"CURRENT_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(CURRENT_);}
"D1"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(D1_);}
"D1_RS"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(D1_RS_);}
"D2"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(D2_);}
"D2_RS"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(D2_RS_);}
"D3"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(D3_);}
"D3_RS"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(D3_RS_);}
"DARCEY"				{print_to_log(); strcpy(yylval.string, yytext);   return(DARCEY_);}
"DATA"|"Data"|"data"			{print_to_log(); strcpy(yylval.string, yytext);   return(DATA_);}
"DARCY_CONTINUOUS"|"DARCY_CONTINUOUS_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(DARCY_CONTINUOUS_);}
"DARCY_FICKIAN"				{print_to_log(); strcpy(yylval.string, yytext);   return(DARCEY_FICKIAN_);}
"DEFORM"				{print_to_log(); strcpy(yylval.string, yytext);   return(DEFORM_);}
"DILATATION"|"Dilatation"|"dilatation"	{print_to_log(); strcpy(yylval.string, yytext);   return(DILATATION_);}
"Dilation"|"DILATION"|"dilation"	{print_to_log(); strcpy(yylval.string, yytext);   return(DILATION_);} 
"DISCONTINUOUS_VELO"|"DISCONTINUOUS_VELO_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(DISCONTINUOUS_VELO_);}
"DISTNG"|"DISTNG_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(DISTNG_);}
"DMX" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(DMY_RS_);}	 
"DMX" 	 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DMX_);}	 
"DMY" 	 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DMY_);}	 
"DMY_RS" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DMY_RS_);}	 
"DMZ" 	 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DMZ_);}	 
"DMZ_RS" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DMZ_RS_);}	 
"DN" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(DN_);}	
"Domain"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DOMAIN_);}	 
"DT1" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(DT1_);}	 
"DT2" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(DT2_);}	 
"DX"|"DX_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(DX_);} 
"DXDISTNG"|"DXDISTNG_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(DXDISTNG_);} 
"DX_RS"|"DX_RS_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(DX_RS_);} 
"DY"|"DY_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(DY_);} 
"DYDISTNG"|"DYDSITNG_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(DYDISTNG_);} 
"DY_RS"|"DY_RS_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(DY_RS_);} 
"DZ"|"DZ_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(DZ_);} 
"DZDISTNG"|"DZDISTNG_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(DZDISTNG_);} 
"DZ_RS"|"DZ_RS_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(DZ_RS_);} 
"Database"|"DATABASE"|"database"	{print_to_log(); strcpy(yylval.string, yytext);   return(DATABASE_);}
"Dumeps"				{print_to_log(); strcpy(yylval.string, yytext);   return(DUMPS_);}
"Datum"|"DATUM"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DATUM_);} 
"Debug" 				{print_to_log(); strcpy(yylval.string, yytext);   return(DEBUG_);} 
"Default"|"default"|"DEFAULT"		{print_to_log(); strcpy(yylval.string, yytext);   return(DEFAULT_);}	
"Density" 				{print_to_log(); strcpy(yylval.string, yytext);   return(DENSITY_);}
"Derivatives"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DERIVATIVES_);}  
"Description"				{print_to_log(); strcpy(yylval.string, yytext);   return(DESCRIPTION_);}	
"Diffusion"			 	{print_to_log(); strcpy(yylval.string, yytext);   return(DIFFUSION_);}
"Diffusivity" 				{print_to_log(); strcpy(yylval.string, yytext);   return(DIFFUSIVITY_);} 
"Disable""|disable"|"DISABLE"		{print_to_log(); strcpy(yylval.string, yytext);   return(DISABLE_);}
"Discontinuous"				{print_to_log(); strcpy(yylval.string, yytext);   return(DISCONTINUOUS_);}
"Divergence"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(DIVERGENCE_);} 
"ELECTROOSMOTIC"                        {print_to_log(); strcpy(yylval.string, yytext);   return(ELECTROOSMOTIC_);}
"END"|"End"|"end"			{print_to_log(); strcpy(yylval.string, yytext);   return(END_);}  
"ENTHALPY"				{print_to_log(); strcpy(yylval.string, yytext);   return(ENTHALPY_);}
"EPOXY"					{print_to_log(); strcpy(yylval.string, yytext);   return(EPOXY_);}
"EQ"  					{print_to_log(); strcpy(yylval.string, yytext);   return(EQ_);}
"EVP"					{print_to_log(); strcpy(yylval.string, yytext);   return(EVP_);}
"EVSS_F"				{print_to_log(); strcpy(yylval.string, yytext);   return(EVSS_F_);}
"EVSS_G"				{print_to_log(); strcpy(yylval.string, yytext);   return(EVSS_G_);}
"EXODUS"				{print_to_log(); strcpy(yylval.string, yytext);   return(EXODUS_);}
"EXPONENTIAL"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(EXPONENTIAL_);} 
"Electrical"|"electrical"|"ELECTRICAL"	{print_to_log(); strcpy(yylval.string, yytext);   return(ELECTRICAL_);}
"ELECTRONEUTRALITY_FICKIAN"		{print_to_log(); strcpy(yylval.string, yytext);   return(ELECTRONEUTRALITY_FICKIAN_);}
"ELECTRONEUTRALITY_SM"			{print_to_log(); strcpy(yylval.string, yytext);   return(ELECTRONEUTRALITY_SM_);}
"Element"  				{print_to_log(); strcpy(yylval.string, yytext);   return(ELEMENT_);} 
"Epsilon" 				{print_to_log(); strcpy(yylval.string, yytext);   return(EPSILON_);}
"Equation"|"EQUATION"|"equation"	{print_to_log(); strcpy(yylval.string, yytext);   return(EQUATION_);} 
"Equations"			 	{print_to_log(); strcpy(yylval.string, yytext);   return(EQUATIONS_);}
"Expansion"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(EXPANSION_);}
"Exponent"|"EXPONENT"|"exponent"	{print_to_log(); strcpy(yylval.string, yytext);   return(EXPONENT_);}
"External"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(EXTERNAL_);} 
"F"|"F_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(F_);}
"FEM"					{print_to_log(); strcpy(yylval.string, yytext);   return(FEM_);}
"FICKIAN"				{print_to_log(); strcpy(yylval.string, yytext);   return(FICKIAN_);}
"FINAL"|"Final"|"final"				{print_to_log(); strcpy(yylval.string, yytext);   return(FINAL_);}
"FILL"|"Fill"|"fill"			{print_to_log(); strcpy(yylval.string, yytext);   return(FILL_);}
"FILLED_EPOXY"				{print_to_log(); strcpy(yylval.string, yytext);   return(FILLED_EPOXY_);}
"FILL_CA"|"FILL_CA_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(FILL_CA_);}
"FILL_INLET"|"FILL_INLET_CA"		{print_to_log(); strcpy(yylval.string, yytext);   return(FILL_INLET_);}
"FIX"|"FIX_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(FIX_);}
"FLORY"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(FLORY_);}
"DELTA_S"|"Delta_s"|"delta_s"		{print_to_log(); strcpy(yylval.string, yytext);   return(DELTA_S_);}
"FLOW_GRADV"|"FLOW_GRADV_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(FLOW_GRADV_);}
"FLOW_HYDROSTATIC"|"FLOW_HYDROSTATIC_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(FLOW_HYDROSTATIC_);}
"FLOW_PRESSURE"|"FLOW_PRESSURE_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(FLOW_PRESSURE_);}
"FLOW_PRESS_USER"|"FLOW_PRESS_USER_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(FLOW_PRESS_USER_);}
"FLOW_STRESSNOBC"|"FLOW_STRESSNOBC_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(FLOW_STRESSNOBC_);}
"FLUID_SOLID"|"FLUID_SOLID_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(FLUID_SOLID_);}
"FLUID_SOLID_RS"|"FLUID_SOLID_RS_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(FLUID_SOLID_RS_);}
"FORCE"|"FORCE_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_);} 
"FORCE_RS"|"FORCE_RS_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_RS_);} 
"FORCE_USER"|"FORCE_USER_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_USER_);}
"FORCE_USER_RS"|"FORCE_USER_RS_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_USER_RS_);}
"FREE_VOL"				{print_to_log(); strcpy(yylval.string, yytext);   return(FREE_VOL_);}
"Field" 				{print_to_log(); strcpy(yylval.string, yytext);   return(FIELD_);}
"File"|"file"|"FILE"			{print_to_log(); strcpy(yylval.string, yytext);   return(FILE_);}
"inline"|"Inline"|"INLINE"		{print_to_log(); strcpy(yylval.string, yytext);   return(INLINE_);}
"First"|"first"|"FIRST"			{print_to_log(); strcpy(yylval.string, yytext);   return(FIRST_);} 
"Flowing"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(FLOWING_);}
"Fluctuations"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(FLUCTUATIONS_);}  
"Fluxiness" 				{print_to_log(); strcpy(yylval.string, yytext);   return(FLUXINESS_);}  
"Fluxlines"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(FLUXLINES_);} 
"Formulation"				{print_to_log(); strcpy(yylval.string, yytext);   return(FORMULATION_);}
"Frac"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(FRAC_);}
"Free"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(FREE_);}
"Frequency"|"FREQUENCY"|"frequency"	{print_to_log(); strcpy(yylval.string, yytext);   return(FREQUENCY_);}
"Function"				{print_to_log(); strcpy(yylval.string, yytext);   return(FUNCTION_);}
"Fusion" 				{print_to_log(); strcpy(yylval.string, yytext);   return(FUSION_);}
"G11"|"G11_BC" 				{print_to_log(); strcpy(yylval.string, yytext);   return(G11_);}
"G12"|"G12_BC"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(G12_);}	 
"G13"|"G13_BC" 				{print_to_log(); strcpy(yylval.string, yytext);   return(G13_);}
"G21"|"G21_BC"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(G21_);}	 
"G22"|"G22_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(G22_);}
"G23"|"G23_BC"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(G23_);}	 
"G31"|"G31_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(G31_);}
"G32"|"G32_BC"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(G32_);}	 
"G33"|"G33_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(G33_);}
"GALERKIN"				{print_to_log(); strcpy(yylval.string, yytext);   return(GALERKIN_);}
"GD" 					{print_to_log(); strcpy(yylval.string, yytext);   return(GD_);}
"GD_CIRC"|"GD_CIRC_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(GD_CIRC_);}
"GD_CONST"|"GD_CONST_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(GD_CONST_);}
"GD_LINEAR"|"GD_LINEAR_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(GD_LINEAR_);}
"GD_PARAB"|"GD_PARAB_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(GD_PARAB_);}
"GD_POLYN"|"GD_POLYN_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(GD_POLYN_);}
"GD_TABLE"|"GD_TABLE_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(GD_TABLE_);}
"GD_TIME"|"GD_TIME_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(GD_TIME_);}
"GEMO"|"GEOM_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(GEOM_);}
"GEMOX"|"GEOMX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(GEOMX_);}
"GEMOY"|"GEOMY_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(GEOMY_);}
"GEMOZ"|"GEOMZ_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(GEOMZ_);}
"GIESEKUS"				{print_to_log(); strcpy(yylval.string, yytext);   return(GIESEKUS_);}
"GUESS"|"Guess"|"guess"			{print_to_log(); strcpy(yylval.string, yytext);   return(GUESS_);}
"Gas" 					{print_to_log(); strcpy(yylval.string, yytext);   return(GAS_);} 
"Gel"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(GEL_);}
"General"				{print_to_log(); strcpy(yylval.string, yytext);   return(GENERAL_);}
"HOOKEAN_PSTRAIN"			{print_to_log(); strcpy(yylval.string, yytext);   return(HOOKEAN_PSTRAIN_);}
"HOOKEAN_PSTRESS"			{print_to_log(); strcpy(yylval.string, yytext);   return(HOOKEAN_PSTRESS_);}
"HYDRODYNAMIC"				{print_to_log(); strcpy(yylval.string, yytext);   return(HYDRODYNAMIC_);}
"HYDROSTATIC_SYMM"|"HYDROSTATIC_SYMM_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(HYDROSTATIC_SYMM_);} 
"Heat" 					{print_to_log(); strcpy(yylval.string, yytext);   return(HEAT_);}
"High" 					{print_to_log(); strcpy(yylval.string, yytext);   return(HIGH_);}
"IDIM"|"idim"|"Idim"			{print_to_log(); strcpy(yylval.string, yytext);   return(IDIM_);}
"II"					{print_to_log(); strcpy(yylval.string, yytext);   return(II_);}
"ID"|"id"|"Id"				{print_to_log(); strcpy(yylval.string, yytext);   return(ID_);}
"INCOMP_STRAIN"				{print_to_log(); strcpy(yylval.string, yytext);   return(INCOMP_PSTRAIN_);}
"INCOMP_STRESS"				{print_to_log(); strcpy(yylval.string, yytext);   return(INCOMP_PSTRESS_);}
"Inertia"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(INERTIA_);}
"Initial"|"INITIAL"|"initial"		{print_to_log(); strcpy(yylval.string, yytext);   return(INITIAL_);}  
"Initialize" 				{print_to_log(); strcpy(yylval.string, yytext);   return(INITIALIZE_);} 
"Integration"|"integration"|"INTEGRATION" {print_to_log(); strcpy(yylval.string, yytext); return(INTEGRATION_);}
"Invariant"|"INVARIANT"|"invariant"|"INVARIENT"|"Invarient"|"invarient"	{print_to_log(); strcpy(yylval.string, yytext);   return(INVARIANT_);} 
"Iterations"|"ITERATIONS"|"iterations"	{print_to_log(); strcpy(yylval.string, yytext);   return(ITERATIONS_);} 
"Jacobi"|"jacobi"|"JACOBI"		{print_to_log(); strcpy(yylval.string, yytext);   return(JACOBI_);}
"Bjacobi"|"bjacobi"|"BJACOBI"		{print_to_log(); strcpy(yylval.string, yytext);   return(BJACOBI_);}
"Jacobian"|"JACOBIAN"|"jacobian"	{print_to_log(); strcpy(yylval.string, yytext);   return(JACOBIAN_);}
"KINEMATIC"|"KINEMATIC_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(KINEMATIC_);} 
"KINEMATIC_COLLOC"|"KINEMATIC_COLLOC_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(KINEMATIC_COLLOC_);}
"KINEMATIC_DISC"|"KINEMATIC_DISC_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(KINEMATIC_DISC_);}
"KINEMATIC_EDGE"|"KINEMATIC_EDGE_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(KINEMATIC_EDGE_);}
"KINEMATIC_PETROV"|"KINEMATIC_PETROV_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(KINEMATIC_PETROV_);} 
"KINEMATIC_SPECIES"|"KINEMATIC_SPECIES_BC"  {print_to_log(); strcpy(yylval.string, yytext);   return(KINEMATIC_SPECIES_);}
"KIN_CHEM"|"KIN_CHEM_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(KIN_CHEM_);}
"KIN_DISPLACEMENT"|"KIN_DISPLACEMENT_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(KIN_DISPLACEMENT_);}
"KIN_LEAK"|"KIN_LEAK_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(KIN_LEAK_);} 
"Krylov"|"krylov"|"KRYLOV" 		{print_to_log(); strcpy(yylval.string, yytext);   return(KRYLOV_);} 
"LAMBDA" 				{print_to_log(); strcpy(yylval.string, yytext);   return(LAMBDA_);}
"LATENT_HEAT"|"LATENT_HEAT_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(LATENT_HEAT_);} 
"LATENT_HEAT_INTERNAL"|"LATENT_HEAT_INTERNAL_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(LATENT_HEAT_INTERNAL_);}
"LINEAR"|"Linear"|"linear"		{print_to_log(); strcpy(yylval.string, yytext);   return(LINEAR_);}
"Lagrangian"|"LAGRANGIAN"		{print_to_log(); strcpy(yylval.string, yytext);   return(LAGRANGIAN_);}
"Lame"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(LAME_);}
"Latent"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(LATENT_);}
"Law" 					{print_to_log(); strcpy(yylval.string, yytext);   return(LAW_);}
"Linear"			 	{print_to_log(); strcpy(yylval.string, yytext);   return(LINEAR_);} 
"Liq"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(LIQ_);}
"Liquid"			 	{print_to_log(); strcpy(yylval.string, yytext);   return(LIQUID_);}
"Liquidus" 				{print_to_log(); strcpy(yylval.string, yytext);   return(LIQUIDUS_);}
"Low"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(LOW_);}
"LOCA"|"Loca"|"loca" 			{print_to_log(); strcpy(yylval.string, yytext);   return(LOCA_);}
"matrix"|"Matrix"|"MATRIX"		{print_to_log(); strcpy(yylval.string, yytext);   return(MATRIX_);} 
"MAT"|"mat"				{print_to_log(); strcpy(yylval.string, yytext);   return(MAT_);}
"MOM"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOM_);}
"MOVING_CA"				{print_to_log(); strcpy(yylval.string, yytext);   return(MOVING_CA_);}
"MOVING_PLANE"|"MOVING_PLANE_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(MOVING_PLANE_);}
"MU" 					{print_to_log(); strcpy(yylval.string, yytext);   return(MU_);}
"Mapping"  				{print_to_log(); strcpy(yylval.string, yytext);   return(MAPPING_);} 
"Mass"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(MASS_);}  
"Materials"				{print_to_log(); strcpy(yylval.string, yytext);   return(MATERIALS_);}
"Maximum"|"maximum"|"MAXIMUM"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(MAXIMUM_);} 
"Mean" 					{print_to_log(); strcpy(yylval.string, yytext);   return(MEAN_);}  
"Mesh" 					{print_to_log(); strcpy(yylval.string, yytext);   return(MESH_);} 
"Minimum"|"minimum"|"MINIMUM"		{print_to_log(); strcpy(yylval.string, yytext);   return(MINIMUM_);}  
"Mobility"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOBILITY_);}
"Model"					{print_to_log(); strcpy(yylval.string, yytext);   return(MODEL_);}
"Modified"|"modified"|"MODIFIED"	{print_to_log(); strcpy(yylval.string, yytext);   return(MODIFIED_);}
"Motion"  				{print_to_log(); strcpy(yylval.string, yytext);   return(MOTION_);} 
"Moving"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(MOVING_);} 
"N" 					{print_to_log(); strcpy(yylval.string, yytext);   return(N_);}
"N2" 					{print_to_log(); strcpy(yylval.string, yytext);   return(N2_);}
"N3" 					{print_to_log(); strcpy(yylval.string, yytext);   return(N3_);}
"NAME"|"Name"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(NAME_);}
"NC"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(NC_);}
"NEWTONIAN"				{print_to_log(); strcpy(yylval.string, yytext);   return(NEWTONIAN_);}
"NO"|"no"|"No"				{print_to_log(); strcpy(yylval.string, yytext);   return(NO_);}
"NONE"|"none"|"None"			{print_to_log(); strcpy(yylval.string, yytext);   return(NONE_);}
"NONLINEAR"				{print_to_log(); strcpy(yylval.string, yytext);   return(NONLINEAR_);}
"NOPOLYMER"				{print_to_log(); strcpy(yylval.string, yytext);   return(NOPOLYMER_);}
"NORM_FLOCE_RS"|"NORM_FORCE_RS_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(NORM_FORCE_RS_);} 
"NORM_FORCE"|"NORM_FORCE_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(NORM_FORCE_);} 
"NO_SLIP"|"NO_SLIP_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(NO_SLIP_);}
"NO_SLIP_RS"|"NO_SLIP_RS_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(NO_SLIP_RS_);}
"NS"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(NS_);}
"Navier" 				{print_to_log(); strcpy(yylval.string, yytext);   return(NAVIER_);}
"Network" 				{print_to_log(); strcpy(yylval.string, yytext);   return(NETWORK_);}
"Neumann"|"NEUMANN"|"neumann"		{print_to_log(); strcpy(yylval.string, yytext);   return(NEUMANN_);}
"Newton"|"NEWTON"|"newton" 		{print_to_log(); strcpy(yylval.string, yytext);   return(NEWTON_);} 
"Normalized"|"normalized"|"NORMALIZED" 	{print_to_log(); strcpy(yylval.string, yytext);   return(NORMALIZED_);}  
"Number"|"number"|"Number"		{print_to_log(); strcpy(yylval.string, yytext);   return(NUMBER_);}
"HUNTING"|"Hunting"|"hunting"		{print_to_log(); strcpy(yylval.string, yytext);   return(HUNTING_);}
"OF"|"of"|"Of"				{print_to_log(); strcpy(yylval.string, yytext);   return(OF_);}
"RECORD"|"Record"|"record"		{print_to_log(); strcpy(yylval.string, yytext);   return(RECORD_);}
"MODES"|"Modes"|"modes"			{print_to_log(); strcpy(yylval.string, yytext);   return(MODES_);}
"OLDROYDB"				{print_to_log(); strcpy(yylval.string, yytext);   return(OLDROYDB_);}
"Orthogonalization"|"ORTHOGONALIZATION"|"orthogonalization"	{print_to_log(); strcpy(yylval.string, yytext);   return(ORTHOGONALIZATION_);}  
"P" 	 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(P_);}	 
"P"|"P_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(P_);}
"P0"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(P0_);}
"PATH"|"Path"|"path"			{print_to_log(); strcpy(yylval.string, yytext);   return(PATH_);}
"P1"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(P1_);}
"pf"|"PF" 				{print_to_log(); strcpy(yylval.string, yytext);   return(PF_);}
"PARTIALLY_WETTING"			{print_to_log(); strcpy(yylval.string, yytext);   return(PARTIALLY_WETTING_);}
"PENETRATION"|"PENETRATION_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(PENETRATION_);}
"PLANE"|"PLANE_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(PLANE_);}
"PLANEX"|"PLANEX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(PLANEX_);}  
"PLANEY"|"PLANEY_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(PLANEY_);} 
"PLANEZ"|"PLANEZ_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(PLANEZ_);}
"POROUS_BRINKMAN"			{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_BRINKMAN_);}
"POROUS_CONV"|"POROUS_CONV_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_CONV_);}
"POROUS_FLUX"|"POROUS_FLUX_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_FLUX_);}
"POROUS_GAS"|"POROUS_GAS_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_GAS_);}
"POROUS_KIN"|"POROUS_KIN_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_KIN_);}
"POROUS_PRESSURE"|"POROUS_PRESSURE_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_PRESSURE_);}
"POROUS_SATURATED"			{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_SATURATED_);}
"POROUS_TWO_PHASE"			{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_TWO_PHASE_);}
"POROUS_UNSATURATED"			{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_UNSATURATED_);}
"POWER_LAW"				{print_to_log(); strcpy(yylval.string, yytext);   return(POWER_LAW_);}
"PQ1"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(PQ1_);}
"PQ2"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(PQ2_);}
"PRESSURE_USER"|"PRESSURE_USER_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(PRESSURE_USER_);}
"PSD_SEXP"				{print_to_log(); strcpy(yylval.string, yytext);   return(PSD_SEXP_);}
"PSD_VOL"				{print_to_log(); strcpy(yylval.string, yytext);   return(PSD_VOL_);}
"PSD_WEXP"				{print_to_log(); strcpy(yylval.string, yytext);   return(PSD_WEXP_);}
"PSPG"|"PSPG_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(PSPG_);} 
"PTT"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(PTT_);}
"PU"|"PU_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(PU_);}
"PV"|"PV_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(PV_);}
"PVX" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(PVX_);}	 
"PVY" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(PVY_);}	 
"PVZ" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(PVZ_);}	 
"PW"|"PW_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(PW_);}
"P_EQUIL"|"P_EQUIL_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(P_EQUIL_);}
"Packing"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(PACKING_);}
"Permeability"|"permeability"|"PERMEABILITY"	{print_to_log(); strcpy(yylval.string, yytext);   return(PERMEABILITY_);}
"Plastic"				{print_to_log(); strcpy(yylval.string, yytext);   return(PLASTIC_);}
"Point"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(POINT_);}
"Polymer"				{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_);}
"Porosity" 				{print_to_log(); strcpy(yylval.string, yytext);   return(POROSITY_);}
"Post"|"POST"|"post"			{print_to_log(); strcpy(yylval.string, yytext);   return(POST_);}
"Power"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(POWER_);}
"Preconditioner"|"PRECONDITIONER"|"preconditioner"	{print_to_log(); strcpy(yylval.string, yytext);   return(PRECONDITIONER_);}  
"Preconditioning" 			{print_to_log(); strcpy(yylval.string, yytext);   return(PRECONDITIONS_);}  
"Printing"|"PRINTING"|"printing"	{print_to_log(); strcpy(yylval.string, yytext);   return(PRINTING_);} 
"Problem"				{print_to_log(); strcpy(yylval.string, yytext);   return(PROBLEM_);}
"Processing"|"PROCESSING"|"processing"		{print_to_log(); strcpy(yylval.string, yytext);   return(PROCESSING_);}
"Pseudo-Solid"				{print_to_log(); strcpy(yylval.string, yytext);   return(PSEUDO_SOLID_);}
"Q1"					{print_to_log(); strcpy(yylval.string, yytext);   return(Q1_);}
"Q1_D"					{print_to_log(); strcpy(yylval.string, yytext);   return(Q1_D_);}
"Q2"					{print_to_log(); strcpy(yylval.string, yytext);   return(Q2_);}
"Q2_D"					{print_to_log(); strcpy(yylval.string, yytext);   return(Q2_D_);}
"Q2_D_LSA"				{print_to_log(); strcpy(yylval.string, yytext);   return(Q2_D_LSA_);}
"Q2_LSA"				{print_to_log(); strcpy(yylval.string, yytext);   return(Q2_LSA_);}
"QCONV"|"QCONV_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(QCONV_);} 
"QRAD"|"QRAD_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(QRAD_);}
"QSIDE"|"QSIDE_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(QSIDE_);} 
"QUAD_GP"   				{print_to_log(); strcpy(yylval.string, yytext);   return(QUAD_GP_);}  
"QUSER"|"QUSER_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(QUSER_);}
"Q_VELO_SLIP"|"Q_VELO_SLIP_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(Q_VELO_SLIP_);} 
"RAOULT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(RAOULT_);}
"REP_FORCE"|"REP_FORCE_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(REP_FORCE_);}
"REP_FORCE_RS"|"REP_FORCE_RS_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(REP_FORCE_RS_);}
"R_ENERGY" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_ENERGY_);}	 
"R_FILL" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_FILL_);}	 
"R_GRADIENT11" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_GRADIENT11_);}	 
"R_GRADIENT12" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_GRADIENT12_);}	 
"R_GRADIENT13" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_GRADIENT13_);}	 
"R_GRADIENT21" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_GRADIENT21_);}	 
"R_GRADIENT22" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_GRADIENT22_);}	 
"R_GRADIENT23" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_GRADIENT23_);}	 
"R_GRADIENT31" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_GRADIENT31_);}	 
"R_GRADIENT32" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_GRADIENT32_);}	 
"R_GRADIENT33" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_GRADIENT33_);}	 
"R_MASS"  	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MASS_);}	 
"R_MASS_SURF"  	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MASS_SURF_);}	 
"R_MESH1"  	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MESH1_);}	 
"R_MESH2"   	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MESH2_);}	 
"R_MESH3"  	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MESH3_);}	 
"R_MESH_NORMAL"  			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MESH_NORMAL_);}	 
"R_MESH_TANG1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MESH_TANG1_);}	 
"R_MESH_TANG2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MESH_TANG2_);}	 
"R_MOMENTUM1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MOMENTUM1_);}	 
"R_MOMENTUM2"  	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MOMENTUM2_);}	 
"R_MOMENTUM3"  	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MOMENTUM3_);}	 
"R_MOM_NORMAL" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MOM_NORMAL_);}	 
"R_MOM_TANG1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MOM_TANG1_);}	 
"R_MOM_TANG2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MOM_TANG2_);}	 
"R_PMOMENTUM1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MOMENTUM1_);}	 
"R_PMOMENTUM2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MOMENTUM2_);}	 
"R_PMOMENTUM3" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_MOMENTUM3_);}	 
"R_POTENTIAL" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_POTENTIAL_);}	 
"R_PRESSURE"   	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_PRESSURE_);}	 
"R_SHEAR_RATE" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SHEAR_RATE_);}	 
"R_SOLID1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SOLID1_);}	 
"R_SOLID2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SOLID2_);}	 
"R_SOLID3" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SOLID3_);}	 
"R_STRESS11"   	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS11_);}	 
"R_STRESS11_1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS11_1_);}	 
"R_STRESS11_2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS11_2_);}	 
"R_STRESS11_3" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS11_3_);}	 
"R_STRESS11_4" 				{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS11_4_);}	 
"R_STRESS11_5" 				{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS11_5_);}	 
"R_STRESS11_6" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS11_6_);}	 
"R_STRESS11_7" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS11_7_);}	 
"R_STRESS12"|"R_STRESS12_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS12_);}	 
"R_STRESS12_1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS12_1_);}	 
"R_STRESS12_2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS12_2_);}	 
"R_STRESS12_3" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS12_3_);}	 
"R_STRESS12_4" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS12_4_);}	 
"R_STRESS12_5" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS12_5_);}	 
"R_STRESS12_6" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS12_6_);}	 
"R_STRESS12_7" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS12_7_);}	 
"R_STRESS13" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS13_);}	 
"R_STRESS13_1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS13_1_);}	 
"R_STRESS13_2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS13_2_);}	 
"R_STRESS13_3" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS13_3_);}	 
"R_STRESS13_4" 				{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS13_4_);}	 
"R_STRESS13_5" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS13_5_);}	 
"R_STRESS13_6" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS13_6_);}	 
"R_STRESS13_7" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS13_7_);}	 
"R_STRESS22" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS22_);}	 
"R_STRESS22_1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS22_1_);}	 
"R_STRESS22_2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS22_2_);}	 
"R_STRESS22_3" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS22_3_);}	 
"R_STRESS22_4" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS22_4_);}	 
"R_STRESS22_5" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS22_5_);}	 
"R_STRESS22_6" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS22_6_);}	 
"R_STRESS22_7" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS22_7_);}	 
"R_STRESS23"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS23_);}	 
"R_STRESS23_1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS23_1_);}	 
"R_STRESS23_2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS23_2_);}	 
"R_STRESS23_3" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS23_3_);}	 
"R_STRESS23_4" 				{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS23_4_);}	 
"R_STRESS23_5" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS23_5_);}	 
"R_STRESS23_6" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS23_6_);}	 
"R_STRESS23_7" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS23_7_);}	 
"R_STRESS33" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS33_);}	 
"R_STRESS33_1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS33_1_);}	 
"R_STRESS33_2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS33_2_);}	 
"R_STRESS33_3" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS33_3_);}	 
"R_STRESS33_4" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS33_4_);}	 
"R_STRESS33_5" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS33_5_);}	 
"R_STRESS33_6" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS33_6_);}	 
"R_STRESS33_7" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_STRESS33_7_);}	 
"R_Species_0" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_0_);}	 
"R_Species_1" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_1_);}	 
"R_Species_10" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_10_);}	 
"R_Species_11" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_11_);}	 
"R_Species_12" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_12_);}	 
"R_Species_13" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_13_);}	 
"R_Species_14" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_14_);}	 
"R_Species_15" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_15_);}	 
"R_Species_16" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_16_);}	 
"R_Species_17" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_17_);}	 
"R_Species_18" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_18_);}	 
"R_Species_19" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_19_);}	 
"R_Species_2" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_2_);}	 
"R_Species_20" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_20_);}	 
"R_Species_21" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_21_);}	 
"R_Species_22" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_22_);}	 
"R_Species_23" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_23_);}	 
"R_Species_24" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_24_);}	 
"R_Species_25" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_25_);}	 
"R_Species_26" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_26_);}	 
"R_Species_27" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_27_);}	 
"R_Species_28" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_28_);}	 
"R_Species_29" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_29_);}	 
"R_Species_3" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_3_);}	 
"R_Species_4" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_4_);}	 
"R_Species_5" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_5_);}	 
"R_Species_6" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_6_);}	 
"R_Species_7" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_7_);}	 
"R_Species_8" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_8_);}	 
"R_Species_9" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_SPECIES_9_);}	 
"R_VolFracPh_0"  			{print_to_log(); strcpy(yylval.string, yytext);   return(R_VOLFRACPH_0_);}	 
"R_VolFracPh_1"  			{print_to_log(); strcpy(yylval.string, yytext);   return(R_VOLFRACPH_1_);}	 
"R_VolFracPh_2"  			{print_to_log(); strcpy(yylval.string, yytext);   return(R_VOLFRACPH_2_);}	 
"R_VolFracPh_3"  			{print_to_log(); strcpy(yylval.string, yytext);   return(R_VOLFRACPH_3_);}	 
"R_VolFracPh_4"  			{print_to_log(); strcpy(yylval.string, yytext);   return(R_VOLFRACPH_4_);}	 
"R_Y0" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y0_);}	 
"R_Y1" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y1_);}	 
"R_Y10" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y10_);}	 
"R_Y11" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y11_);}	 
"R_Y12" 				{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y12_);}	 
"R_Y13" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y13_);}	 
"R_Y14" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y14_);}	 
"R_Y15" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y15_);}	 
"R_Y16" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y16_);}	 
"R_Y17" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y17_);}	 
"R_Y18" 				{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y18_);}	 
"R_Y19" 				{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y19_);}	 
"R_Y2" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y2_);}	 
"R_Y20"  	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y20_);}	 		 
"R_Y21" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y21_);}	 
"R_Y22" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y22_);}	 
"R_Y23" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y23_);}	 
"R_Y24" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y24_);}	 
"R_Y25" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y25_);}	 
"R_Y26" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y26_);}	 
"R_Y27" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y27_);}	 
"R_Y28" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y28_);}	 
"R_Y29" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y29_);}	 
"R_Y3" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y3_);}	 
"R_Y4" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y4_);}	 
"R_Y5" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y5_);}	 
"R_Y6" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y6_);}	 
"R_Y7" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y7_);}	 
"R_Y8" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y8_);}	 
"R_Y9" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(R_Y9_);}	 
"Rate"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(RATE_);}
"Ratio"|"RATIO"|"ratio"	 		{print_to_log(); strcpy(yylval.string, yytext);   return(RATIO_);}  
"Real"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(REAL_);} 
"Reference"|"REFERENCE"|"reference"	{print_to_log(); strcpy(yylval.string, yytext);   return(REFERENCE_);}
"Reform"|"REFORM"|"reform" 		{print_to_log(); strcpy(yylval.string, yytext);   return(REFORM_);}
"Rel "		 			{print_to_log(); strcpy(yylval.string, yytext);   return(REL_);}
"Residuals" 				{print_to_log(); strcpy(yylval.string, yytext);   return(RESIDUALS_);}  
"S"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(S_);}
"S11"|"S11_BC"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(S11_);}	 
"S11_1"|"S11_1_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S11_1_);}
"S11_2"|"S11_2_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S11_2_);}
"S11_3"|"S11_3_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S11_3_);}
"S11_4"|"S11_4_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S11_4_);}
"S11_5"|"S11_5_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S11_5_);}
"S11_6"|"S11_6_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S11_6_);}
"S11_7"|"S11_7_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S11_7_);}
"S12"|"S12_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(S12_);}
"S12_1"|"S12_1_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S12_1_);}
"S12_2"|"S12_2_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S12_2_);}
"S12_3"|"S12_3_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S12_3_);}
"S12_4"|"S12_4_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S12_4_);}	 
"S12_5"|"S12_5_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S12_5_);}
"S12_6"|"S12_6_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S12_6_);}	 
"S12_7"|"S12_7_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S12_7_);}
"S13"|"S13_BC" 				{print_to_log(); strcpy(yylval.string, yytext);   return(S13_);}	 
"S13_1"|"S13_1_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S13_1_);}
"S13_2"|"S13_2_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S13_2_);}	 
"S13_3"|"S13_3_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S13_3_);}
"S13_4"|"S13_4_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S13_4_);}	 
"S13_5"|"S13_5_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S13_5_);}
"S13_6"|"S13_6_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S13_6_);}	 
"S13_7"|"S13_7_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S13_7_);}
"S21"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(S21_);}
"S22"|"S22_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(S22_);}
"S22_1"|"S22_1_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S22_1_);}	 
"S22_2"|"S22_2_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S22_2_);}
"S22_3"|"S22_3_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S22_3_);}	 
"S22_4"|"S22_4_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S22_4_);}
"S22_5"|"S22_5_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S22_5_);}	 
"S22_6"|"S22_6_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S22_6_);}
"S22_7"|"S22_7_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S22_7_);}	 
"S23"|"S23_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(S23_);}
"S23_1"|"S23_1_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S23_1_);}
"S23_2"|"23_2_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S33_2_);}
"S23_3"|"S23_3_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S23_3_);}	 
"S23_4"|"S23_4_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S23_4_);}
"S23_5"|"S23_5_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S23_5_);}	 
"S23_6"|"S23_6_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(S23_6_);}
"S23_7"|"S23_7_BC" 	 		{print_to_log(); strcpy(yylval.string, yytext);   return(S23_7_);}	 
"S33"|"S33_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(S33_);}
"S33_1"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(S33_1_);}	 
"S33_2"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(S33_2_);}	 
"S33_3"					{print_to_log(); strcpy(yylval.string, yytext);   return(S33_3_);}
"S33_4" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(S33_4_);}	 
"S33_5"					{print_to_log(); strcpy(yylval.string, yytext);   return(S33_5_);}
"S33_6"					{print_to_log(); strcpy(yylval.string, yytext);   return(S33_6_);}	 
"S33_7"					{print_to_log(); strcpy(yylval.string, yytext);   return(S33_7_);}
"SC"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(SC_);}
"SEED" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SEED_);}
"SH"|"SH_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(SH_);}
"SINUSOIDAL"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(SINUSOIDAL_);} 
"SLOPE"|"SLOPE_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(SLOPE_);} 
"SLOPEX"|"SLOPEX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(SLOPEX_);} 
"SLOPEY"|"SLOPEY_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(SLOPEY_);} 
"SLOPEZ"|"SLOPEZ_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(SLOPEZ_);} 
"SOLID_FLUID"|"SOLID_FLUID_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_FLUID_);}
"SOLID_FLUID_RS"|"SOLID_FLUID_RS_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_FLUID_RS_);}
"SOLN"					{print_to_log(); strcpy(yylval.string, yytext);   return(SOLN_);}
"SP"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(SP_);}
"SPHERICAL" 				{print_to_log(); strcpy(yylval.string, yytext);   return(SPHERICAL_);}
"FRACTION"|"Fraction"|"fraction"	{print_to_log(); strcpy(yylval.string, yytext);   return(FRACTION_);}
"SPLINE"|"SPLINE_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPLINE_);}
"SPLINEX"|"SPLINEX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPLINEX_);}  
"SPLINEY"|"SPLINEY_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPLINEY_);} 
"SPLINZ"|"SPLINEZ_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPLINEZ_);}
"SS"|"ss"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SS_);}
"STRONG_FILL_CA"|"STRONG_FILL_CA_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(STRONG_FILL_CA_);}
"STRONG_INT_EDGE" 			{print_to_log(); strcpy(yylval.string, yytext);   return(STRONG_INT_EDGE_);}
"SUM_TO_ONE"				{print_to_log(); strcpy(yylval.string, yytext);   return(SUM_TO_ONE_);}
"SUPG"					{print_to_log(); strcpy(yylval.string, yytext);   return(SUPG_);}
"SURFDOMAIN_MASS_FRACTION"|"SURFDOMAIN_MASS_FRACTION_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(SURFCOMAIN_MASS_FRACTION_);}
"SURFDOMAIN_ENERGY_BAL"|"SURFDOMAIN_ENERGY_BAL_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(SURFCOMAIN_ENERGY_BAL_FRACTION_);}
"SURFDOMAIN_STEFAN_FLOW"|"SURFDOMAIN_STEFAN_FLOW_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(SURFCOMAIN_STEFAN_FLOW_);}
"SURFTANG"|"SURFTANG_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(SURFTANG_);}
"SURFTANG_EDGE"|"SURFTANG_EDGE_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(SURFTANG_EDGE_);}
"SURFTANG_SCALAR"|"SURFTANG_SCALAR_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(SURFTANG_SCALAR_);}
"SURFTANG_SCALAR_EDGE"|"SURFTANG_SCALAR_EDGE_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(SURFTANG_SCALAR_EDGE_);}
"SUSPENSION"|"Suspension"|"suspension"	{print_to_log(); strcpy(yylval.string, yytext);   return(SUSPENSION_);}
"SUSPENSION_PM"|"Suspension_pm"|"suspension_pm"	{print_to_log(); strcpy(yylval.string, yytext);   return(SUSPENSION_PM_);}
"SWIRLING"				{print_to_log(); strcpy(yylval.string, yytext);   return(SWIRLING_);}
"Saturation" 				{print_to_log(); strcpy(yylval.string, yytext);   return(SATURATION_);}
"Scaling"|"scaling"|"SCALING"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(SCALING_);} 
"Second"|"SECOND"|"second"		{print_to_log(); strcpy(yylval.string, yytext);   return(SECOND_);}  
"Sensitivities"|"sensitivities"|"SENSITIVITIES"		{print_to_log(); strcpy(yylval.string, yytext);   return(SENSITIVITIES_);}
"Shear"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(SHEAR_);} 
"Size"|"size"|"Size"	 		{print_to_log(); strcpy(yylval.string, yytext);   return(SIZE_);} 
"Solid"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_);}
"Solidus"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLIDUS_);}
"Solution"|"solution"|"SOLUTION"	{print_to_log(); strcpy(yylval.string, yytext);   return(SOLUTION_);}
"Solve"|"SOLVE"|"solve"			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLVE_);}
"EIGEN"|"Eigen"|"eigen"			{print_to_log(); strcpy(yylval.string, yytext);   return(EIGEN_);}
"EIGENSOLVER"|"Eigensolver"|"eigensolver"	{print_to_log(); strcpy(yylval.string, yytext);   return(EIGENSOLVER_);}
"Solvent"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(SOLVENT_);}
"Solver"|"SOLVER"|"solver"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(SOLVER_);} 
"SPHERE" 				{print_to_log(); strcpy(yylval.string, yytext);   return(SPHERE_);}
"Sp_0" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_0_);}
"Sp_0" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_0_);}	 
"Sp_1" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_1_);}	 
"Sp_10" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_10_);}	 
"Sp_11" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_11_);}	 
"Sp_12" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_12_);}	 
"Sp_13" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_13_);}	 
"Sp_14" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_14_);}	 
"Sp_15" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_15_);}	 
"Sp_16" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_16_);}	 
"Sp_17" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_17_);}	 
"Sp_18" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_18_);}	 
"Sp_19" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_19_);}	 
"Sp_2" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_2_);}	 
"Sp_20" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_20_);}	 
"Sp_21" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_21_);}	 
"Sp_22" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_22_);}	 
"Sp_23"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_23_);}	 
"Sp_24" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_24_);}	 
"Sp_25" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_25_);}	 
"Sp_26" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_26_);}	 
"Sp_27" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_27_);}	 
"Sp_28" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_28_);}	 
"Sp_29" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SP_29_);}	 
"Sp_3" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_3_);}	 
"Sp_4" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_4_);}	 
"Sp_5" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_5_);}	 
"Sp_6" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_6_);}	 
"Sp_7" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_7_);}	 
"Sp_8" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_8_);}	 
"Sp_9" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SP_9_);}
"PROPERTY"|"Property"|"property"	{print_to_log(); strcpy(yylval.string, yytext);   return(PROPERTY_);}	 
"Specifications"|"Specifications"|"specifications"	{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIFICATIONS_);}
"Stabilization"|"stabilization"|"STABILIZATION"		{print_to_log(); strcpy(yylval.string, yytext);   return(STABILIZATION_);} 
"Stability"|"stability"|"STABILITY"	{print_to_log(); strcpy(yylval.string, yytext);   return(STABILITY_);} 
"Steps"|"steps"|"STEPS"	 		{print_to_log(); strcpy(yylval.string, yytext);   return(STEPS_);} 
"Stokes"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(STOKES_);} 
"Strain" 				{print_to_log(); strcpy(yylval.string, yytext);   return(STRAIN_);} 
"Stream"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(STREAM_);}  
"Streamwise"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(STREAMWISE_);}  
"Stress"|"stress"			{print_to_log(); strcpy(yylval.string, yytext);   return(STRESS_);}
"Stride"|"STRIDE"|"stride"		{print_to_log(); strcpy(yylval.string, yytext);   return(STRIDE_);} 
"RECYCLE"|"Recycle"|"recycle"		{print_to_log(); strcpy(yylval.string, yytext);   return(RECYCLE_);} 
"Subcycle"|"SUBCYCLE"|"subcycle"	{print_to_log(); strcpy(yylval.string, yytext);   return(SUBCYCLE_);}\
"Subdomain"|"SUBDOMAIN"|"subdomain"	{print_to_log(); strcpy(yylval.string, yytext);   return(SUBDOMAIN_);}
"Surface"|"SURFACE"|"surface"		{print_to_log(); strcpy(yylval.string, yytext);   return(SURFACE_);}
"SURF"|"Surf"|"surf"			{print_to_log(); strcpy(yylval.string, yytext);   return(SURF_);}
"System"		  		{print_to_log(); strcpy(yylval.string, yytext);   return(SYSTEM_);} 
"SUBINDEX"|"Subindex"|"subindex"	{print_to_log(); strcpy(yylval.string, yytext);   return(SUBINDEX_);} 
"T"|"T_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(T_);} 
"T1" 					{print_to_log(); strcpy(yylval.string, yytext);   return(T1_);}
"T2" 					{print_to_log(); strcpy(yylval.string, yytext);   return(T2_);}
"IDEAL_GAS"|"Ideal_Gas"|"ideal_gas"	{print_to_log(); strcpy(yylval.string, yytext);   return(IDEAL_GAS_);}
"TABLE"|"TABLE_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(TABLE_);}
"TABLE_WICS"|"TABLE_WICS_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(TABLE_WICS_);}
"TABLE_WICV"|"TABLE_WICV_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(TABLE_WICV_);}
"TENSOR"				{print_to_log(); strcpy(yylval.string, yytext);   return(TENSOR_);}
"THERMAL"|"Thermal"|"thermal"		{print_to_log(); strcpy(yylval.string, yytext);   return(THERMAL_);}
"THERMAL_BATTERY"|"Thermal_battery"|"thermal_battery"  {print_to_log(); strcpy(yylval.string, yytext);   return(THERMAL_BATTERY_);}
"TNRMLSIDE"				{print_to_log(); strcpy(yylval.string, yytext);   return(TNRMLSIDE_);} 
"TOTAL_ALE"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(TOTAL_ALE_);}
"TSHRSIDE"				{print_to_log(); strcpy(yylval.string, yytext);   return(TSHRSIDE_);} 
"T_MELT"|"T_MELT_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(T_MELT_);}
"T_USER"|"T_USER_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(T_USER_);}
"Tangent"|"TANGENT"|"tangent"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(TANGENT_);} 
"Temperature"|"TEMPERATURE"|"temperature"  {print_to_log(); strcpy(yylval.string, yytext);   return(TEMPERATURE_);}
"Tension"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(TENSION_);}
"Tensor" 				{print_to_log(); strcpy(yylval.string, yytext);   return(TENSOR_);} 
"FLOAT"|"Float"|"float"			{print_to_log(); strcpy(yylval.string, yytext);   return(THE_WORD_FLOAT_);}
"Third"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(THIRD_);}  
"Time"|"TIME"|"time" 			{print_to_log(); strcpy(yylval.string, yytext);   return(TIME_);}
"U"|"U_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(U_);} 
"U1"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(U1_);}
"U2"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(U2_);}
"U3"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(U3_);}
"UMF"|"umf"				{print_to_log(); strcpy(yylval.string, yytext);   return(UMF_);}
"umff"|"UMFF"				{print_to_log(); strcpy(yylval.string, yytext);   return(UMFF_);}
"UNSEEDED"				{print_to_log(); strcpy(yylval.string, yytext);   return(UNSEEDED_);}
"USER"					{print_to_log(); strcpy(yylval.string, yytext);   return(USER_);}
"USER_GEN"				{print_to_log(); strcpy(yylval.string, yytext);   return(USER_GEN_);}
"UUSER"|"UUSER_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(UUSER_);}
"UVARY"|"UVARY_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(UVARY_);} 
"Unreacted"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(UNREACTED_);}
"User-Defined"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(USER_DEFINED_);} 
"V"|"V_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(V_);} 
"VAN_GENUCHTEN"				{print_to_log(); strcpy(yylval.string, yytext);   return(VAN_GENUCHTEN_);}
"VAR_CA_EDGE"|"VAR_CA_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VAR_CA_EDGE_);}
"VAR_CA_USER"|"VAR_CA_USER_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VAR_CA_USER_);}
"VELO_NORMAL"|"VELO_NORMAL_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_NORMAL_);}
"VELO_NORMAL_DISC"|"VELO_NORMAL_DISC_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_NORMAL_DISC_);}
"VELO_NORMAL_EDGE"|"VELO_NORMAL_EDGE_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_NORMAL_EDGE_);}
"VELO_NORMAL_EDGE_INT"|"VELO_NORMAL_EDGE_INT_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_NORMAL_EDGE_INT_);}
"VELO_NORM_COLLOC"|"VELO_NORM_COLLOC_BC"  		{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_NORM_COLLOC_);}
"VELO_SLIP_FILL"|"VELO_SLIP_FILL_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_SLIP_FILL_);}
"VELO_SLIP_ROT"|"VELO_SLIP_ROT_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_SLIP_ROT_);}
"VELO_SLIP_SOLID"|"VELO_SLIP_SOLID_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_SLIP_SOLID_);}
"VELO_SLIP"|"VELO_SLIP_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_SLIP_);}
"VELO_TANGENT"|"VELO_TANGENT_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_TANGENT_);}
"VELO_TANGENT_3D"|"VELO_TANGENT_3D_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_TANGENT_3D_);}
"VELO_TANGENT_EDGE"|"VELO_TANGENT_EDGE_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_TANGENT_EDGE_);}
"VELO_TANGENT_EDGE_INT"|"VELO_TANGENT_EDGE_INT_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_TANGENT_EDGE_INT_);}
"VELO_TANGENT_SOLID"|"VELO_TANGENT_SOLID_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VELO_TANGENT_SOLID_);}
"VERTEX"				{print_to_log(); strcpy(yylval.string, yytext);   return(VERTEX_);}
"VFP_0" 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(VFP_0_);}	 
"VFP_1" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VFP_1_);}	 
"VFP_2" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VFP_2_);}	 
"VFP_3" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VFP_3_);}	 
"VFP_4" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VFP_4_);}	 
"VL_EQUIL"|"VL_EQUIL_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VL_EQUIL_);}
"VL_POLY"|"VL_POLY_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(VL_POLY_);}
"VN" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(VN_);}	 
"VNORM_LEAK"|"VNORM_LEAK_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VNORM_LEAK_);} 
"VN_POROUS"|"VN_POROUS_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VN_POROUS_);}
"VOLT"|"VOLT_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(VOLT_);}
"VP_EQUIL"|"VP_EQUIL_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(VP_EQUIL_);}
"VT1" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(VT1_);}	 
"VT2" 	 				{print_to_log(); strcpy(yylval.string, yytext);   return(VT2_);}
"VUSER"|"VUSER_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(VUSER_);}
"VVARY"|"VVARY_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(VVARY_);} 
"VX" 	 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(VX_);}	 
"VY" 	 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(VY_);}	 
"VZ" 	 	 			{print_to_log(); strcpy(yylval.string, yytext);   return(VZ_);}	 
"Vapor"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(VAPOR_);}
"Vaporization" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VAPORIZATION_);}
"Vector"|"VECTOR"|"vector"		{print_to_log(); strcpy(yylval.string, yytext);   return(VECTOR_);}
"Velocity"|"VELOCIYT"|"velocity"	{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY_);}
"Viscosity"|"viscosity"|"VISCOSITY"	{print_to_log(); strcpy(yylval.string, yytext);   return(VISCOSITY_);}
"SHIFTS"|"Shifts"|"shifts"		{print_to_log(); strcpy(yylval.string, yytext);   return(SHIFTS_);}
"WAVE"|"Wave"|"wave"			{print_to_log(); strcpy(yylval.string, yytext);   return(WAVE_);}
"NUMBERS"|"Numbers"|"numbers"		{print_to_log(); strcpy(yylval.string, yytext);   return(NUMBERS_);}			
"Vol"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(VOL_);}
"Volume"|"VOLUME"|"volume"		{print_to_log(); strcpy(yylval.string, yytext);   return(VOLUME_);}
"W"|"W_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(W_);}
"WEAK_INT_EDGE"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(WEAK_INT_EDGE_);}
"WETTING"				{print_to_log(); strcpy(yylval.string, yytext);   return(WETTING_);}
"WHITE-METZNER"				{print_to_log(); strcpy(yylval.string, yytext);   return(WHITE_METZNER_);}
"WLF"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(WLF_);}
"WUSER"|"WUSER_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(WUSER_);}
"WVARY"|"WVARY_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(WVARY_);} 
"Weight"|"WEIGHT"|"weight"		{print_to_log(); strcpy(yylval.string, yytext);   return(WEIGHT_);}
"Weighting"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(WEIGHTING_);}
"Write"					{print_to_log(); strcpy(yylval.string, yytext);   return(WRITE_);}
"X" 					{print_to_log(); strcpy(yylval.string, yytext);   return(X_);}
"XDIM"|"Xdim"|"xdim"			{print_to_log(); strcpy(yylval.string, yytext);   return(XDIM_);}
"Xi"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(XI_);}
"Y"|"Y_BC"				{print_to_log(); strcpy(yylval.string, yytext);   return(Y_);}
"FORCE_X"|"force_x"			{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_X_);}	
"FORCE_Y"|"force_y"			{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_Y_);}
"FORCE_Z"|"force_z"			{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_Z_);}
"FORCE_NORMAL"|"force_normal"		{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_NORMAL_);}
"FORCE_TANGENT1"|"force_tangent1"	{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_TANGENT1_);}
"FORCE_TANGENT2"|"force_tangent2"	{print_to_log(); strcpy(yylval.string, yytext);   return(FORCE_TANGENT2_);}
"VOLUME_FLUX"|"volume_flux"		{print_to_log(); strcpy(yylval.string, yytext);   return(VOLUME_FLUX_);}
"HEAT_FLUX"|"heat_flux"			{print_to_log(); strcpy(yylval.string, yytext);   return(HEAT_FLUX_);}
"SPECIES_FLUX"|"species_flux"		{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_FLUX_);}
"Y0" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y0_);}	 
"Y1" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y1_);}	 
"Y10" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y10_);}	 
"Y11" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y11_);}	 
"Y12" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y12_);}	 
"Y13" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y13_);}	 
"Y14" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y14_);}	 
"Y15" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y15_);}	 
"Y16" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y16_);}	 
"Y17" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y17_);}	 
"Y18" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y18_);}	 
"Y19" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y19_);}	 
"Y2" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y2_);}	 
"Y20" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y20_);}	 
"Y21" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y21_);}	 
"Y22" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y22_);}	 
"Y23" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y23_);}	 
"Y24" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y24_);}	 
"Y25" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y25_);}	  
"Y26" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y26_);}	 
"Y27" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y27_);}	 
"Y28" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y28_);}	 
"Y29" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y29_);}	 
"Y2_ELECTRONEURALITY"|"Y2_ELECTRONEUTRALITY_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(Y2_ELECTRONEUTRALITY_);}   
"Y3" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y3_);}	 
"Y4" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y4_);}	 
"Y5" 		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y5_);}	 
"Y6" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y6_);}	 
"Y7" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y7_);}	 
"Y8" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y8_);}	 
"Y9" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Y9_);}	 
"YES"|"yes"|"Yes"			{print_to_log(); strcpy(yylval.string, yytext);   return (YES_);}
"FLUXES"|"fluxes"|"Fluxes"		{print_to_log(); strcpy(yylval.string, yytext);   return (FLUXES_);}
"FLUX"|"flux"|"Flux"			{print_to_log(); strcpy(yylval.string, yytext);   return (FLUX_);}
"FLUX_SENS"|"flux_sens"|"Flux_Sens"	{print_to_log(); strcpy(yylval.string, yytext);   return (FLUX_SENS_);}
"VC"|"vc"|"Vc"				{print_to_log(); strcpy(yylval.string, yytext);   return (VC_);}
"FC"|"fc"|"Fc"				{print_to_log(); strcpy(yylval.string, yytext);   return (FC_);}
"YFLUX"|"YFLUX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_);}
"YFLUX_BV"|"YFLUX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_BV_);}
"YFLUX_HOR"|"YFLUX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_HOR_);}
"YFLUX_ORR"|"YFLUX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_ORR_);}
"YFLUX_H2O_ANODE"|"YFLUX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_H2O_ANODE_);}
"YFLUX_H2O_CATHODE"|"YFLUX_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_H2O_CATHODE_);}
"YFLUX_CONST"|"YFLUX_CONST_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_CONST_);}
"YFLUX_EQUIL"|"YFLUX_EQUIL_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_EQUIL_);}
"YFLUX_SUS"|"YFLUX_SUS_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_SUS_);}
"YFLUX_USER"|"YFLUX_USER_BC"		{print_to_log(); strcpy(yylval.string, yytext);   return(YFLUX_USER_);}
"YTOTALFLUX_CONST"|"YTOTALFLUX_CONST_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(YTOTALFLUX_CONST_);}  /*  RSL 6/7/00  */
"YUSER"|"YUSER_BC"			{print_to_log(); strcpy(yylval.string, yytext);   return(YUSER_);}
"YX"					{print_to_log(); strcpy(yylval.string, yytext);   return(YX_);}
"YZ"					{print_to_log(); strcpy(yylval.string, yytext);   return(YZ_);}
"Y_DISCONTINUOUS"|"Y_DISCONTINUOUS_BC"	{print_to_log(); strcpy(yylval.string, yytext);   return(Y_DISCONTINUOUS_);}
"Yield"					{print_to_log(); strcpy(yylval.string, yytext);   return(YEILD_);}
"Z" 					{print_to_log(); strcpy(yylval.string, yytext);   return(Z_);}
"ZX" 					{print_to_log(); strcpy(yylval.string, yytext);   return(ZX_);}
"ZY" 					{print_to_log(); strcpy(yylval.string, yytext);   return(ZY_);}
"auxiliary"|"AUXILIARY"|"Auxiliary"	{print_to_log(); strcpy(yylval.string, yytext);   return(AUXILIARY_);}  
"bicgstab"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(BICGSTAB_);}
"bilu"|"BILU"|"Bilu" 			{print_to_log(); strcpy(yylval.string, yytext);   return(BILU_);}
"cg"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(CG_);}
"cgs"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(CGS_);}
"atm"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(ATM_);}
"torr"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(TORR_);}
"filter"|"FILTER"|"Filter"		{print_to_log(); strcpy(yylval.string, yytext);   return(FILTER_);}
"concentration"|"Concentration"|"CONCENTRATION"		{print_to_log(); strcpy(yylval.string, yytext);   return(CONCENTRATION_);}
"continuity"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(CONTINUITY_);}
"CONTINUATION"|"Continuation"|"continuation"		{print_to_log(); strcpy(yylval.string, yytext);   return(CONTINUATION_);}
"convection"			 	{print_to_log(); strcpy(yylval.string, yytext);   return(CONVECTION_);}
"correction"|"CORRECTION"|"Correction"	{print_to_log(); strcpy(yylval.string, yytext);   return(CORRECTION_);} 
"cylindrical"|"CYLINDRICAL"		{print_to_log(); strcpy(yylval.string, yytext);   return(CYLINDRICAL_);}
"deformable"			 	{print_to_log(); strcpy(yylval.string, yytext);   return(DEFORMABLE_);}  
"delta_t"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DELTA_T_);}  
"density"				{print_to_log(); strcpy(yylval.string, yytext);   return(DENSITY_);}
"drop"|"DROP"|"Drop"			{print_to_log(); strcpy(yylval.string, yytext);   return(DROP_);}  
"energy"|"Energy" 			{print_to_log(); strcpy(yylval.string, yytext);   return(ENERGY_);}
"error"|"Error"|"ERROR"	 		{print_to_log(); strcpy(yylval.string, yytext);   return(ERROR_);}  
"factor"|"Factor"|"FACTOR" 		{print_to_log(); strcpy(yylval.string, yytext);   return(FACTOR_);} 
"factorization"|"FACTORIZATION"|"Factorization"		{print_to_log(); strcpy(yylval.string, yytext);   return(FACTORIZATION_);} 
"format"|"FORMAT"|"Format"		{print_to_log(); strcpy(yylval.string, yytext);   return(FORMAT_);}
"front"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(FRONT_);}
"gmres"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(GMRES_);}
"gradient11"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(GRADIENT11_);}
"gradient12"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(GRADIENT12_);}
"gradient13"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(GRADIENT13_);}
"gradient21"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(GRADIENT21_);}
"gradient22"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(GRADIENT22_);}
"gradient23"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(GRADIENT23_);}
"gradient31"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(GRADIENT31_);}
"gradient32"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(GRADIENT32_);}
"gradient33"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(GRADIENT33_);}
"ilu"|"ILU"|"Ilu" 			{print_to_log(); strcpy(yylval.string, yytext);   return(ILU_);}
"rilu"|"RILU"|"Rilu" 			{print_to_log(); strcpy(yylval.string, yytext);   return(RILU_);}
"ilut"|"ILUT"|"Ilut" 			{print_to_log(); strcpy(yylval.string, yytext);   return(ILUT_);}
"idd"|"ICC"|"Icc" 			{print_to_log(); strcpy(yylval.string, yytext);   return(ICC_);}
"in" 					{print_to_log(); strcpy(yylval.string, yytext);   return(IN_);}
"TRACKING"|"Tracking"|"tracking"	{print_to_log(); strcpy(yylval.string, yytext);   return(TRACKING_);}
"LENGTH"|"Lenght"|"length"		{print_to_log(); strcpy(yylval.string, yytext);   return(LENGTH_);}
"SCALE"|"Scale"|"scale"			{print_to_log(); strcpy(yylval.string, yytext);   return(SCALE_);}
"IMPLICIT"|"Implicit"|"implicit"	{print_to_log(); strcpy(yylval.string, yytext);   return(IMPLICIT_);}
"EMBEDDED"|"Embedded"|"embedded"	{print_to_log(); strcpy(yylval.string, yytext);   return(EMBEDDED_);}
"INITIALIZATION"|"Initialization"|"initialization"	{print_to_log(); strcpy(yylval.string, yytext);   return(INITIALIZATION_);}
"METHOD"|"Method"|"method"		{print_to_log(); strcpy(yylval.string, yytext);   return(METHOD_);}
"PROJECT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(PROJECT_);}
"EXO_READ"				{print_to_log(); strcpy(yylval.string, yytext);   return(EXO_READ_);}
"SURFACES"|"Surfaces"|"surfaces"	{print_to_log(); strcpy(yylval.string, yytext);   return(SURFACES_);}
"CONTROL"|"Control"|"control"		{print_to_log(); strcpy(yylval.string, yytext);   return(CONTROL_);}
"WIDTH"|"Width"|"width"			{print_to_log(); strcpy(yylval.string, yytext);   return(WIDTH_);}
"RENORMALIZATION"|"Renormalization"|"renormalization"		{print_to_log(); strcpy(yylval.string, yytext);   return(RENORMALIZATION_);}
"HUYGENS" 				{print_to_log(); strcpy(yylval.string, yytext);   return(HUYGENS_);}
"HUYGENS_CONSTRAINED" 			{print_to_log(); strcpy(yylval.string, yytext);   return(HUYGENS_CONSTRAINED_);}
"intermediate"|"INTERMEDIATE"|"Intermediate"	{print_to_log(); strcpy(yylval.string, yytext);   return(INTERMEDIATE_);}
"INTERFACE"|"Interface"|"interface"	{print_to_log(); strcpy(yylval.string, yytext);   return(INTERFACE_);}
"isoparametric" 			{print_to_log(); strcpy(yylval.string, yytext);   return(ISOPARAMETRIC_);}
"level"|"Level"|"LEVEL"			{print_to_log(); strcpy(yylval.string, yytext);   return(LEVEL_);}
"level_set"|"Level_Set"|"LEVEL_SET"	{print_to_log(); strcpy(yylval.string, yytext);   return(LEVEL_SET_);}
"set"|"SET"|"Set"			{print_to_log(); strcpy(yylval.string, yytext);   return(SET_);}
"foam"|"FOAM"|"Foam"			{print_to_log(); strcpy(yylval.string, yytext);   return(FOAM_);}
"ls"|"LS"|"Ls"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(LS_);}
"lu"|"LU"|"Lu"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(LU_);}
"ma28"|"MA28"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MA28_);}
"media"|"Media"|"MEDIA"			{print_to_log(); strcpy(yylval.string, yytext);   return(MEDIA_);}
"mesh1"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(MESH1_);}
"mesh2"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(MESH2_);}
"mesh3"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(MESH3_);}
"mom_solid1"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOM_SOLID1_);}
"mom_solid2"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOM_SOLID2_);}
"mom_solid3" 				{print_to_log(); strcpy(yylval.string, yytext);   return(MOM_SOLID3_);}
"momentum1"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOMENTUM1_);}
"momentum2"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOMENTUM2_);}
"momentum3"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOMENTUM3_);}
"norm"			 		{print_to_log(); strcpy(yylval.string, yytext);   return(NORM_);} 
"normal"		 		{print_to_log(); strcpy(yylval.string, yytext);   return(NORMAL_);}  
"of" 					{print_to_log(); strcpy(yylval.string, yytext);   return(OF_);}
"order"|"ORDER"|"Order"	 		{print_to_log(); strcpy(yylval.string, yytext);   return(ORDER_);} 
"output"|"Output"|"OUTPUT"		{print_to_log(); strcpy(yylval.string, yytext);   return(OUTPUT_);}  
"overlap"			 	{print_to_log(); strcpy(yylval.string, yytext);   return(OVERLAP_);} 
"parameter"|"Parameter"|"PARAMETER"	{print_to_log(); strcpy(yylval.string, yytext);   return(PARAMETER_);}
"parameters"|"Parameters"|"PARAMETERS"	{print_to_log(); strcpy(yylval.string, yytext);   return(PARAMETERS_);}
"phase"				 	{print_to_log(); strcpy(yylval.string, yytext);   return(PHASE_);} 
"pmomentum1"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(PMOMENTUM1_);}
"pmomentum2"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(PMOMENTUM2_);}
"pmomentum3"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(PMOMENTUM3_);}
"polynomial"|"Polynomial"|"POLYNOMIAL"	{print_to_log(); strcpy(yylval.string, yytext);   return(POLYNOMIAL_);} 
"porous"|"Porous"|"POROUS"		{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_);}
"pressure"|"Pressure"|"PRESSURE"	{print_to_log(); strcpy(yylval.string, yytext);   return(PRESSURE_);} 
"processors"|"PROCESSORS"|"Processors"	{print_to_log(); strcpy(yylval.string, yytext);   return(PROCESSORS_);} 
"residual"|"Residual"|"RESIDUAL"	{print_to_log(); strcpy(yylval.string, yytext);   return(RESIDUAL_);} 
"results"				{print_to_log(); strcpy(yylval.string, yytext);   return(RESULTS_);}
"reuse"|"Resue"|"REUSE"			{print_to_log(); strcpy(yylval.string, yytext);   return(REUSE_);} 
"recalc"|"Recalc"|"RECALC"		{print_to_log(); strcpy(yylval.string, yytext);   return(RECALC_);}  
"calc"|"Calc"|"CALC"			{print_to_log(); strcpy(yylval.string, yytext);   return(CALC_);}   
"shear_rate"|"SHEAR_RATE"		{print_to_log(); strcpy(yylval.string, yytext);   return(SHEAR_RATE_);}
"species"|"Species"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_);}
"species_bulk"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_BULK_);}
"species_surf"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_SURF_);}
"species_unknown"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNKNOWN_);}
"SPECIES_UNK_0"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_0_);}
"SPECIES_UNK_1"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_1_);}
"SPECIES_UNK_2"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_2_);}
"SPECIES_UNK_3"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_3_);}
"SPECIES_UNK_4"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_4_);}
"SPECIES_UNK_5"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_5_);}
"SPECIES_UNK_6"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_6_);}
"SPECIES_UNK_7"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_7_);}
"SPECIES_UNK_8"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_8_);}
"SPECIES_UNK_9"				{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_9_);}
"SPECIES_UNK_10"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_10_);}
"SPECIES_UNK_11"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_11_);}
"SPECIES_UNK_12"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_12_);}
"SPECIES_UNK_13"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_13_);}
"SPECIES_UNK_14"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_14_);}
"SPECIES_UNK_15"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_15_);}
"SPECIES_UNK_16"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_16_);}
"SPECIES_UNK_17"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_17_);}
"SPECIES_UNK_18"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_18_);}
"SPECIES_UNK_19"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_19_);}
"SPECIES_UNK_20"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_20_);}
"SPECIES_UNK_21"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_21_);}
"SPECIES_UNK_22"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_22_);}
"SPECIES_UNK_23"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_23_);}
"SPECIES_UNK_23"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_24_);}
"SPECIES_UNK_25"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_25_);}
"SPECIES_UNK_26"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_26_);}
"SPECIES_UNK_27"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_27_);}
"SPECIES_UNK_28"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_28_);}
"SPECIES_UNK_29"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_29_);}
"SPECIES_UNK_LAST"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_UNK_LAST_);}
"stages"|"STAGES"|"Stages"		{print_to_log(); strcpy(yylval.string, yytext);   return(STAGES_);}
"step"|"Step"|"STEP" 			{print_to_log(); strcpy(yylval.string, yytext);   return(STEP_);}  
"LIMIT"|"Limit"|"limit"			{print_to_log(); strcpy(yylval.string, yytext);   return(LIMIT_);}
"storage"|"Storage"|"STORAGE"		{print_to_log(); strcpy(yylval.string, yytext);   return(STORAGE_);}
"stress11"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(STRESS11_);}
"stress12"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(STRESS12_);}
"stress13"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(STRESS13_);}
"stress22"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(STRESS22_);}
"stress23"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(STRESS23_);}
"stress33"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(STRESS33_);}
"subparametric" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SUBPARAMETRIC_);}
"subspace"|"Subspace"|"SUBSPACE"	{print_to_log(); strcpy(yylval.string, yytext);   return(SUBSPACE_);}  
"sym_GS"|"SYM_GS" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SYM_GS_);}
"sym_row_sum"|"SYM_ROW_SUM"		{print_to_log(); strcpy(yylval.string, yytext);   return(SYM_ROW_SUM_);}
"row_sym"|"ROW_SYM" 			{print_to_log(); strcpy(yylval.string, yytext);   return(ROW_SYM_);}
"sym_diag"|"SYM_DIAG" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SYM_DIAG_);}
"tfqmr"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(TFQMR_);}
"time"|"Time"				{print_to_log(); strcpy(yylval.string, yytext);   return(TIME_);}
"tolerance"|"Tolerance"|"TOLERANCE"	{print_to_log(); strcpy(yylval.string, yytext);   return(TOLERANCE_);}  
"type"|"Type"|"TYPE"			{print_to_log(); strcpy(yylval.string, yytext);   return(TYPE_);}  
"vbr"|"VBR"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(VBR_);}
"msr"|"MSR"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MSR_);}
"vectors"|"VECTORS"|"Vectors"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(VECTORS_);}
"voltage"|"VOLTAGE" 			{print_to_log(); strcpy(yylval.string, yytext);   return(VOLTAGE_);}
"y12m"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(Y12M_);}
"dom_decomp"|"DOM_DECOMP"		{print_to_log(); strcpy(yylval.string, yytext);   return(DOM_DECOMP_);}	  
"VELOCITY1"|"Velocity1"|"velocity1"	{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY1_);}  
"VELOCITY2"|"Velocity2"|"velocity2"	{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY2_);}  
"VELOCITY3"|"Velocity3"|"velocity3"	{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY3_);} 
"VOLF_PHASE_0"|"volf_phase_0"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(VOLF_PHASE_0_);}
"VOLF_PHASE_1"|"volf_phase_1"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(VOLF_PHASE_1_);}
"VOLF_PHASE_2"|"volf_phase_2"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(VOLF_PHASE_2_);}
"VOLF_PHASE_3"|"volf_phase_3"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(VOLF_PHASE_3_);}
"VOLF_PHASE_4"|"volf_phase_4"	 	{print_to_log(); strcpy(yylval.string, yytext);   return(VOLF_PHASE_4_);}
"VOLF_PHASE_LAST"|"volf_phase_last" 	{print_to_log(); strcpy(yylval.string, yytext);   return(VOLF_PHASE_LAST_);}
"MASS_FRACTION"|"mass_fraction"|"Mass_Fraction"		{print_to_log(); strcpy(yylval.string, yytext);   return(MASS_FRACTION_);}  
"MESH_DISPLACEMENT1" 			{print_to_log(); strcpy(yylval.string, yytext);   return( MESH_DISPLACEMENT1_);} 
"MESH_DISPLACEMENT2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(MESH_DISPLACEMENT2_);}  
"MES_DISPALCEMENT3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(MESH_DISPLACEMENT3_);}  
"POLYMER_STRESS11" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS11_);}  
"POLYMER_STRESS12" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS12_);}  
"POLYMER_STRESS13" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS13_);}  
"POLYMER_STRESS22" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS22_);}  
"POLYMER_STRESS23" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS23_);}  
"POLYMER_STRESS33" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS33_);}  
"SOLID_DISPLACEMENT1" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_DISPLACEMENT1_);}  
"DMX_RS" 				{print_to_log(); strcpy(yylval.string, yytext);   return(DMX_RS_);}  
"SOLID_DISPLACEMENT2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_DISPLACEMENT2_);}  
"SOLID_DISPLACEMENT3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_DISPLACEMENT3_);} 
"VALUE"|"Value"|"value"			{print_to_log(); strcpy(yylval.string, yytext);   return(VALUE_);}  
"VELOCITY_GRADIENT11" 			{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY_GRADIENT11_);}  
"VELOCITY_GRADIENT12" 			{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY_GRADIENT12_);}  
"VELOCITY_GRADIENT22" 			{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY_GRADIENT22_);}  
"VELOCITY_GRADIENT13" 			{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY_GRADIENT13_);}  
"VELOCITY_GRADIENT23"			{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY_GRADIENT23_);} 
"VEOLCITY_GRADIENT31" 			{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY_GRADIENT31_);} 
"VEOLCITY_GRADIENT32" 			{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY_GRADIENT32_);}  
"VEOLCITY_GRADIENT33" 			{print_to_log(); strcpy(yylval.string, yytext);   return(VELOCITY_GRADIENT33_);}  
"PVELOCITY1"|"pvelocity1"		{print_to_log(); strcpy(yylval.string, yytext);   return(PVELOCITY1_);} 
"PVELOCITY2"|"pvelocity2"		{print_to_log(); strcpy(yylval.string, yytext);   return(PVELOCITY2_);} 
"PVELOCITY3"|"pvelocity3"		{print_to_log(); strcpy(yylval.string, yytext);   return(PVELOCITY3_);} 
"POLYMER_STRESS11_1" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS11_1_);}
"POLYMER_STRESS11_2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS12_1_);}
"POLYMER_STRESS22_1" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS22_1_);}
"POLYMER_STRESS13_1" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS13_1_);}
"POLYMER_STRESS23_1" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS23_1_);}
"POLYMER_STRESS33_1" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS33_1_);}
"POLYMER_STRESS11_2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS11_2_);}
"POLYMER_STRESS12_2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS12_2_);}
"POLYMER_STRESS13_2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS13_2_);}
"POLYMER_STRESS23_2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS23_2_);}
"POLYMER_STRESS33_2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS33_2_);}
"POLYMER_STRESS11_3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS11_3_);}
"POLYMER_STRESS12_3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS12_3_);}
"POLYMER_STRESS22_3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS22_3_);}
"POLYMER_STRESS13_3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS13_3_);}
"POLYMER_STRESS23_3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS23_3_);}
"POLYMER_STRESS33_3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS33_3_);}
"POLYMER_STRESS11_4" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS11_4_);}
"POLYMER_STRESS12_4" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS12_4_);}
"POLYMER_STRESS13_4" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS13_4_);}
"POLYMER_STRESS23_4" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS23_4_);}
"POLYMER_STRESS33_4" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS33_4_);}
"POLYMER_STRESS11_5" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS11_5_);}
"POLYMER_STRESS12_5" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS12_5_);}
"POLYMER_STRESS22_5" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS22_5_);}
"POLYMER_STRESS13_5" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS13_5_);}
"POLYMER_STRESS23_5" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS23_5_);}
"POLYMER_STRESS33_5" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS33_5_);}
"POLYMER_STRESS11_6" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS11_6_);}
"POLYMER_STRESS12_6" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS12_6_);}
"POLYMER_STRESS22_6" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS22_6_);}
"POLYMER_STRESS13_6" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS13_6_);}
"POLYMER_STRESS23_6" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS23_6_);}
"POLYMER_STRESS33_6" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS33_6_);}
"POLYMER_STRESS11_7" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS11_7_);}
"POLYMER_STRESS12_7"			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS12_7_);}
"POLYMER_STRESS22_7" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS22_7_);}
"POLYMER_STRESS13_7" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS13_7_);}
"POLYMER_STRESS23_7" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS23_7_);}
"POLYMER_STRESS33_7" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POLYMER_STRESS33_7_);}
"POR_LIQ_PRES"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(POR_LIQ_PRES_);}	
"POR_GAS_PRES"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(POR_GAS_PRES_);}		
"POR_POROSITY"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(POR_POROSITY_);}			
"POR_SATURATION" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POR_SATURATION_);}	
"POR_LAST" 				{print_to_log(); strcpy(yylval.string, yytext);   return(POR_LAST_);}	            
"VORT_DIR1" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VORT_DIR1_);}	           
"VORT_DIR2" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VORT_DIR2_);}	           
"VORT_DIR3" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VORT_DIR3_);}	            
"VORT_LAMBDA" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VORT_LAMBDA_);}	
"Species_Conc_0" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_0_);}  
"Species_Conc_1" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_1_);}  
"Species_Conc_2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_2_);}  
"Species_Conc_3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_3_);}  
"Species_Conc_4" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_4_);}  
"Species_Conc_5" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_5_);}  
"Species_Conc_6" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_6_);}  
"Species_Conc_7" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_7_);}  
"Species_Conc_8" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_8_);}  
"Species_Conc_9" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_9_);}  
"Species_Conc_10" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_10_);}
"Species_Conc_11" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_11_);}
"Species_Conc_12"			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_12_);}
"Species_Conc_13" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_13_);}
"Species_Conc_14" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_14_);}
"Species_Conc_15" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_15_);}
"Species_Conc_16" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_16_);}
"Species_Conc_17" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_17_);}
"Species_Conc_18" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_18_);}
"Species_Conc_19" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_19_);}
"Species_Conc_20" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_20_);}
"Species_Conc_21" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_21_);}
"Species_Conc_22" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_22_);}
"Species_Conc_23" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_23_);}
"Species_Conc_24" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_24_);}
"Species_Conc_25" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_25_);}
"Species_Conc_26" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_26_);}
"Species_Conc_27" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_27_);}
"Species_Conc_28" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_28_);}
"Species_Conc_29" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_CONC_29_);}
"VolFracPh0" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VOLFRACPH_0_);}  
"VolFracPh1" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VOLFRACPH_0_);}  
"VolFracPh2" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VOLFRACPH_0_);}  
"VolFracPh3" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VOLFRACPH_0_);}  
"VolFracPh3" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VOLFRACPH_0_);}  
"MESH_POSITION1" 			{print_to_log(); strcpy(yylval.string, yytext);   return(MESH_POSITION1_);} 
"MESH_POSITION2" 			{print_to_log(); strcpy(yylval.string, yytext);   return(MESH_POSITION2_);} 
"MESH_POSITION3" 			{print_to_log(); strcpy(yylval.string, yytext);   return(MESH_POSITION3_);} 
"VEL_NORM" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VEL_NORM_);}  
"D_VEL1_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_VEL1_DT_);} 
"UDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(UDOT_);}  
"D_VEL2_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_VEL2_DT_);}  
"VDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(VDOT_);}  
"D_VEL3_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_VEL3_DT_);}  
"WDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(WDOT_);}  
"D_T_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_T_DT_);}  
"TDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(TDOT_);}  
"D_C_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_C_DT_);}  
"TDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(TDOT_);}  
"D_C_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_C_DT_);}  
"CDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(CDOT_);}  
"D_X1_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_X1_DT_);}  
"XDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(XDOT_);}  
"D_X2_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_X2_DT_);}  
"VDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(VDOT_);}  
"D_X3_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_X3_DT_);}  
"ZDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(ZDOT_);}  
"D_S_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_S_DT_);}  
"SDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(SDOT_);}  
"D_P_DT" 				{print_to_log(); strcpy(yylval.string, yytext);   return(D_P_DT_);}  
"PDOT" 					{print_to_log(); strcpy(yylval.string, yytext);   return(PDOT_);}  
"QUADRATIC" 				{print_to_log(); strcpy(yylval.string, yytext);   return(QUADRATIC_);}  
"XY"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(XY_);} 
"XZ"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(XZ_);}
"EVAPORATION"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(EVAPORATION_);}
"Nodeset"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(NODESET_);}
"EVAPORATION"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(EVAPORATION_);}
"DISSOLUTION"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DISSOLUTION_);}
"LIQUID_VAPOR"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(LIQUID_VAPOR_);}
"SOLID_LIQUID"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_LIQUID_);}
"zero"|"ZERO"|"Zero" 			{print_to_log(); strcpy(yylval.string, yytext);   return(ZERO_);}
"random"|"RANDOM" 			{print_to_log(); strcpy(yylval.string, yytext);   return(RANDOM_);}
"one"|"ONE"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(ONE_);}
"sines"|"SINES"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SINES_);}
"read"|"READ"|"Read" 			{print_to_log(); strcpy(yylval.string, yytext);   return(READ_);}
"read_exoii"|"READ_EXOII"		{print_to_log(); strcpy(yylval.string, yytext);   return(READ_EXOII_);}
"read_exoii_file"|"READ_EXOii_FILE"	{print_to_log(); strcpy(yylval.string, yytext);   return(READ_EXOII_FILE_);}
"R0"|"r0"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(R0_);}
"RHS"|"rhs"|"Rhs" 			{print_to_log(); strcpy(yylval.string, yytext);   return(RHS_);}
"ANORM"|"Anorm"|"anorm"			{print_to_log(); strcpy(yylval.string, yytext);   return(ANORM_);}
"NORM"|"Norm"|"norm"			{print_to_log(); strcpy(yylval.string, yytext);   return(NORM_);}
"SOL"|"Sol"|"sol" 			{print_to_log(); strcpy(yylval.string, yytext);   return(SOL_);}
"NOSCLAED"|"Noscaled"|"noscaled"	{print_to_log(); strcpy(yylval.string, yytext);   return(NOSCALED_);}
"WEIGHTED"|"Weighted"|"weighted"	{print_to_log(); strcpy(yylval.string, yytext);   return(WEIGHTED_);}
"ALL"|"All"|"all"			{print_to_log(); strcpy(yylval.string, yytext);   return(ALL_);}
"LAST"|"Last"|"last"			{print_to_log(); strcpy(yylval.string, yytext);   return(LAST_);}
"WARNINGS"|"Warnings"|"warnings"	{print_to_log(); strcpy(yylval.string, yytext);   return(WARNINGS_);}
"GRAPH"|"Graph"|"graph"			{print_to_log(); strcpy(yylval.string, yytext);   return(GRAPH_);}
"FILLIN"|"Fillin"|"fillin"		{print_to_log(); strcpy(yylval.string, yytext);   return(FILLIN_);}
"FULL"|"Full"|"full"			{print_to_log(); strcpy(yylval.string, yytext);   return(FULL_);}
"DIAG"|"Diag"|"diag"			{print_to_log(); strcpy(yylval.string, yytext);   return(DIAG_);}
"STANDARD"|"Standard"|"standard"	{print_to_log(); strcpy(yylval.string, yytext);   return(STANDARD_);}
"SYMMETRIC"|"Symmetric"|"symmetric"	{print_to_log(); strcpy(yylval.string, yytext);   return(SYMMETRIC_);}
"RAND"|"Rand"|"rand"			{print_to_log(); strcpy(yylval.string, yytext);   return(RAND_);}
"RESID"|"Resid"|"resid"			{print_to_log(); strcpy(yylval.string, yytext);   return(RESID_);}
"REORDER"|"Reorder"|"reorder"		{print_to_log(); strcpy(yylval.string, yytext);   return(REORDER_);}
"RCM"|"Rcm"|"rcm"			{print_to_log(); strcpy(yylval.string, yytext);   return(RCM_);}
"SAVE"|"Save"|"save"			{print_to_log(); strcpy(yylval.string, yytext);   return(SAVE_);}
"RELAX"|"Relax"|"relax"			{print_to_log(); strcpy(yylval.string, yytext);   return(RELAX_);}
"threshold"|"THRESHOLD"|"Threshold"	{print_to_log(); strcpy(yylval.string, yytext);   return(THRESHOLD_);}
"ABSOLUTE"|"Absolute"|"absolute	"	{print_to_log(); strcpy(yylval.string, yytext);   return(ABSOLUTE_);}
"RELATIVE"|"Relative"|"relative"	{print_to_log(); strcpy(yylval.string, yytext);   return(RELATIVE_);}
"classic"|"Classic"|"CLASSIC"		{print_to_log(); strcpy(yylval.string, yytext);   return(CLASSIC_);}
"TOTAL"|"Total"|"total"			{print_to_log(); strcpy(yylval.string, yytext);   return(TOTAL_);}
"DIFFUSE"|"Diffuse"|"diffuse"		{print_to_log(); strcpy(yylval.string, yytext);   return(DIFFUSE_);}
"VISCOPLASTIC"|"Viscoplastic"|"viscoplastic"	{print_to_log(); strcpy(yylval.string, yytext);   return(VISCOPLASTIC_);}
"DEF"|"Def"|"def"			{print_to_log(); strcpy(yylval.string, yytext);   return(DEF_);}
"GRAD"|"Grad"|"grad"			{print_to_log(); strcpy(yylval.string, yytext);   return(GRAD_);}
"ZZ"|"Zz"|"zz"				{print_to_log(); strcpy(yylval.string, yytext);   return(ZZ_);}
"SOLVENTS"|"Solvents"|"solvents"	{print_to_log(); strcpy(yylval.string, yytext);   return(SOLVENTS_);}
"DARCY"|"Darcy"|"darcy"			{print_to_log(); strcpy(yylval.string, yytext);   return(DARCY_);}
"GRID"|"Grid"|"grid"			{print_to_log(); strcpy(yylval.string, yytext);   return(GRID_);}
"PECLET"|"Peclet"|"peclet"		{print_to_log(); strcpy(yylval.string, yytext);   return(PECLET_);}
"VORTICITY"|"Vorticity"|"vorticity"	{print_to_log(); strcpy(yylval.string, yytext);   return(VORTICITY_);}
"CROSS-STREAM"|"Cross-Stream"|"cross-stream"	{print_to_log(); strcpy(yylval.string, yytext);   return(CROSSSTREAM_);}
"ELEM"|"Elem"|"elem"				{print_to_log(); strcpy(yylval.string, yytext);   return(ELEM_);}
"AVERAGE_CONC"|"Average_conc"|"average_conc"	{print_to_log(); strcpy(yylval.string, yytext);   return(AVERAGE_CONC_);}
"SURF_DISSIP"|"surf_dissip"|"Surf_dissip"	{print_to_log(); strcpy(yylval.string, yytext);   return(SURF_DISSIP_);}
"AREA"|"Area"|"area"				{print_to_log(); strcpy(yylval.string, yytext);   return(AREA_);}
"VOL_REVOLUTION"|"Vol_revolution"|"vol_revolution"	{print_to_log(); strcpy(yylval.string, yytext);   return(VOL_REVOLUTION_);}
"TORQUE"|"Torque"|"torque"		{print_to_log(); strcpy(yylval.string, yytext);   return(TORQUE_);}
"PROFILE"|"Profile"|"profile"		{print_to_log(); strcpy(yylval.string, yytext);   return(PROFILE_);}
"SOLID_POSITION1"			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_POSITION1_);}
"SOLID_POSITION2"			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_POSITION2_);}
"SOLID_POSITION3"			{print_to_log(); strcpy(yylval.string, yytext);   return(SOLID_POSITION3_);}
"X_RS"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(X_RS_);}
"Y_RS"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(Y_RS_);}
"Z_RS"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(Z_RS_);}
"VDX"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(VDX_);}
"VDY"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(VDY_);}
"VDZ"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(VDZ_);}
"VLAMBDA" 				{print_to_log(); strcpy(yylval.string, yytext);   return(VLAMBDA_);}
"P_SAT"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(P_SAT_);}
"P_POR"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(P_POR_);}
"P_GAS"	 				{print_to_log(); strcpy(yylval.string, yytext);   return(P_GAS_);}
"USER_POST"				{print_to_log(); strcpy(yylval.string, yytext);   return(USER_POST_);}
"POR_POROSITY"|"POR_PORSITY"		{print_to_log(); strcpy(yylval.string, yytext);   return(POR_PORSITY_);}
"ERROR_ZZ_VEL"				{print_to_log(); strcpy(yylval.string, yytext);   return(ERROR_ZZ_VEL_);}
"ERROR_ZZ_Q" 				{print_to_log(); strcpy(yylval.string, yytext);   return(ERROR_ZZ_Q_);}
"ERROR_ZZ_P" 				{print_to_log(); strcpy(yylval.string, yytext);   return(ERROR_ZZ_P_);}
"SURFACE_VECTORS"			{print_to_log(); strcpy(yylval.string, yytext);   return(SURFACE_VECTORS_);}
"LAGRANGE_CONVECTION"			{print_to_log(); strcpy(yylval.string, yytext);   return(LAGRANGE_CONVECTION_);}
"STRAIN_TENSOR"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(STRAIN_TENSOR_);}
"TIME_DERIVATIVES"			{print_to_log(); strcpy(yylval.string, yytext);   return(TIME_DERIVATIVES_);}
"REAL_STRESS_TENSOR"			{print_to_log(); strcpy(yylval.string, yytext);   return(REAL_STRESS_TENSOR_);}
"CONDUCTION_VECTORS"			{print_to_log(); strcpy(yylval.string, yytext);   return(CONDUCTION_VECTORS_);}
"DIFFUSION_VECTORS"			{print_to_log(); strcpy(yylval.string, yytext);   return(DIFFUSION_VECTORS_);}
"MM_RESIDUALS"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MM_RESIDUALS_);}
"NS_RESIDUALS"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(NS_RESIDUALS_);}
"DIV_VELOCITY"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DIV_VELOCITY_);}
"THIRD_INVAR_STRAIN"			{print_to_log(); strcpy(yylval.string, yytext);   return(THIRD_INVAR_STRAIN_);}
"SEC_INVAR_STRAIN"	 		{print_to_log(); strcpy(yylval.string, yytext);   return(SECOND_INVAR_STRAIN_);}
"FIRST_INVAR_STRAIN"	 		{print_to_log(); strcpy(yylval.string, yytext);   return(FIRST_INVAR_STRAIN_);}
"FILL_CONT"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(FILL_CONT_);}
"PRESSURE_CONT"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(PRESSURE_CONT_);}
"MEAN_SHEAR"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MEAN_SHEAR_);}
"CROSS_STREAM_SHEAR"			{print_to_log(); strcpy(yylval.string, yytext);   return(CROSS_STREAM_SHEAR_);}
"STREAM_NORMAL_STRESS"			{print_to_log(); strcpy(yylval.string, yytext);   return(STREAM_NORMAL_STRESS_);}
"PARTICLES"|"Particles"|"particles"	{print_to_log(); strcpy(yylval.string, yytext);   return(PARTICLES_);}
"DISSIPATION"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(DISSIPATION_);}
"SPECIES_MASS"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIES_MASS_);}
"HEAT_ENERGY"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(HEAT_ENERGY_);}
"MOMENTUMX"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOMENTUMX_);}
"MOMENTUMY"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOMENTUMY_);}
"MOMENTUMZ"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(MOMENTUMZ_);}
"STRESS_TRACE"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(STRESS_TRACE_);}
"POSITIVE_FILL"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(POSITIVE_FILL_);}
"NEGATIVE_FILL"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(NEGATIVE_FILL_);}
"NONLINEAR_PLANE_STRAIN"		{print_to_log(); strcpy(yylval.string, yytext);   return(NONLINEAR_PLANE_STRAIN_);}
"INCOMP_3D"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(INCOMP_3D_);}
"plasticity"|"PLASTICITY"|"Plasticity"  {print_to_log(); strcpy(yylval.string, yytext);   return(PLASTICITY_);}
"EVP_HYPER"|"evp_hyper"			{print_to_log(); strcpy(yylval.string, yytext);   return(EVP_HYPER_);}
"DENSE_POWER_LAW" 			{print_to_log(); strcpy(yylval.string, yytext);   return(DENSE_POWER_LAW_);}
"SHEAR_HARDEN" 				{print_to_log(); strcpy(yylval.string, yytext);   return(SHEAR_HARDEN_);}
"NO_MODEL"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(NO_MODEL_);}
"CARREAU_WLF"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(CARREAU_WLF_);}
"CARREAU_SUSPENSION" 			{print_to_log(); strcpy(yylval.string, yytext);   return(CARREAU_SUSPENSION_);}
"POWERLAW_SUSPENSION" 			{print_to_log(); strcpy(yylval.string, yytext);   return(POWERLAW_SUSPENSION_);}
"SAT_DES"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(SAT_DES_);}
"COMP_LEVEL"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(COMP_LEVEL_);}
"CAP_PRES"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(CAP_PRES_);}
"BILINEAR"	 			{print_to_log(); strcpy(yylval.string, yytext);   return(BILINEAR_);} 
"NON_LINEAR"				{print_to_log(); strcpy(yylval.string, yytext);   return(NON_LINEAR_);}
"constant2"|"CONSTANT2"|"Constant2"     {print_to_log(); strcpy(yylval.string, yytext);   return(CONSTANT2_);}
"wlf"|"WLF"|"Wlf"			{print_to_log(); strcpy(yylval.string, yytext);   return(WLF_);}
"NS"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(NS_);} 
"EB"		 			{print_to_log(); strcpy(yylval.string, yytext);   return(EB_);} 
"surfaces"|"Surfaces"|"SURFACES"	{print_to_log(); strcpy(yylval.string, yytext);   return(WLF_);}		 			
"BIGX1"   				{print_to_log(); strcpy(yylval.string, yytext);   return(BIGX1_);} 
"X1"     				{print_to_log(); strcpy(yylval.string, yytext);   return(X1_);} 
"BIGX2"     				{print_to_log(); strcpy(yylval.string, yytext);   return(BIGX2_);} 
"X2"      				{print_to_log(); strcpy(yylval.string, yytext);   return(X2_);} 
"BIGX3"     				{print_to_log(); strcpy(yylval.string, yytext);   return(BIGX3_);} 
"X3"       				{print_to_log(); strcpy(yylval.string, yytext);   return(X3_);} 
"CURVATURE"    				{print_to_log(); strcpy(yylval.string, yytext);   return(CURVATURE_);}
"H"       				{print_to_log(); strcpy(yylval.string, yytext);   return(H_);}  
"ISOSURFACE"  				{print_to_log(); strcpy(yylval.string, yytext);   return(ISOSURFACE_);}  
"chemkin_mat"|"CHEMKIN_MAT"|"Chemkin_mat"	{print_to_log(); strcpy(yylval.string, yytext);   return(CHEMKIN_MAT_);}  
"goma_mat"|"GOMA_MAT"|"Goma_mat"	{print_to_log(); strcpy(yylval.string, yytext);   return(GOMA_MAT_);}  
"STEFAN_MAXWELL_CHARGED"       		{print_to_log(); strcpy(yylval.string, yytext);   return(STEFAN_MAXWELL_CHARGED_);}  
"FICKIAN_CHARGED"			{print_to_log(); strcpy(yylval.string, yytext);   return(FICKIAN_CHARGED_);}  
"ELECTRODE_KINETICS"       		{print_to_log(); strcpy(yylval.string, yytext);   return(ELECTRODE_KINETICS_);}  
"EPOXY_DEA"       			{print_to_log(); strcpy(yylval.string, yytext);   return(EPOXY_DEA_);} 
"JOULE"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(JOULE_);} 
"VISC_DISS"       			{print_to_log(); strcpy(yylval.string, yytext);   return(VISC_DISS_);} 
"Navier-Stokes"|"NAVIER-STOKES"		{print_to_log(); strcpy(yylval.string, yytext);   return(NAVIER_STOKES_);} 
"SUSPEND"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(SUSPEND_);} 
"BOUSSINESQ"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(BOUSSINESQ_);} 
"HUGGINS-HUGGINS"|"Flory-Huggins"    	{print_to_log(); strcpy(yylval.string, yytext);   return(FLORY_HUGGINS_);} 
"BOUSS_JXB"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(BOUSS_JXB_);} 
"Non-volatile"|"NON-VOLATILE"|"non-volatile"	{print_to_log(); strcpy(yylval.string, yytext);   return(NON_VOLATILE_);} 
"Specific"|"SPECIFIC"   	    	{print_to_log(); strcpy(yylval.string, yytext);   return(SPECIFIC_);} 
"BOUSS_JXB"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(BOUSS_JXB_);} 
"molar"|"Molar"|"MOLAR"			{print_to_log(); strcpy(yylval.string, yytext);   return(MOLAR_);} 
"Source"|"source"|"SOURCE"			{print_to_log(); strcpy(yylval.string, yytext);   return(SOURCE_);} 
"NON-CONDENSABLE"|"non-condensable"|"Non-condensable"	{print_to_log(); strcpy(yylval.string, yytext);   return(NON_CONDENSABLE_);} 
"molecular"|"Molecular"|"MOLECULAR"			{print_to_log(); strcpy(yylval.string, yytext);   return(MOLECULAR_);} 
"POROUS_PART_SAT"   	    		{print_to_log(); strcpy(yylval.string, yytext);   return(POROUS_PART_SAT_);} 
"continuous"|"Continuous"|"CONTINUOUS"	{print_to_log(); strcpy(yylval.string, yytext);   return(CONTINUOUS_);} 
"SOLIDIFICATION"   	    		{print_to_log(); strcpy(yylval.string, yytext);   return(SOLIDIFICATION_);} 
"KOZENY_CARMAN"   	    		{print_to_log(); strcpy(yylval.string, yytext);   return(KOZENY_CARMAN_);}
"TENSOR"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(TENSOR_);}  
"flowingliquid"|"FlowingLiquid"|"FLOWINGLIQUID"	{print_to_log(); strcpy(yylval.string, yytext);   return(FLOWINGLIQUID_);} 
"molten_glass"|"Molten_Glass"|"MOLTEN_GLASS"	{print_to_log(); strcpy(yylval.string, yytext);   return(MOLTEN_GLASS_);} 
"DARCY_FICKIAN"   	    		{print_to_log(); strcpy(yylval.string, yytext);   return(DARCY_FICKIAN_);} 
"GENERALIZED_FICKIAN"   	    	{print_to_log(); strcpy(yylval.string, yytext);   return(GENERALIZED_FICKIAN_);} 
"STEFAN_MAXWELL"   	    		{print_to_log(); strcpy(yylval.string, yytext);   return(STEFAN_MAXWELL_);} 
"STEFAN_MAXWELL_VOLUME"   	    	{print_to_log(); strcpy(yylval.string, yytext);   return(STEFAN_MAXWELL_VOLUME_);} 
"HYDRO"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(HYDRO_);} 
"SUSP_BAL"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(SUSP_BAL_);} 
"GENERALIZED_FREE_VOL"   	    	{print_to_log(); strcpy(yylval.string, yytext);   return(GENERALIZED_FREE_VOL_);} 
"GENERALIZED"   	    		{print_to_log(); strcpy(yylval.string, yytext);   return(GENERALIZED_);} 
"KELVIN"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(KELVIN_);}
"FLAT"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(FLAT_);}
"ANTOINE"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(ANTOINE_);}
"RIEDEL"   	    			{print_to_log(); strcpy(yylval.string, yytext);   return(RIEDEL_);}
"PROJECTED_CARTESIAN"  			{print_to_log(); strcpy(yylval.string, yytext);   return(PROJECTED_CARTESIAN_);}
"Really meant number of bulk species equations in line above" {print_to_log(); strcpy(yylval.string, yytext);   return(REALLY_MEANT_);}

"Example"	  			{print_to_log(); strcpy(yylval.string, yytext);   return(EXAMPLE_TOKEN2_);} /*EXAMPLE CARD TOKEN 1 DEFINITION*/
"Card"		  			{print_to_log(); strcpy(yylval.string, yytext);   return(EXAMPLE_TOKEN2_);} /*EXAMPLE CARD TOKEN 2 DEFINITION*/




[+-]?[0-9]+				{print_to_log(); strcpy(yylval.string, yytext); return (INTEGER_);}

[+-]?(([0-9]*\.[0-9]*)|([0-9]+([eE][-+]?[0-9]+)?)|(([0-9]*\.[0-9]*)([eE][-+]?[0-9]+)?))	{print_to_log(); strcpy(yylval.string, yytext);   return (FLOAT_);}

^\n					{
					/*print_to_log();*/		/* blank line */ 
					/*line_number++;*/
					/*if (ProcID == 0) fprintf(parser_log, "%d ", line_number);*/
					/*return(BLANK_LINE_);*/
					
					}
					
^(\n)+					{
					/*print_to_log();*/		/* blank line */ 
					/*line_number++;*/
					/*if (ProcID == 0) fprintf(parser_log, "%d ", line_number);*/
					
					
				/*	return(BLANK_LINE_);*/
					}					

^(" "|\t)(" "|\t)*\n			{
					/*print_to_log();*/		/* another blank line */
					/*line_number++;
					if (ProcID == 0) fprintf(parser_log, "%d ", line_number);*/
					
										
				/*	return(BLANK_LINE_);*/
					}					

^[-|*|$|%|#|/].*\n       		{
					print_to_log(); /* a commnet */
					line_number++;
					if (ProcID == 0) fprintf(parser_log, "%d ", line_number);
					}	/* commnet line */
			

^(" "|\t)*[-|*|$|%|#|/].*\n		{
					print_to_log(); /* a commnet */
					line_number++;
					if (ProcID == 0) fprintf(parser_log, "%d ", line_number);
					}	/* commnet line */
					
			
^(" "|\t)*"Material Data".*\n		{
					print_to_log();
					line_number++;
					if (ProcID == 0) fprintf(parser_log, "%d ", line_number);
					}	/* commnet line */



[a-zA-Z_][0-9a-zA-Z_@\.]*		{print_to_log(); strcpy(yylval.string, yytext);  return (STRING_);}


\n					{ 
					  count_of_characters_in_this_line = 0;
					  if(error_found_in_last_line)
					  { 
					    if (ProcID == 0) 
					    {	  
					      /*fprintf(parser_log,"\n <<< ERROR %d Line %d: near \"%s\" >>>\n", error_number, line_number, debug_clue);
					      fprintf(error_log,"\n<<< ERROR %d Line %d: near \"%s\" >>>", error_number, line_number, debug_clue);						  
					      fprintf(stdout,"\n <<< ERROR %d Line %d: near \"%s\" >>>", error_number, line_number, debug_clue);*/
					      sprintf(msg,"Error in line %d near %s",line_number, debug_clue);
					      declare_error(msg);
					    }  
					    error_found_in_last_line = 0;
					    /*if (ProcID == 0) fprintf(parser_log,"%d ",line_number);*/
					  } 	
					  			    
					  return(CR_);
					}

					
									
\t					{
					if (ProcID == 0) fprintf(parser_log, " ");
					count_of_characters_in_this_line++;
					while ( (count_of_characters_in_this_line % 8) > 0)
						{ 
						if (ProcID == 0) fprintf(parser_log," ");
						count_of_characters_in_this_line++;
						}
					}
   
.					{print_to_log();}

%%

/*
<COMMENT>([a-zA-Z_]|.)?([0-9a-zA-Z_@\.]|.)?\n				{printf("%s", yytext);   BEGIN 0;}
*/

myinput(char *buf, int max)
{ 	  
  if(mystring && strlen(mystring))
    {

      *buf = *mystring++;

      return 1;
	/*
      if (len < strlen(mystring))
	{	
	  memcpy(buf, mystring++, sizeof(char));
	  len++;
	  return 1;
	}
      else 
	return 0;
	*/
    }
  else
    {
     /* fflush(stdout);
      printf("Uh oh! I was passed a null string!\n");*/
      fflush(stdout);
      return 0;
    }
}

void print_to_log()
{
 	if (ProcID == 0) fprintf(parser_log, "%s", yytext);
	count_of_characters_in_this_line = count_of_characters_in_this_line + strlen(yytext);
	if ( (strcmp(yytext,"\n")) && (strcmp(yytext," ")))
	{
		strcpy(debug_clue,yytext);
	}

}





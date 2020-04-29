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
 *$Id: mm_porous_EOS.c,v 5.0 2006-07-06 16:18:51 edwilke Exp $
 */


/* GOMA include files */

#include "rf_fem_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_fill_porous.h"

/*********** R O U T I N E S   I N   T H I S   F I L E *************************
*
*       NAME                    TYPE                    CALLED_BY
*    ------------             ---------               --------------
*
*     load_enthalpy            void        
*
******************************************************************************/
/******************************************************************************/

/***************************************************************************/


void
load_enthalpy(double saturation, double pressure)

    /*************************************************************************
     * load_enthalpy -- calculate enthalpy of liquid and gas phases in a porous
     *                    media from the unknowns used in the problem
     *
     *  input:   Assume that load_fv and load_fv_grads and load_fv_mesh_derivs
     *  -----    have all been called already, so unknowns and their
     *           sensitivies are known at this gauss point.
     *
     *  output:  calculates the enthalpy and its first and second
     *  -----    derivatives with respect to all the problem unknowns
     *           For transient calculations, the enthalpy at the
     *           old time step must be recalculated as well, in order
     *           to be used in the capacitance term.
     *
     *      mp->enthalpy[0] -> liquid (assumes no dissolved air)
     *      mp->enthalpy[1] -> water in the gas phase
     *      mp->enthalpy[2] -> air in the gas phase
     *      mp->d_enthalpy
     *      mp->d_d_enthalpy
     *	    mp_old->enthalpy
     ************************************************************************/
{
   /*  Local automatic variables  */
   
   double Pg, Pl, T, Pvsat, dPvsat_dT, d_dPvsat_dT;
   double Pg_old, Pl_old, T_old;
   double dh_dP, dh_dT, d_dha_dT, d_dh_dP[2], d_dh_dT[2];


  /* 
   *  Compute liquid enthalpy, ignoring enthalpy of air dissolved in solution
   *  Note: Using the liquid  pressure here.
   */

   Pg = pressure;
   Pl = fv->p_liq;
   T = fv->T;
   Pl_old = fv_old->p_liq;
   if (pd->e[pg->imtrx][R_POR_GAS_PRES]) Pg_old = fv_old->p_gas; else Pg_old = pressure;
   T_old = fv_old->T;

   if ( saturation != 0.0) {
     pmv->enthalpy[0] = h_water_compressed_EOS(Pl, T, &dh_dP, &dh_dT,
                                               d_dh_dP, d_dh_dT);
     /* fill only the non-zero terms */
     pmv->d_enthalpy[0][POR_LIQ_PRES] = dh_dP;
     pmv->d_d_enthalpy[0][POR_LIQ_PRES][POR_LIQ_PRES] = d_dh_dP[0];

     pmv->d_enthalpy[0][POR_TEMP] = dh_dT;
     pmv->d_d_enthalpy[0][POR_LIQ_PRES][POR_TEMP] = d_dh_dP[1];
     pmv->d_d_enthalpy[0][POR_TEMP][POR_LIQ_PRES] = d_dh_dT[0];
     pmv->d_d_enthalpy[0][POR_TEMP][POR_TEMP] = d_dh_dT[1];

     if (pd->TimeIntegration == TRANSIENT) {
         pmv_old->enthalpy[0] = h_water_compressed_EOS(Pl_old, T_old, &dh_dP, &dh_dT,
                                                       d_dh_dP, d_dh_dT);  
         pmv_old->d_enthalpy[0][POR_LIQ_PRES] = dh_dP;                                              
         pmv_old->d_enthalpy[0][POR_TEMP] = dh_dT;
     }
   }

  /*
   *  Compute enthalpy of liquid in the gas phase.
   */
   if ( saturation < 1.0) {

     Pvsat = P_water_sat_EOS(T, &dPvsat_dT, &d_dPvsat_dT);
     /* P = MIN (P, Pvsat); */
     if(Pg >= Pvsat) {

        /* Use Pvsat as pressure; note h is f(T) only here! */

        pmv->enthalpy[1] = h_water_superheat_EOS(Pvsat, T, &dh_dP, &dh_dT,
                                                 d_dh_dP, d_dh_dT);  

        /* fill only the non-zero terms */
        pmv->d_enthalpy[1][POR_TEMP] = dh_dT + dh_dP*dPvsat_dT;
        pmv->d_d_enthalpy[1][POR_TEMP][POR_TEMP] = d_dh_dT[1] + d_dh_dT[0]*dPvsat_dT +
                                                   d_dh_dP[1]*dPvsat_dT + dh_dP*d_dPvsat_dT;
     }
     else {

        /* Use Pg as pressure; note h is f(Pg,T) here! */

        pmv->enthalpy[1] = h_water_superheat_EOS(Pg, T, &dh_dP, &dh_dT,
                                                 d_dh_dP, d_dh_dT);  
        /* fill only the non-zero terms */

        pmv->d_enthalpy[1][POR_TEMP] = dh_dT;
        pmv->d_d_enthalpy[1][POR_TEMP][POR_TEMP] = d_dh_dT[1];

         pmv->d_enthalpy[1][POR_GAS_PRES] = dh_dP;
         pmv->d_d_enthalpy[1][POR_GAS_PRES][POR_GAS_PRES] = d_dh_dP[0];
         pmv->d_d_enthalpy[1][POR_GAS_PRES][POR_TEMP] = d_dh_dP[1];
         pmv->d_d_enthalpy[1][POR_TEMP][POR_GAS_PRES] = d_dh_dT[0];
     }

     if (pd->TimeIntegration == TRANSIENT) {
       Pvsat = P_water_sat_EOS(T_old, &dPvsat_dT, &d_dPvsat_dT);
       /* P = MIN (P, Pvsat); */
       if(Pg_old >= Pvsat) {

         /* Use Pvsat as pressure; note h is f(T) only here! */
 
         pmv_old->enthalpy[1] = h_water_superheat_EOS(Pvsat, T_old, &dh_dP, &dh_dT,
                                                  d_dh_dP, d_dh_dT);  
         /* fill only the non-zero terms */
         pmv_old->d_enthalpy[1][POR_TEMP] = dh_dT + dh_dP*dPvsat_dT;
       }
       else {

        /* Use Pg as pressure; note h is f(Pg,T) here! */

          pmv_old->enthalpy[1] = h_water_superheat_EOS(Pg_old, T_old, &dh_dP, &dh_dT,
                                                       d_dh_dP, d_dh_dT);  
          /* fill only the non-zero terms */

          pmv_old->d_enthalpy[1][POR_TEMP] = dh_dT;

          pmv_old->d_enthalpy[1][POR_GAS_PRES] = dh_dP;
          
       }
     }
   }

  /*
   *  Compute enthalpy of air in the gas phase. Function of T only.
   */

   if ( saturation < 1.0 ) {
     pmv->enthalpy[2] =  h_air_EOS(T, &dh_dT, &d_dha_dT);

     /* fill only the non-zero terms */
     pmv->d_enthalpy[2][POR_LIQ_PRES] = 0.0; 
     pmv->d_enthalpy[2][POR_GAS_PRES] = 0.0; 
     pmv->d_enthalpy[2][POR_TEMP] = dh_dT; 
     pmv->d_d_enthalpy[2][POR_TEMP][POR_TEMP] = d_dha_dT;

     if (pd->TimeIntegration == TRANSIENT) { 
         pmv_old->enthalpy[2] =  h_air_EOS(T_old, &dh_dT, &d_dha_dT);
         pmv_old->d_enthalpy[2][POR_TEMP] = dh_dT; 
     }
   }

}  /* end   load_enthalpy  */


/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double calc_rho_gas(
         double Pg, 
         double T,
         double *Pvsat,
         double *dPvsat_dT,
         double *drho_gas_dPg,
         double *drho_gas_dT,
         double *d_drho_gas_dPg,
         double *d_drho_gas_dT,
         double *rho_wg_sat,
         double *drho_wgsat_dT,
         double *d_drho_wgsat_dT
)
{
 double Pa;
 double d_dPvsat_dT;
 double rho_gas, rho_air;
 double drho_wg_dPv, drho_wg_dT, d_drho_wg_dPv[2], d_drho_wg_dT[2]; 
 double drho_air_dP, drho_air_dT;
 double d_drho_air_dP[2], d_drho_air_dT[2];


/*  Calculate the gas density based on P(gas) = P(sat_vapor) + P(air)
 *  and thus rho(gas) = rho(sat_vapor) + rho(air)
 */

 *Pvsat = P_water_sat_EOS(T, dPvsat_dT, &d_dPvsat_dT);
  Pa = Pg - *Pvsat;

 *rho_wg_sat = rho_sat_water_vap_EOS(*Pvsat, T, &drho_wg_dPv, &drho_wg_dT,
                                    d_drho_wg_dPv, d_drho_wg_dT);

/* Note that rho_wg_sat is a function of Pvsat, which is a function of T only */

 *drho_wgsat_dT   = drho_wg_dT + drho_wg_dPv*(*dPvsat_dT);
 *d_drho_wgsat_dT = d_drho_wg_dT[1] + d_drho_wg_dT[0]*(*dPvsat_dT) +
                    d_drho_wg_dPv[1]*(*dPvsat_dT) + drho_wg_dPv*d_dPvsat_dT;

 rho_air = rho_air_EOS(Pa, T, *dPvsat_dT, d_dPvsat_dT,
                       &drho_air_dP, &drho_air_dT, d_drho_air_dP, d_drho_air_dT);


/* Note that d(rho_air)/dPg = d(rho_air)/dPa since dPa/dPg = 1 */

 rho_gas = *rho_wg_sat + rho_air;

 *drho_gas_dPg = drho_air_dP;
 *drho_gas_dT  = *drho_wgsat_dT + drho_air_dT;

 d_drho_gas_dPg[0] = d_drho_air_dP[0];
 d_drho_gas_dPg[1] = d_drho_air_dP[1];

 d_drho_gas_dT[0] = d_drho_air_dT[0];
 d_drho_gas_dT[1] = *d_drho_wgsat_dT + d_drho_air_dT[1];

 return ( rho_gas );
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double calc_Ywg(
         double rho_wg_sat,
         double drho_wgsat_dT,
         double d_drho_wgsat_dT, 
         double rho_gas,
         double drho_gas_dPg,
         double drho_gas_dT,
         double *d_drho_gas_dPg,
         double *d_drho_gas_dT,
         double *dYwg_dPg,
         double *dYwg_dT,
         double *d_dYwg_dPg,
         double *d_dYwg_dT
)
{
  double Ywg;

  Ywg = rho_wg_sat / rho_gas;

  *dYwg_dPg = -(Ywg/rho_gas)*drho_gas_dPg;
  *dYwg_dT  = drho_wgsat_dT/rho_gas - (Ywg/rho_gas)*drho_gas_dT; 

  d_dYwg_dPg[0] = -((*dYwg_dPg)/rho_gas)*drho_gas_dPg + 
                   (Ywg/rho_gas/rho_gas)*drho_gas_dPg*drho_gas_dPg -
                   (Ywg/rho_gas)*d_drho_gas_dPg[0];
  d_dYwg_dPg[1] = -((*dYwg_dT)/rho_gas)*drho_gas_dPg + 
                   (Ywg/rho_gas/rho_gas)*drho_gas_dT*drho_gas_dPg -
                   (Ywg/rho_gas)*d_drho_gas_dPg[1];
  d_dYwg_dT[0]  = -drho_wgsat_dT/rho_gas/rho_gas*drho_gas_dPg -
                   (*dYwg_dPg)/rho_gas*drho_gas_dT +
                   Ywg/rho_gas/rho_gas*drho_gas_dPg*drho_gas_dT -
                   Ywg/rho_gas*d_drho_gas_dT[0];
  d_dYwg_dT[1]  =  d_drho_wgsat_dT/rho_gas -
                   drho_wgsat_dT/rho_gas/rho_gas*drho_gas_dT -
                   (*dYwg_dT)/rho_gas*drho_gas_dT +
                   Ywg/rho_gas/rho_gas*drho_gas_dT*drho_gas_dT -
                   Ywg/rho_gas*d_drho_gas_dT[1];


 return ( Ywg );
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double rho_air_EOS(
          double Pa, 
          double T,
          double dPv_dT,
          double d_dPv_dT,
          double *drhoa_dP,
          double *drhoa_dT,
          double d_drhoa_dP[2],
          double d_drhoa_dT[2]
)

/*************************************************************************
 *
 * Equation of state for air density assuming ideal gas (as in FEHM)
 *
 * Units:  temperature, T   -  degrees Celsius
 *         air partial pressure,  Pa   -  Pascals
 *
 *         return density in kg/m^3
 *
 *************************************************************************/
  
{

  double C, TK, TK_sq, TK_cu;

   /* convert to MPa  */
 
  /*return ( 1.292864*(273.15/TK)*(Pa/0.101325) ); */

  C = 1.292864 * 273.15 / 0.101325 / 1.0e6;
  TK = T + 273.15;
  TK_sq = TK * TK;
  TK_cu = TK_sq * TK;

  *drhoa_dP = C / TK;
  *drhoa_dT = -C*Pa/TK_sq - C*dPv_dT/TK; 
  d_drhoa_dP[0] = 0.0;
  d_drhoa_dP[1] = -C/TK_sq;
  d_drhoa_dT[0] = -C/TK_sq;
  d_drhoa_dT[1] = 2.0*C*Pa/TK_cu + 2.0*C*dPv_dT/TK_sq - C*d_dPv_dT/TK;
  
 return ( C*Pa/TK );

}  
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double rho_sat_water_vap_EOS(
          double P_sat,
          double T,
          double *drho_wgsat_dP,
          double *drho_wgsat_dT,
          double d_drho_wgsat_dP[2],
          double d_drho_wgsat_dT[2]
)

/*************************************************************************
 *
 * Equation of state for density of saturated water vapor (as in FEHM)
 *
 * Units:  Pressure of saturated vapor at T, P_sat - Pascals
 *         temperature, T   -  degrees Celsius
 *
 *         return density in kg/m^3
 *
 *
 *************************************************************************/
  
{
 double f, g, P, dfdP, dfdT, dgdP, dgdT;
 double gsq, gcb;
 double df_d_dT[2], df_d_dP[2], dg_d_dT[2], dg_d_dP[2];

/* The following coefficients are for the pressure range .001 to 20 MPa *
 * and temperatures .5 to 360 degrees C                                 */

  static double Y[] = {0.13299942e-04,
                       0.10000000e+01,
                      -0.10000000e+01,
                      -0.56746973e-02,
                      -0.32791354e-06,
                       0.21636240e-08,
                      -0.38485869e-11,
                       0.40896880e-01,
                       0.27696827e-02,
		      -0.94741649e-04}; 
  static double Z[] = {0.12789230e+00,
                      -0.28996744e+00,
                       0.26873883e-02,
                       0.33783903e-04,
                       0.55690966e-02,
                       0.72603809e-05,
                      -0.44323127e-07,
                       0.49878874e-03,
                      -0.13186635e-04,
		       0.72041771e-06};
                      

/* The following coefficients are for the pressure range .001 to 20 MPa
   and temperatures 15 to 360 degrees C ...                         

  static double Y[] = {0.15089524e-05,
                       0.10000000e+01,
                      -0.10000000e+01,
                      -0.16676705e-02,
                       0.40111210e-07,
                       0.25625316e-10,
                      -0.40479650e-12,
                       0.43379623e-01,
                       0.24991800e-02,
                      -0.94755043e-04}; 
  static double Z[] = {0.12636224e+00,
                      -0.30463489e+00,
                       0.27981880e-02,
                       0.51132337e-05,
                       0.59318010e-02,
                       0.80972509e-05,
                      -0.43798358e-07,
                       0.53046787e-03,
                      -0.84916607e-05,
                       0.48444919e-06};   ...   end of "high-pressure" coef   */

  P = P_sat / 1.0e6;  /* function expects pressure in MPa  */

  /* return ( eval_rho_poly (Y,P_sat,T) 
           / eval_rho_poly (Z,P_sat,T) ); */

  f = eval_poly (Y,P,T);
  g = eval_poly (Z,P,T);
  gsq = g*g;
  gcb = gsq*g;


  dfdP = eval_poly_dP (Y,P,T,df_d_dP);
  dgdP = eval_poly_dP (Z,P,T,dg_d_dP);

  dfdT = eval_poly_dT (Y,P,T,df_d_dT);
  dgdT = eval_poly_dT (Z,P,T,dg_d_dT);

  *drho_wgsat_dP = (dfdP/g - f/g/g*dgdP) / 1.0e6;
  *drho_wgsat_dT = dfdT/g - f/g/g*dgdT;


  d_drho_wgsat_dP[0] = (df_d_dP[0]/g - 2.0*dfdP/gsq*dgdP + 2.0*f/gcb*dgdP*dgdP 
                           - f/gsq*dg_d_dP[0]) / 1.0e12;
  d_drho_wgsat_dP[1] = (df_d_dP[1]/g - dfdP*dgdT/gsq - dfdT/gsq*dgdP + 2.0*f/gcb*dgdT*dgdP
                           - f/gsq*dg_d_dP[1]) / 1.0e6;
  d_drho_wgsat_dT[0] = (df_d_dT[0]/g - dfdT/gsq*dgdP - dfdP/gsq*dgdT 
                           + 2.0*f/gcb*dgdP*dgdT - f/gsq*dg_d_dT[0]) / 1.0e6;
  d_drho_wgsat_dT[1] = df_d_dT[1]/g - 2.0*dfdT/gsq*dgdT + 2.0*f/gcb*dgdT*dgdT 
                           - f/gsq*dg_d_dT[1];

  return (  f/g  ); 

} 


/****************************************************************************************/

double P_water_sat_EOS(
          double T,
          double *dPvsat_dT,
          double *d_dPvsat_dT
)

/*************************************************************************
 *
 * Equation of state for vapor pressure of water (as in FEHM)
 *
 * Units:  temperature, T   -  degrees Celsius
 *
 *         Return pressure in Pascals
 *
 *
 *************************************************************************/
  
{
  double f, g, dfdT, dgdT, df_d_dT, dg_d_dT;
  double gsq, gcb;
  static int  order = 4;
  static double Y[] = {0.71725602e-03,
                       0.22607516e-04,
                       0.26178556e-05,
                      -0.10516335e-07,
                       0.63167028e-09};
  static double Z[] = {0.10000000e+01,
                      -0.22460012e-02,
                       0.30234492e-05,
                      -0.32466525e-09,
                       0.0}; 
 
  /* return ( eval_pv_poly (order,Y,T) / eval_pv_poly (order,Z,T) * 1.0e6 );*/

  f = eval_pv_poly (order,Y,T);
  g = eval_pv_poly (order,Z,T);
  gsq = g*g;
  gcb = gsq*g;

  dfdT = eval_pv_poly_dT (order,Y,T, &df_d_dT);
  dgdT = eval_pv_poly_dT (order,Z,T, &dg_d_dT);

  *dPvsat_dT = (dfdT/g - f/g/g*dgdT) * 1.0e6;


  *d_dPvsat_dT = (df_d_dT/g - 2.0*dfdT*dgdT/gsq + 2.0*f/gcb*dgdT*dgdT - f/gsq*dg_d_dT) * 1.0e6;

  return (f / g * 1.e6);

} 

double eval_pv_poly (int order,
                     double coef[],
                     double T
)
 
     /****************************************************************
      *  Evaluate the polynomial for vapor pressure of water,        *
      *  as defined in FEHM                                          *
      ****************************************************************/
{

   int i;
   double sum, T_power = 1.0;

   sum = coef[0];
   for (i=1; i<=order; i++) {
        T_power = T_power*T;
        sum += coef[i]*T_power;
   }

   return sum;
}
 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double eval_pv_poly_dT (int order,
                        double coef[],
                        double T,
                        double *df_d_dT
)
 
     /****************************************************************
      *  Evaluate the derivative of the polynomial in eval_pv_poly   *
      *  wrt T                                                       *
      ****************************************************************/
{

   int i;
   double sum1, sum2, T_power = 1.0;

   sum1 = coef[1];

   for (i=2; i<=order; i++) {
        T_power = T_power*T;
        sum1 += i*coef[i]*T_power;
   }

   T_power = 1.0;
   sum2 = 2.0*coef[2];
   for (i=3; i<=order; i++) {
        T_power = T_power*T;
        sum2 += i*(i-1)*coef[i]*T_power;
   }
   *df_d_dT = sum2;

   return sum1;
}
 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef NOT_USED
double rho_sat_water_liq_EOS(
          double P_sat,
          double T,
          double *drhol_dP,
          double *drhol_dT,
          double drhol_d_dP[],
          double drhol_d_dT[]
)

/*************************************************************************
 *
 * Equation of state for density of saturated liquid water  (as in FEHM)
 *
 * Units:  Pressure of saturated liquid at T, P_sat - Pascals
 *         temperature, T   -  degrees Celsius
 *
 *         return density in kg/m^3
 *
 *
 *************************************************************************/
  
{

 double f, g, P, dfdP, dfdT, dgdP, dgdT;
 double gsq, gcb;
 double df_d_dT[2], df_d_dP[2], dg_d_dT[2], dg_d_dP[2];

/* The following coefficients are for the pressure range .001 to 20 MPa *
 * and temperatures .5 to 360 degrees C                                 */

  static double Y[] = {0.10000000e+01,
                      -0.50430220e-01,
                       0.12147449e-02,
                      -0.29566543e-04,
                       0.11719555e-01,
                      -0.10272834e-03,
                       0.16483547e-06,
                       0.74802254e-03,
                       0.17552861e-05,
		      -0.16978281e-05}; 
  static double Z[] = {0.10020170e-02,
                      -0.52711077e-04,
                       0.14548166e-05,
                      -0.36472636e-07,
                       0.11718816e-04,
                      -0.93182060e-07,
                       0.12768238e-09,
                       0.72200359e-06,
                       0.18887078e-08,
		      -0.14167944e-08};
                      

/* The following coefficients are for the pressure range .001 to 110 MPa
   and temperatures 15 to 360 degrees C ...                         

  static double Y[] = {0.10000000e+01,
                       0.17472599e-01,
                      -0.20443098e-04,
                      -0.17442012e-06,
                       0.49564109e-02,
                      -0.40757664e-04,
                       0.50676664e-07,
                       0.50330978e-04,
                       0.33914814e-06,
                      -0.18383009e-06}; 
  static double Z[] = {0.10009476e-02,
                       0.16812589e-04,
                      -0.24582622e-07,
                      -0.17014984e-09,
                       0.48841156e-05,
                      -0.32967985e-07,
                       0.28619380e-10,
                       0.53249055e-07,
                       0.30456698e-09,
                      -0.12221899e-09};   ...   end of high-pressure coef   */

  P = P_sat / 1.0e6;  /* function expects pressure in MPa  */

  /*  return ( eval_poly (Y,P_sat,T) 
           / eval_poly (Z,P_sat,T) );
   */

  f = eval_poly (Y,P,T);
  g = eval_poly (Z,P,T);
  gsq = g*g;
  gcb = gsq*g;


  dfdP = eval_poly_dP (Y,P,T, df_d_dP);
  dgdP = eval_poly_dP (Z,P,T, dg_d_dP);

  dfdT = eval_poly_dT (Y,P,T, df_d_dT);
  dgdT = eval_poly_dT (Z,P,T, dg_d_dT);

  *drhol_dP = dfdP/g - f/gsq*dgdP;
  *drhol_dT = dfdT/g - f/gsq*dgdT;

  *drhol_dP = *drhol_dP / 1.0e6;
  *drhol_dT = *drhol_dT;

  drhol_d_dP[0] = (df_d_dP[0]/g - 2.0*dfdP/gsq*dgdP + 2.0*f/gcb*dgdP*dgdP 
                             - f/gsq*dg_d_dP[0]) / 1.0e12;
  drhol_d_dP[1] = (df_d_dP[1]/g - dgdT/gsq + dfdT/gsq*dgdP + 2.0*f/gcb*dgdT*dgdP
                             - f/gsq*dg_d_dP[1]) / 1.0e6;
  drhol_d_dT[0] = (df_d_dT[0]/g - dfdT/gsq*dgdP - dfdP/gsq*dgdT 
                             + 2.0*f/gcb*dgdP*dgdT - f/gsq*dg_d_dT[0]) / 1.0e6;
  drhol_d_dT[1] = df_d_dT[1]/g + 2.0*dfdT/gsq*dgdT + 2.0*f/gcb*dgdT 
                             - f/gsq*dg_d_dT[1];

  return (  f/g  ); 
} 

#endif

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double h_air_EOS(
         double T,
         double *dha_dT,
         double *d_dha_dT
)

/*************************************************************************
 *
 * Equation of state for enthalpy of air (as in FEHM)
 *
 * Units:  temperature, T   -  degrees Celsius
 *
 *         return enthalpy in J/kg
 *************************************************************************/
  
{
  double f, dfdT, d_df_dT;
  static int order = 3;

  static double Y[] = {1003.7,
                       0.025656,
                       4.5457e-4,
		      -2.7107e-7};
  T = T + 273.15;  /* T in K  */
 
  /* return ( eval_pv_poly (order,Y,T) * T ); */

  f    = eval_pv_poly (order,Y,T);
  dfdT = eval_pv_poly_dT (order,Y,T,&d_df_dT);
  *dha_dT = f + dfdT*T;
  *d_dha_dT = dfdT + d_df_dT*T + dfdT;
 
  return ( f * T );
}
  
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double h_water_compressed_EOS(
          double P,
          double T,
          double *dhl_dP,
          double *dhl_dT,
          double d_dhl_dP[2],
          double d_dhl_dT[2]
)

/*************************************************************************
 *
 * Equation of state for enthalpy of liquid water  (as in FEHM)
 *
 * Units:  Pressure, P - Pascals
 *         temperature, T   -  degrees Celsius
 *
 *         return enthalpy in J/kg
 *
 *
 *************************************************************************/
  
{

 double f, g, dfdP, dfdT, dgdP, dgdT;
 double gsq, gcb;
 double df_d_dT[2], df_d_dP[2], dg_d_dT[2], dg_d_dP[2];

/* The following coefficients are for the pressure range .001 to 20 MPa *
 * and temperatures .5 to 360 degrees C                                 */

  static double Y[] = {-0.28892688e-04,
                        0.10155128e-02,
                        0.38182267e-04,
                        0.29406408e-06,
                        0.42068671e-02,
                       -0.26722745e-04,
                        0.39965615e-07,
                        0.14983417e-03,
                        0.11199162e-05,
		       -0.44963038e-06}; 
  static double Z[] = {0.10000000e+01,
                       0.38028489e-01,
                       0.32800006e-03,
                       0.38164755e-07,
                      -0.62817403e-02,
                       0.87410801e-05,
                       0.18991534e-08,
                      -0.11452490e-03,
                      -0.11341777e-06,
		       0.19903338e-08};
                      

/* The following coefficients are for the pressure range .001 to 110 MPa
   and temperatures 15 to 360 degrees C ...                         

  static double Y[] = {0.25623465e-03,
                       0.10184405e-02,
                       0.22554970e-04,
                       0.34836663e-07,
                       0.41769866e-02,
                      -0.21244879e-04,
                       0.25493516e-07,
                       0.89557885e-04,
                       0.10855046e-06,
                      -0.21720560e-06}; 
  static double Z[] = {0.10000000e+01,
                       0.23513278e-01,
                       0.48716386e-04,
                      -0.19935046e-08,
                      -0.50770309e-02,
                       0.57780287e-05,
                       0.90972916e-09,
                      -0.58981537e-04,
                      -0.12990752e-07,
                       0.45872518e-08};   ...   end of high-pressure coef   */

  P = P / 1.0e6;  /* function expects pressure in MPa; 
                     returns enthalpy in MJ/kg - convert to J at return  */

  /*  return ( eval_poly (Y,P,T) 
           / eval_poly (Z,P,T) * 1.0e6 );
  */

  f = eval_poly (Y,P,T);
  g = eval_poly (Z,P,T);
  gsq = g*g;
  gcb = gsq*g;


  dfdP = eval_poly_dP (Y,P,T,df_d_dP);
  dgdP = eval_poly_dP (Z,P,T,dg_d_dP);

  dfdT = eval_poly_dT (Y,P,T,df_d_dT);
  dgdT = eval_poly_dT (Z,P,T,dg_d_dT);

  *dhl_dP = dfdP/g - f/gsq*dgdP;
  *dhl_dT = (dfdT/g - f/gsq*dgdT) * 1.0e6;

  d_dhl_dP[0] = (df_d_dP[0]/g - 2.0*dfdP/gsq*dgdP + 2.0*f/gcb*dgdP*dgdP 
                           - f/gsq*dg_d_dP[0]) / 1.0e6;
  d_dhl_dP[1] = df_d_dP[1]/g - dgdT/gsq + dfdT/gsq*dgdP + 2.0*f/gcb*dgdT*dgdP
                           - f/gsq*dg_d_dP[1];
  d_dhl_dT[0] = (df_d_dT[0]/g - dfdT/gsq*dgdP - dfdP/gsq*dgdT 
                           + 2.0*f/gcb*dgdP*dgdT - f/gsq*dg_d_dT[0]);
  d_dhl_dT[1] = df_d_dT[1]/g + 2.0*dfdT/gsq*dgdT + 2.0*f/gcb*dgdT 
                           - f/gsq*dg_d_dT[1] * 1.0e6;

  return ( f/g * 1.0e6 );
} 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double h_water_superheat_EOS(
          double P,
          double T,
          double *dhwg_dP, 
          double *dhwg_dT,
          double d_dhwg_dP[2],
          double d_dhwg_dT[2]
)

/*************************************************************************
 *
 * Equation of state for enthalpy of  water vapor (as in FEHM)
 *
 * Units:  Pressure, P - Pascals
 *         temperature, T   -  degrees Celsius
 *
 *         return enthalpy in J/kg
 *
 *
 *************************************************************************/
  
{

 double f, g, dfdP, dfdT, dgdP, dgdT;
 double gsq, gcb;
 double df_d_dT[2], df_d_dP[2], dg_d_dT[2], dg_d_dP[2];

/* The following coefficients are for the pressure range .001 to 20 MPa *
 * and temperatures .5 to 360 degrees C                                 */

  static double Y[] = {0.49023415e+00,
                      -0.10000000e+01,
                       0.24474474e-01,
                       0.23476073e-03,
                       0.86459576e-02,
                       0.38256791e-04,
                       0.19190905e-07,
                       0.16237610e-02,
                      -0.74126396e-04,
		      -0.78086106e-06}; 
  static double Z[] = {0.19602927e+00,
                      -0.35954866e+00,
                       0.54884993e-02,
                       0.58496026e-04,
                       0.33114850e-02,
                       0.12829588e-04,
                      -0.20053974e-08,
                       0.78784157e-03,
                      -0.18512345e-04,
		      -0.52896691e-06};
                      

/* The following coefficients are for the pressure range .001 to 20 MPa
   and temperatures 15 to 360 degrees C ...                         

  static double Y[] = {0.31290881e+00,
                      -0.10000000e+01,
                       0.25748596e-01,
                       0.38846142e-03,
                       0.11319298e-01,
                       0.20966376e-04,
                       0.74228083e-08,
                       0.19206133e-02,
                      -0.10372453e-03,
                       0.59104245e-07}; 
  static double Z[] = {0.12511319e+00,
                      -0.36061317e+00,
                       0.58668929e-02,
                       0.99059715e-04,
                       0.44331611e-02,
                       0.50902084e-05,
                      -0.10812602e-08,
                       0.90918809e-03,
                      -0.26960555e-04,
                      -0.36454880e-06};   ...   end of high-pressure coef   */

  P = P/1.0e6;  /* function expects pressure in MPa;
                   returns enthalpy in MJ - convert to J  */

  /* return ( eval_poly (Y,P,T) 
           / eval_poly (Z,P,T) * 1.0e6 );
   */

  f = eval_poly (Y,P,T);
  g = eval_poly (Z,P,T);
  gsq = g*g;
  gcb = gsq*g;


  dfdP = eval_poly_dP (Y,P,T,df_d_dP);
  dgdP = eval_poly_dP (Z,P,T,dg_d_dP);

  dfdT = eval_poly_dT (Y,P,T,df_d_dT);
  dgdT = eval_poly_dT (Z,P,T,dg_d_dT);

  *dhwg_dP = dfdP/g - f/g/g*dgdP;
  *dhwg_dT = dfdT/g - f/g/g*dgdT;

  *dhwg_dP = *dhwg_dP;
  *dhwg_dT = *dhwg_dT * 1.0e6;

  d_dhwg_dP[0] = (df_d_dP[0]/g - 2.0*dfdP/gsq*dgdP + 2.0*f/gcb*dgdP*dgdP 
                           - f/gsq*dg_d_dP[0]) / 1.0e6;
  d_dhwg_dP[1] = df_d_dP[1]/g - dgdT/gsq + dfdT/gsq*dgdP + 2.0*f/gcb*dgdT*dgdP
                           - f/gsq*dg_d_dP[1];
  d_dhwg_dT[0] = (df_d_dT[0]/g - dfdT/gsq*dgdP - dfdP/gsq*dgdT 
                           + 2.0*f/gcb*dgdP*dgdT - f/gsq*dg_d_dT[0]);
  d_dhwg_dT[1] = (df_d_dT[1]/g + 2.0*dfdT/gsq*dgdT + 2.0*f/gcb*dgdT 
                           - f/gsq*dg_d_dT[1]) * 1.0e6;
  return (  f/g * 1.0e6 ); 


} 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double eval_poly ( double coef[],
                       double P,
                       double T
)
 
     /****************************************************************
      *  Evaluate the given polynomial                               *
      ****************************************************************/
{

   int i;
   double sum, T_power = 1.0, P_power = 1.0;
   double PT;

   sum = coef[0];
   for (i=1; i<4; i++) {
        P_power = P_power*P;
        T_power = T_power*T;
        sum += coef[i]*P_power;
        sum += coef[i+3]*T_power;
   }
   PT = P*T;
   sum += coef[7]*PT;
   sum += coef[8]*PT*P;
   sum += coef[9]*PT*T;

   return sum;
}
 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double eval_poly_dP ( double coef[],
                       double P,
                       double T,
                       double d_dP_poly[]
)
 
     /****************************************************************
      *  Evaluate the the derivative of the given polynomial wrt P   *
      *  as well as the second derivatives stored as:                *
      *     d_dP_poly[0] is d/dP[dpoly/dP],                          *
      *     d_dP_poly[1] is d/dT[dpoly/dP]                           *
      ****************************************************************/
{
   double sum, T_power = 1.0, P_power = 1.0;

   sum = coef[1];

   P_power = P_power*P;
   T_power = T_power*T;
   sum += 2.*coef[2]*P_power;
   sum += coef[7]*T_power;

   P_power = P_power*P;
   T_power = T_power*T;
   sum += 3.*coef[3]*P_power;
   sum += coef[9]*T_power;

   sum += 2.*coef[8]*P*T;

   d_dP_poly[0] = 2.0*coef[2] + 6.0*coef[3]*P + 2.0*coef[8]*T;
   d_dP_poly[1] = coef[7] + 2.0*coef[8]*P + 2.0*coef[9]*T;

   return sum;

}
 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double eval_poly_dT ( double coef[],
                       double P,
                       double T,
                       double d_dT_poly[]
)
 
     /****************************************************************
      *  Evaluate the derivative of the given polynomial wrt T       *
      *  as well as the second derivatives stored as:                *
      *     d_dT_poly[0] is d/dP[dpoly/dT],                          *
      *     d_dT_poly[1] is d/dT[dpoly/dT]                           *
      ****************************************************************/
{
   double sum, T_power = 1.0, P_power = 1.0;

   sum = coef[4];

   P_power = P_power*P;
   T_power = T_power*T;
   sum += coef[7]*P_power;
   sum += 2.*coef[5]*T_power;

   P_power = P_power*P;
   T_power = T_power*T;
   sum += coef[8]*P_power;
   sum += 3.*coef[6]*T_power;

   sum += 2*coef[9]*P*T;

   d_dT_poly[0] = coef[7] + 2.0*coef[8]*P + 2.0*coef[9]*T;
   d_dT_poly[1] = 2.0*coef[5] + 6.0*coef[6]*T + 2.0*coef[9]*P;

   return sum;

}
 

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef NOT_USED
static double P_henry_wa(double T)
/*************************************************************************
 *
 * Pressure (Henry's Law)
 *
 * Units:  temperature, T   -  degrees Celsius
 *         pressure,    P   -  Pascals
 *
 *         returns Henry's constant for air, P_henry Pascals
 *
 *************************************************************************/
  
{
   static double P_henry = 6.207e9;  /*  Units in Pascals; constant
                                      *   taken from TOUGH
                                      */
  
  
 
  return  P_henry ;
}  
#endif
/*****************************************************************************/
/*****************************************************************************/


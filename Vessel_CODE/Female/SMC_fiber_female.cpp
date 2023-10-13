/*
A COMPUTATIONAL MODEL PREDICTS SEX-SPECIFIC RESPONSES TO CALCIUM CHANNEL BLOCKERS IN MESENTERIC VASCULAR SMOOTH MUSCLE
Copyright (c) 2023 Gonzalo Hernandez-Hernandez < ghernandezh@ucdavis.edu>
Clancy Lab
University of California, Davis, Ca, USA
Vessel female Model
Usage:	icc -qopenmp SMC_fiber_female.cpp

*Parallel computing
OpenMP is an application programming interface that supports multi-platform shared-memory multiprocessing programming in C++

*/




#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;
//Constants
double F=	96.4853415;
double R=	8.314472;
double T =	295.2;

//Cell Size
double C_m = 16.0;              // Cell membrane capacitance - pF


//Cell Volumes
double VOL_cell = 1.0*1000;     // Cytosolic Volume - 1 pL = 1000 um^3
double Vol_SRcen = 0.05*1000;   // SR central volume - 0.1 pL = 1000 um^3
double Vol_jun = 0.005*1000 ;   // Junctional volume - 0.01 pL = 1000 um^3

//Intracellular concentrations
//All concentrations are in units of millimolar (mM)
double Ca_in;
double Ca_out = 2.0;
double z_Ca = 2.0;

double Na_in;
double Na_out =	130;
double z_Na =	1.0;

double K_in;
double K_out = 5.0;
double z_K = 1.0;

double Cl_in;
double Cl_out =	141;
double z_Cl =	-1.0;

const double BCL = 500000;
const int stimuli = 1;
//const int tl = 20;  // length of the tissue
const double base_dt = 0.005;
const int fold = 1.0/base_dt;
const double dx = 0.001; //Calculated in C
const double Df = 0.00001714*1;//0.0001;
const int tl = 400; // length of the tissue

typedef struct cell {
  double v, v_new, dvdt;
  double dL;
  double dF;
  double p_K  ;
  double p_K15 ;
  double Na_in, Na_in_new;
  double K_in, K_in_new ;
  double Cl_in, Cl_in_new;
  double Ca_in, Ca_in_new;
  double Ca_SRcen,Ca_SRcen_new;
  double Ca_jun, Ca_jun_new;
  double BUF_1,BUF_1_new ;
  double CSQ_SRcen1,CSQ_SRcen1_new ;
  double BUF_jun1,BUF_jun1_new ;
  double R_10, R_10_new ;
  double xa ;
  double xab ;
  double p_cc;

 double I_CaL;
 double I_PMCA;
 double I_NCX;
 double I_BK;
 double I_Ca_leak;
 double I_Na_leak;
 double I_K_leak;
 double I_Cl_leak;

 double I_NSC_Na;
 double I_NSC_K;

 double I_NaK;
 double I_Kv21;
 double I_Kv15;
 double I_Cl;
 double J_SERCA;
 double J_rel;

 double I_ALL_leak;
 double scaling10 ;
} Cell;

typedef struct simState {
    double t, tt;
    int tstep, counter, beat;
    Cell cellData[1][tl];
} SimState;

SimState theState;
SimState *S = &theState;

//All currents are in pA
double Calculate_I_Total(Cell *theCell, double dt, double I_inj);
double Calculate_CaL(Cell *theCell, double dt);
double Calculate_PMCA(Cell *theCell);
double Calculate_I_NCX(Cell *theCell);
double Calculate_I_NaK(Cell *theCell);
double Calculate_BK(Cell *theCell, double dt);
double Calculate_Ca_leak(Cell *theCell);
double Calculate_Na_leak(Cell *theCell);
double Calculate_K_leak(Cell *theCell);
double Calculate_Cl_leak(Cell *theCell);

double Calculate_NSC_Na(Cell *theCell);
double Calculate_NSC_K(Cell *theCell);

double Calculate_I_Kv21(Cell *theCell, double dt);
double Calculate_I_Kv15(Cell *theCell, double dt);
double Calculate_I_Cl(Cell *theCell, double dt);

//fluxes - all units are mM/ms
double Calculate_J_SERCA(Cell *theCell);
//double Calculate_J_RyR();
  //double ITA = [17];
int main () {




    char name[30];
    FILE *output;

    int ll;

    double V;
    Cell *theCell;

    double dt = base_dt;
    double I_Total= 0;
    double dv;
    double cycle_length = BCL;

    double I_inj;
    int done = 0;


    S->counter=0;

    for (ll=0; ll<=(tl-1); ll+=1){
      double ITA[] ={-44.5857,0.000764932,0.764938,0.0296891,0.0145411,5.75111,144.244,42,8.67168e-05,0.105953,0.000163333,0.00726092,0.0350835,0.0133054,0.000264944,0.00319834};


          theCell = &(S->cellData[0][ll]);
          theCell->v_new =ITA[0];
          theCell->dL=ITA[1];
          theCell->dF=ITA[2];
          theCell->p_K =ITA[3];
          theCell->p_K15 = ITA[4];

          theCell->Na_in_new =ITA[5];
          theCell->K_in_new =	ITA[6];
          theCell->Cl_in_new = ITA[7];
          theCell->Ca_in_new = ITA[8];

          theCell->Ca_SRcen_new= ITA[9];
          theCell->Ca_jun_new= ITA[10];

          theCell->BUF_1_new = ITA[11];
          theCell->CSQ_SRcen1_new = ITA[12];
          theCell->BUF_jun1_new = ITA[13];

          theCell->R_10_new= ITA[14];
          //theCell->xa = 4.6E-5;
          theCell->xab = ITA[15];
          theCell->p_cc = ITA[16];

theCell->scaling10 = 80;

         //  if (ll ==0){  theCell->scaling10 = 100;}
         //  if (ll ==1){  theCell->scaling10 = 30;}
         // if(ll ==2){  theCell->scaling10 = 50;}
         //  if (ll ==3){  theCell->scaling10 = 10;}
         //  if (ll ==4){  theCell->scaling10 = 90;}
         //  if (ll ==5){  theCell->scaling10 = 30;}
         //  if (ll ==6){  theCell->scaling10 = 70;}
         //  if (ll ==7){  theCell->scaling10 = 10;}
         //  if (ll ==8){  theCell->scaling10 = 20;}
         //  if (ll ==9){  theCell->scaling10 = 60;}
         //  if (ll ==10){  theCell->scaling10 = 100;}
         //  if (ll ==11){  theCell->scaling10 = 20;}
         //  if(ll ==12){  theCell->scaling10 = 80;}
         //  if (ll ==13){  theCell->scaling10 = 50;}
         //  if (ll ==14){  theCell->scaling10 = 90;}
         //  if (ll ==15){  theCell->scaling10 = 10;}
         //  if (ll ==16){  theCell->scaling10 = 70;}
         //  if (ll ==17){  theCell->scaling10 = 80;}
         //  if (ll ==18){  theCell->scaling10 = 20;}
         //  if (ll ==19){  theCell->scaling10 = 60;}


    }


    S->beat=1;
    S->t = 0.0;
    S->tt = 0.0;
    S->tstep = 0;

    output = fopen("data_FEMALE_FINAL.txt", "w");

    while (!done) {




      double theta, UN,ND, randomN, csc;



#pragma omp parallel for private( ll, theCell)

        for (ll=0; ll<=(tl-1); ll+=1){
      theCell = &(S->cellData[0][ll]);
      theCell->v = theCell->v_new;
      theCell->Na_in  = theCell->Na_in_new ;
      theCell->K_in = theCell->K_in_new ;
      theCell->Cl_in = theCell ->Cl_in_new;
      theCell->Ca_in = theCell->Ca_in_new;
      theCell->Ca_SRcen  = theCell->Ca_SRcen_new;
      theCell->Ca_jun = theCell->Ca_jun_new;
      theCell->BUF_1 = theCell->BUF_1_new ;
      theCell->CSQ_SRcen1  = theCell->CSQ_SRcen1_new ;
      theCell->BUF_jun1 = theCell->BUF_jun1_new ;
      theCell->R_10 = theCell-> R_10_new;
        }

#pragma omp parallel for private( ll, theCell, I_inj, dv, V, I_Total, theta, UN, ND, randomN, csc)


        for (ll=0; ll<=(tl-1); ll+=1){

        theta = drand48();
        UN =  drand48();
        ND = -2*log(UN);
        randomN =  sqrt(ND)*cos(2*3.1416*theta);
        csc = 0.1* randomN*sqrt(dt);

            theCell = &(S->cellData[0][ll]);


// if (S->tt>=10000 && S->tt<=20000)
//   { theCell->scaling10 = 20;}
// else if (S->tt>=30000 && S->tt<=40000)
//   {theCell->scaling10 = 20; }
// else {theCell->scaling10 = 20;}


            V = theCell->v;

            I_Total=Calculate_I_Total(theCell, dt, I_inj);


            //dv=dt*(-I_Total) ;

            if(ll>0 && ll<(tl-1)) { dv=dt*(-I_Total+Df*(S->cellData[0][ll-1].v - 2*theCell->v + S->cellData[0][ll+1].v )/(dx*dx)) + csc; }
            else if (ll==0) { dv=dt*(-I_Total+Df*(-theCell->v + S->cellData[0][ll+1].v )/(dx*dx))+ csc; }
            else if (ll==(tl-1)) { dv=dt*(-I_Total+Df*(-theCell->v + S->cellData[0][ll-1].v )/(dx*dx))+ csc ; }

            theCell->v_new = theCell->v+dv  ;
            theCell->Ca_in;


        } // End for loop


        if(  ( (S->counter) % (fold*1) == 0 ) ) {

            fprintf (output, "%10.f", S->tt);
            for (ll = 0; ll < tl; ll += 1) {
              //  fprintf (output, "\t%8.2f", S->cellData[0][ll].v);
              fprintf (output, "\t%8.10f", S->cellData[0][ll].v);

            }
            fprintf (output, "\n");

      }



        S->t += dt;
        S->tt += dt;
        S->counter += 1;

        if(S->t >= cycle_length) {
          S->t = 0;
            S->beat++;
           cout << "now the beat is: " << S->beat << ", BCL = " << BCL << endl;
    }
      if (S->tt>BCL) { done = 1;}



    }


    return 0;
}



double Calculate_I_Total(Cell *theCell, double dt, double I_inj){
  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);
  double csc = 0.3*randomN*sqrt(dt);

  double I_KvALL, I_ALL_leak, J_rel;

  double I_CaL, I_PMCA, I_NCX, I_BK,I_Ca_leak, I_Na_leak, I_K_leak, I_Cl_leak, I_NaK, I_Kv21, I_Kv15, I_Cl, J_SERCA;

  //Buffers
  //Parameters
  double BUF_T  = 0.3;
  double k_BUF_on  = 22.0 ; //%1/mM*ms
  double k_BUF_off  = 0.077 ; //%1/ms
  double BUF_jun_T  = 0.3; //%mM
  double CSQ_SRcen_T  = 0.3; //mM
  double CSQ_SRper_T  = 0.3; //mM
  double k_CSQ_on  = 2.0 ; //1/mM*ms
  double k_CSQ_off  = 1.60 ;//1/ms

  //J_RyR
  double R_10;
  double Y_RyR_Ca_jun;
  double n_RyR_Ca_jun  = 1.2;
  double Kmr=0.625;
  double Ca_th =0.0010;
  double k_d = 0.1;

  double J_jun_in = 0.0008*(theCell->Ca_jun - theCell->Ca_in);
  theCell-> J_rel = (0.0022*theCell->R_10)*(theCell->Ca_SRcen-theCell->Ca_jun) ;

  theCell->I_CaL = Calculate_CaL(theCell, dt);
  theCell->I_PMCA = Calculate_PMCA(theCell);          //time independent
  theCell->I_NCX = Calculate_I_NCX(theCell);          //time independent
  theCell->I_BK = Calculate_BK(theCell, dt);
  theCell->I_Ca_leak = Calculate_Ca_leak(theCell);    //time independent
  theCell->I_Na_leak = Calculate_Na_leak(theCell);    //time independent
  theCell->I_K_leak = Calculate_K_leak(theCell);      //time independent
  theCell->I_Cl_leak = Calculate_Cl_leak(theCell);    //time independent

  theCell->I_NSC_Na = Calculate_NSC_Na(theCell);      //time independent
  theCell->I_NSC_K  = Calculate_NSC_K(theCell);    //time independent

  theCell->I_NaK = Calculate_I_NaK(theCell);          //time independent
  theCell->I_Kv21 = Calculate_I_Kv21(theCell, dt);
  theCell->I_Kv15 = Calculate_I_Kv15(theCell, dt);
  theCell->I_Cl = Calculate_I_Cl(theCell, dt);

  //Fluxes
  theCell->J_SERCA = Calculate_J_SERCA(theCell);    //time independent
  //J_RyR = Calculate_J_RyR();

  I_KvALL = (theCell->I_Kv21 + theCell->I_Kv15)*1;
  theCell->I_ALL_leak = theCell->I_Ca_leak + theCell->I_Na_leak + theCell->I_K_leak ;

  Y_RyR_Ca_jun = 1.0/(1.0+(pow(0.12/theCell->Ca_SRcen,5)));
  double alpha_inf = (pow(theCell->Ca_jun,n_RyR_Ca_jun )*Y_RyR_Ca_jun )/(pow(theCell->Ca_jun,n_RyR_Ca_jun )*Y_RyR_Ca_jun + pow(0.0650,n_RyR_Ca_jun));
  double tau_alpha = 1.0/(500*(pow(theCell->Ca_jun,n_RyR_Ca_jun )*Y_RyR_Ca_jun + pow(0.0650,n_RyR_Ca_jun)));
  double R_10_dot = (alpha_inf - theCell->R_10)/tau_alpha ;

  double I_Ca_total = (theCell->I_CaL + theCell->I_Ca_leak + theCell->I_PMCA - 2*theCell->I_NCX );
  double I_Na_total = (3*theCell->I_NCX + 3*theCell->I_NaK + theCell->I_NSC_Na + theCell->I_Na_leak);
  double I_K_total = (I_KvALL - 2*theCell->I_NaK + theCell->I_K_leak + theCell->I_BK + + theCell->I_NSC_K);
  double I_Cl_total = (theCell->I_Cl + theCell->I_Cl_leak);

  double ddt_Ca_in = (-I_Ca_total)/(2*F*VOL_cell) - theCell->J_SERCA - (k_BUF_on * theCell->Ca_in * (BUF_T - theCell->BUF_1) - k_BUF_off * theCell->BUF_1) + (1) * J_jun_in;
  double ddt_BUF_1 = (k_BUF_on * theCell->Ca_in * (BUF_T - theCell->BUF_1) - k_BUF_off * theCell->BUF_1);

  double ddt_Ca_SRcen = theCell->J_SERCA * (VOL_cell/Vol_SRcen) - theCell->J_rel * (VOL_cell/Vol_SRcen) -  0*(k_CSQ_on * theCell->Ca_SRcen * (CSQ_SRcen_T - theCell->CSQ_SRcen1) - k_CSQ_off * theCell->CSQ_SRcen1);
  double ddt_CSQ_SRcen1 = (k_CSQ_on * theCell->Ca_SRcen * (CSQ_SRcen_T - theCell->CSQ_SRcen1) - k_CSQ_off * theCell->CSQ_SRcen1);

  double ddt_Ca_jun = theCell->J_rel * (VOL_cell/Vol_jun) - J_jun_in * (VOL_cell/Vol_jun) - 0 * (k_BUF_on * theCell->Ca_jun * (BUF_jun_T - theCell->BUF_jun1)-k_BUF_off + theCell->BUF_jun1);
  double ddt_BUF_jun1 = (k_BUF_on * theCell->Ca_jun * (BUF_jun_T - theCell->BUF_jun1) - k_BUF_off * theCell->BUF_jun1);

  double ddt_Na_in = (-I_Na_total)/(z_Na*F*VOL_cell)*1;
  double ddt_K_in = (-I_K_total)/(z_K*F*VOL_cell)*1;
  double ddt_Cl_in = (I_Cl_total)/(z_Cl*F*VOL_cell)*0;

  //theCell->V= theCell->V + dt * ddt_V;
  theCell->Na_in_new = theCell->Na_in + dt * ddt_Na_in;
  theCell->K_in_new = theCell->K_in + dt * ddt_K_in;
  theCell->Cl_in_new = theCell->Cl_in + dt * ddt_Cl_in;

  theCell->Ca_in_new = theCell->Ca_in + dt * ddt_Ca_in;
  theCell->Ca_SRcen_new = theCell->Ca_SRcen + dt * ddt_Ca_SRcen + 0.00002*randomN ;
  theCell->Ca_jun_new = theCell->Ca_jun + dt * ddt_Ca_jun;

  theCell->BUF_1_new = theCell->BUF_1 + dt * ddt_BUF_1;
  theCell->CSQ_SRcen1_new = theCell->CSQ_SRcen1 + dt * ddt_CSQ_SRcen1;
  theCell->BUF_jun1_new = theCell->BUF_jun1 + dt * ddt_BUF_jun1;

  theCell->R_10_new = theCell->R_10 + dt * R_10_dot   ;

  return (1/C_m) * (theCell->I_CaL + I_KvALL + theCell->I_K_leak + theCell->I_PMCA + theCell->I_NCX + theCell->I_NaK + theCell->I_Ca_leak + theCell->I_BK + theCell->I_Na_leak + theCell->I_NSC_Na + theCell->I_NSC_K  + theCell->I_Cl + theCell->I_Cl_leak );


}

double Calculate_CaL(Cell *theCell, double dt)
{

  double V = theCell->v;
  double J_CaV, dL,dF,tau_dL, tau_dF,dL_bar,dF_bar,alpha_act, beta_act,alpha_ina, beta_ina;
 //alpha_act =0.93321*exp(V/12.9372); //Structural changes
 alpha_act =0.4278*exp(V/10.8837);
 beta_act = 0.1938*exp(V/-11.5189);

 dL_bar = 1/(1+(beta_act/alpha_act));
 tau_dL = 1/(alpha_act + beta_act) +  1.8555;

 alpha_ina =0.0032*exp(V/-44.637);
 beta_ina = 0.0117*exp(V/30.4336)    ;

 dF_bar = 1/(1+(beta_ina/alpha_ina));
 tau_dF = 1/(alpha_ina + beta_ina) +  68.7682;

  theCell->dL = theCell->dL + ((dL_bar-theCell->dL)/tau_dL)*dt;
  theCell->dF = theCell->dF + ((dF_bar-theCell->dF)/tau_dF)*dt;

  if ( V == 0 ) {
       J_CaV = 0.38*1 * theCell->dL*theCell->dF*z_Ca * (F) * ( theCell->Ca_in-Ca_out );
  } else {
      J_CaV =0.38*1 * theCell->dL*theCell->dF* z_Ca *z_Ca *( F) * (F/(R*T)) * V * (( theCell->Ca_in * exp( (z_Ca *F*V)/(R*T)) - Ca_out )/( exp( (z_Ca *F*V)/(R*T) ) - 1 )) ;
    }

  return J_CaV;
}

double Calculate_PMCA(Cell *theCell)
{
  double V = theCell->v;
  //double I_PMCAbar = 4.4460288;
  double I_PMCAbar = 2.8;   // Maximun current (pA)
  double K_mPMCA = 170E-6;            // Half saturation constat in mM

  return I_PMCAbar * pow(theCell->Ca_in,2)/(pow(theCell->Ca_in,2)+pow(K_mPMCA,2));

}

double Calculate_I_NCX(Cell *theCell) {

  double V = theCell->v;
  double phi_R, phi_F, X_NCX;
  double P_NCX = 0.0003;  // pA //Maximum I_NCX
  double gammax = 0.45;

  phi_F= exp(gammax*V*F/(R*T));
  phi_R= exp((gammax-1)*V*F/(R*T));



X_NCX = (pow(theCell->Na_in,3)*Ca_out*phi_F -  pow(Na_out,3)*theCell->Ca_in*phi_R )/(1+0.0003*(pow(Na_out,3)*theCell->Ca_in +pow(theCell->Na_in,3)*Ca_out ));


    return P_NCX*X_NCX;
}

double Calculate_I_NaK(Cell *theCell) {

  double V = theCell->v;
  double N_1, N_2, N_0, N_pow;
  double I_NaK_max = 86.8*1.2;
//  double I_NaK_max = 86.8;         // Maximun current - pA
  double KmNaK_K = 1.6;               // Half saturation concentration for K_out depedency
  double KmNaK_Na = 22;               // Half saturation concentration for Na_in depedency
  double Q10_Nak = 1.87;              // Q10 temperature coefficient

  N_pow = pow(Q10_Nak,((T-309.2)/10));
  N_0 = 1.0/(1 + (0.1245*exp(-0.1*V*F/(R*T))) + (2.19e-3*(exp(Na_out/49.71))*exp(-1.9*V*F/(R*T))) );
  N_1 = pow(K_out,1.1)/(pow(K_out,1.1) + pow(KmNaK_K,1.1));
  N_2 = pow(theCell->Na_in,1.7)/(pow(theCell->Na_in,1.7) + pow(KmNaK_Na,1.7));

  return N_pow * I_NaK_max * N_1 * N_2 * N_0;
}

double Calculate_NSC_Na(Cell *theCell) {
    double V = theCell->v;

    double gns = 0.0123*theCell->scaling10;	 // Maximal conductance for INSCC - nS/pF //10-100 in steps of 10
    double PnsK =	1.3;		     // Permeability of K+
    double PnsNa = 0.9;	       // Permeability of Na+
    double gnsNa =	0.5	;	     // Conductance ratio of Na+ in INSCC
    double gnsK	=0.1;		       // Conductance ratio of K+ in INSCC

    double Enscc = ((R*T/F)*log((PnsK*K_out+PnsNa*Na_out)/(PnsK*theCell->K_in+PnsNa*theCell->Na_in)));

    return (gnsNa) * gns * (V-Enscc);

}

double Calculate_NSC_K(Cell *theCell) {

    double V = theCell->v;
    double gns = 0.0123*theCell->scaling10;	 // Maximal conductance for INSCC - nS/pF //10-100 in steps of 10
    double PnsK =	1.3;		     // Permeability of K+
    double PnsNa = 0.9;	       // Permeability of Na+
    double gnsNa =	0.5	;	     // Conductance ratio of Na+ in INSCC
    double gnsK	=0.1;		       // Conductance ratio of K+ in INSCC
    double Enscc=((R*T/F)*log((PnsK*K_out+PnsNa*Na_out)/(PnsK*theCell->K_in+PnsNa*theCell->Na_in)));

    return (gnsK)*gns*(V - Enscc);

}

double Calculate_Ca_leak(Cell *theCell) {

    double V = theCell->v;
    double G_Ca_b = 0.0031999*0.8;   //Calcium leak conductance [nS] - #####0.003999 in document#####
    double E_Ca = (0.5*R*T/F)*log(Ca_out/theCell->Ca_in);

    return G_Ca_b * (V-E_Ca);

}

double Calculate_Na_leak(Cell *theCell) {
  double V = theCell->v;
  double G_Na_b = 0.001;   //Calcium leak conductance [nS] - #####0.003999 in document#####
  double E_Na = (R*T/F)*log(Na_out/theCell->Na_in);

  return G_Na_b  * (V-E_Na);

}

double Calculate_K_leak(Cell *theCell) {
  double V = theCell->v;
  double G_K_b = 0.009;   //Calcium leak conductance [nS] - #####0.003999 in document#####
  double E_K = (R*T/F)*log(K_out/theCell->K_in);

  return G_K_b* (V-E_K);
}



double Calculate_Cl_leak(Cell *theCell) {  //######

    double V = theCell->v;
    double G_Cl_b =0.0;
    double E_Cl = -R*T/F*log(Cl_out/theCell->Cl_in);

    return G_Cl_b*(V - E_Cl);

}


double Calculate_I_Kv21(Cell *theCell, double dt )
{

  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);

    double V = theCell->v;
    double p_Kbar, tau_p_K, p_K, alpha_x, beta_x;
    double g_K = 2.0;     //Maximun conductance of K2.1 channels - nS

    double E_K = R*T/F*log(K_out/theCell->K_in);
  //  beta_x =0.0368*exp(V/234.6258);
  //  alpha_x =0.1011*exp(V/10.7216);
  //  p_Kbar = 1/(1+(beta_x/alpha_x));
  //  tau_p_K = 1/(alpha_x+beta_x)+11.3301;

  alpha_x =0.0109*exp(V/11.3912);
  beta_x =0.0084*exp(V/241.3324);
  p_Kbar = 1/(1+(beta_x/alpha_x));
  tau_p_K = 1/(alpha_x+beta_x)+46.6311;

    theCell->p_K = theCell->p_K + ((p_Kbar - theCell->p_K)/tau_p_K )*dt + sqrt((dt/tau_p_K)*(theCell->p_K + p_Kbar - 2*theCell->p_K*p_Kbar)/200)*randomN*0;

    return 1*g_K*theCell->p_K*(V-E_K);
}

double Calculate_I_Kv15( Cell *theCell, double dt)
{
  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);

    double V = theCell->v;
    double alpha_15, beta_15,p_Kbar15, tau_p_K15,p_K15;
    double g_K1 = 1.5*1;    //Maximun conductance of K1.5 channels - nS

    double E_K = R*T/F*log(K_out/theCell->K_in);
//   alpha_15 =0.0852*exp(V/13.2254);
//   beta_15 =0.0306*exp(V/540.0872);
//   p_Kbar15 = 1/(1+(beta_15/alpha_15));
//   tau_p_K15 = 1/(alpha_15+beta_15)+10.5378;

alpha_15 =0.0666*exp(V/9.6398);
beta_15 =0.0637*exp(V/123.4370 );
p_Kbar15 = 1/(1+(beta_15/alpha_15));
tau_p_K15 = 1/(alpha_15+beta_15) - 0.0015;


    theCell->p_K15 = theCell->p_K15 + ((p_Kbar15 - theCell->p_K15)/tau_p_K15 )*dt + sqrt((dt/tau_p_K15)*(theCell->p_K15 + p_Kbar15 - 2*theCell->p_K15*p_Kbar15)/200)*randomN*0 ;

    return g_K1*theCell->p_K15*(V-E_K);
}

double Calculate_BK( Cell *theCell, double dt)
{

  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);

    double V = theCell->v;
    double I_BK , J_BK ;
    double BK_T = 8*1;
    double perm_BK = 4.00E-13;

    double gkca, gbka,	xass_z, xass_vh ,xass , xatc , xa, xab, I_BKa , I_BKab;
    double gbkab, xabss_z, xabss_vh,xabss, xabtc;
    double E_K = R*T/F*log(K_out/theCell->K_in);

    xabss_z =(-0.681249/(1.0+((theCell->Ca_jun*1000.0-0.218988)/0.428335)*((theCell->Ca_jun*1000.0-0.218988)/0.428335))+1.40001/(1.0+((theCell->Ca_jun*1000.0+228.71)/684.946)*((theCell->Ca_jun*1000.0+228.71)/684.946))) ;
    xabss_vh =	(8540.23/pow(((theCell->Ca_jun*1000.0+0.401189)/0.00399115),0.668054)) - 109.28;

    xabss =	(1.0/(1.0+exp(-xabss_z*F*(V-xabss_vh)/(R*T))));
    xabtc = (13.8049/(1.0+((V-153.019)/66.4952)*((V-153.019)/66.4952)))	;
    theCell->xab = theCell->xab + ((xabss - theCell->xab)/xabtc)*dt + sqrt((dt/xabtc)*(theCell->xab + xabss - 2*theCell->xab*xabss)/1000)*randomN*0 ;

    if ( V == 0 )
                {   J_BK = 1.E+09 *BK_T * perm_BK * theCell->xab  *(F) * ( theCell->K_in-K_out );}
      else
          {J_BK = 1.E+09 *BK_T * perm_BK * theCell->xab *F*(V*F/(R*T)) * (( theCell->K_in * exp( z_K * (F*V/(R*T)))-K_out ) / ( exp( z_K * (F*V/(R*T)) ) - 1 )) ;}


  //   I_BK = ;

       return  J_BK;
}



double Calculate_I_Cl(Cell *theCell, double dt) {
    double V = theCell->v;
    double P_clca, p_cc_bar, p_cc;
    double tau_pcc = 50;
    double G_clca = 0.80*0;
    double E_Cl = -R*T/F*log(Cl_out/theCell->Cl_in);

    p_cc_bar = pow(theCell->Ca_in,3)/pow(theCell->Ca_in,3) + pow(140,3);
    theCell->p_cc = theCell->p_cc + ((p_cc_bar - theCell->p_cc)/tau_pcc) * dt;

    return G_clca*P_clca*(V-E_Cl);
}



double Calculate_J_SERCA(Cell *theCell) {
  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);


    double V = theCell->v;
    double I_SERCAbar = 0.0000035*1.5;   //Maximun flux via PMCA - mM/ms
    double K_mSERCA = 219E-6;          //Half saturation concentration - mM

    return I_SERCAbar*pow(theCell->Ca_in,2)/(pow(theCell->Ca_in,2) + pow(K_mSERCA,2)) + 0.000001*randomN*0  ;
}

/*
A COMPUTATIONAL MODEL PREDICTS SEX-SPECIFIC RESPONSES TO CALCIUM CHANNEL BLOCKERS IN MESENTERIC VASCULAR SMOOTH MUSCLE
Copyright (c) 2023 Gonzalo Hernandez-Hernandez < ghernandezh@ucdavis.edu>
Clancy Lab
University of California, Davis, Ca, USA
Female Model
Usage:	icc FullCode_Structurex_rand.cpp

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
//Universal Constants
double F = 96.4853415;    // Faraday's constat - C/mmol
double R = 8.314472;      // Ideal gas constant - J/(mol*K)
double T = 295.2;         // Temperature - Kelvin

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

//Main parameter modulating I_NSC currents (I_NSC_K & I_NSC_Na)
//Figure 5, paramater was set to 90
double scaling10 = 90;


const double BCL = 500000;
const int stimuli = 1;
const double base_dt = 0.005;
const int fold = 1.0/base_dt;


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
} Cell;

typedef struct simState {
    double t, tt;
    int tstep, counter, beat;
    Cell cellData[1][1];
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
double Calculate_J_SERCA(Cell *theCell, double dt);

int main () {

  time_t startTime;
	time_t previousTime;
	startTime = time(NULL);
	previousTime = startTime;
  srand48(startTime);

    ofstream myfile("all_ODEs_female.txt");
    ofstream myfile1("Currents_female.txt");


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

//Initial conditions for different membrane potential baselines
//scaling factor  30
//double ITA[] ={-41.0168,0.00145012,0.748035,0.0348826,0.0205698,7.40108,142.593,42,0.00010162,0.0987109,0.000182671,0.00846889,0.032948,0.0147593,0.000237394,0.0043282};
//scaling factor  40
//double ITA[]	={-39.8018,0.00179569,0.71279,0.0441301,0.0227121,8.3065,141.684,42,0.000120774,0.106567,0.000443485,0.0100139,0.0352621,0.032654,0.000894312,0.00802381};
//scaling factor  50
//double ITA[]	={-36.9234,0.00301148,0.72143,0.0434812,0.0293148,9.15314,140.834,42,0.000146498,0.0987891,0.000285037,0.0120513,0.0329691,0.0229765,0.000405987,0.00889779};
//scaling factor  60
//double ITA[]	={-36.6029,0.00317992,0.688147,0.0532145,0.0308449,9.95532,140.028,42,0.000170554,0.10353,0.000487633,0.0139483,0.0343696,0.0354274,0.000910841,0.00831806};
//scaling factor  70
//double ITA[]	={-35.4134,0.00393228,0.678052,0.0571433,0.0345432,10.7366,139.242,42,0.000198181,0.103705,0.000615853,0.0160837,0.03442,0.0429737,0.00121174,0.00811567};
//scaling factor  80
//double ITA[]	={-36.9762,0.0027735,0.706712,0.0489629,0.0221784,11.5106,138.464,42,0.000228372,0.097257,0.00166661,0.0184025,0.0325119,0.117069,0.00318955,0.0675516};
//scaling factor  90
double ITA[]	={-33.0885,0.00596999,0.669711,0.0610689,0.0428124,12.2766,137.693,42,0.000263102,0.102066,0.000794397,0.0209691,0.0339355,0.0544433,0.00155707,0.0107454};



          theCell = &(S->cellData[0][0]);
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
          theCell->xab = ITA[15];
          theCell->p_cc = ITA[16];




    S->beat=1;
    S->t = 0.0;
    S->tt = 0.0;
    S->tstep = 0;

  

    while (!done) {

    double theta = drand48();
    double UN =  drand48();
    double ND = -2*log(UN);
    double randomN =  sqrt(ND)*cos(2*3.1416*theta);
    double csc = 0.2* randomN*sqrt(dt);



            theCell = &(S->cellData[0][0]);
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



            theCell = &(S->cellData[0][0]);

            V = theCell->v;

            I_Total=Calculate_I_Total(theCell, dt, I_inj);

            dv=dt*(-I_Total);

            theCell->v_new = theCell->v+dv + csc;



        if(  ( (S->counter) % (fold*1) == 0 ) ) {


            //Printing variables
            myfile << S->tt << "," << S->cellData[0][0].v << "," << S->cellData[0][0].dL << "," << S->cellData[0][0].dF << "," << S->cellData[0][0].p_K<< "," << S->cellData[0][0].p_K15 << ","  << S->cellData[0][0].Na_in << "," << S->cellData[0][0].K_in <<"," << S->cellData[0][0].Cl_in << "," << S->cellData[0][0].Ca_in <<"," << S->cellData[0][0].Ca_SRcen << "," << S->cellData[0][0].Ca_jun << "," <<  S->cellData[0][0].BUF_1 << "," <<S->cellData[0][0].CSQ_SRcen1 <<"," << S->cellData[0][0].BUF_jun1<< ","  <<  S->cellData[0][0].R_10 << "," << S->cellData[0][0].xab  <<endl;

            myfile1 << S->tt << "\t" <<  S->cellData[0][0].I_CaL  << "\t" << S->cellData[0][0].I_Kv21 << "\t" << S->cellData[0][0].I_Kv15 << "\t" <<S->cellData[0][0].I_Kv21 + S->cellData[0][0].I_Kv15  << "\t" << S->cellData[0][0].I_BK << "\t" <<   S->cellData[0][0].I_NSC_Na  <<"\t" <<   S->cellData[0][0].I_NSC_K  <<"\t" <<   S->cellData[0][0].I_NaK  <<"\t" <<   S->cellData[0][0].I_NCX  << "\t" <<S->cellData[0][0].I_PMCA<< "\t" << S->cellData[0][0].I_Ca_leak<< "\t" << S->cellData[0][0].I_Na_leak  << "\t" << S->cellData[0][0].I_K_leak << "\t" << S->cellData[0][0].I_Cl_leak << "\t" << S->cellData[0][0].I_ALL_leak << "\t" << S->cellData[0][0].J_SERCA<< "\t" << S->cellData[0][0].J_rel <<endl;




      }


        S->tt += dt;
        S->counter += 1;


      if (S->tt>BCL) { done = 1;}

    }

    return 0;
}




double Calculate_I_Total(Cell *theCell, double dt, double I_inj){
  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);
  double csc = 0.2*randomN*sqrt(dt);

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
  theCell-> J_rel = (0.0022*theCell->R_10)*(theCell->Ca_SRcen-theCell->Ca_jun);

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
  theCell->J_SERCA = Calculate_J_SERCA(theCell, dt);    //time independent
  //J_RyR = Calculate_J_RyR();

  I_KvALL = (theCell->I_Kv21 + theCell->I_Kv15);
  theCell->I_ALL_leak = theCell->I_Ca_leak + theCell->I_Na_leak + theCell->I_K_leak ;

  Y_RyR_Ca_jun = 1.0/(1.0+(pow(0.12/theCell->Ca_SRcen,5)));
  double alpha_inf = (pow(theCell->Ca_jun,n_RyR_Ca_jun )*Y_RyR_Ca_jun )/(pow(theCell->Ca_jun,n_RyR_Ca_jun )*Y_RyR_Ca_jun + pow(0.0650,n_RyR_Ca_jun));
  double tau_alpha = 1.0/(500*(pow(theCell->Ca_jun,n_RyR_Ca_jun )*Y_RyR_Ca_jun + pow(0.0650,n_RyR_Ca_jun)));
  double R_10_dot = (alpha_inf - theCell->R_10)/tau_alpha;

  double I_Ca_total = (theCell->I_CaL + theCell->I_Ca_leak + theCell->I_PMCA - 2*theCell->I_NCX);
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
  theCell->Ca_SRcen_new = theCell->Ca_SRcen + dt * ddt_Ca_SRcen + 0.00002*randomN*1 ;
  theCell->Ca_jun_new = theCell->Ca_jun + dt * ddt_Ca_jun;

  theCell->BUF_1_new = theCell->BUF_1 + dt * ddt_BUF_1;
  theCell->CSQ_SRcen1_new = theCell->CSQ_SRcen1 + dt * ddt_CSQ_SRcen1;
  theCell->BUF_jun1_new = theCell->BUF_jun1 + dt * ddt_BUF_jun1;

  theCell->R_10_new = theCell->R_10 + dt * R_10_dot;

  return (1/C_m) * (theCell->I_CaL + I_KvALL + theCell->I_K_leak + theCell->I_PMCA + theCell->I_NCX + theCell->I_NaK + theCell->I_Ca_leak + theCell->I_BK + theCell->I_Na_leak + theCell->I_NSC_Na + theCell->I_NSC_K );

}


double Calculate_CaL(Cell *theCell, double dt)
{

    double V = theCell->v;
    double J_CaV, dL,dF,tau_dL, tau_dF,dL_bar,dF_bar,alpha_act, beta_act,alpha_ina, beta_ina;

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
         J_CaV =1*0.38 *1* theCell->dL*theCell->dF*z_Ca * (F) * ( theCell->Ca_in-Ca_out );
    } else {
        J_CaV = 1*0.38*1*theCell->dL*theCell->dF* z_Ca *z_Ca *( F) * (F/(R*T)) * V * (( theCell->Ca_in * exp( (z_Ca *F*V)/(R*T)) - Ca_out )/( exp( (z_Ca *F*V)/(R*T) ) - 1 )) ;
      }

    return J_CaV;
}


double Calculate_PMCA(Cell *theCell)
{
    double V = theCell->v;
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
    double I_NaK_max = 86.8*1.2;         // Maximun current - pA
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
    double gns = 0.0123*scaling10;	 // Maximal conductance for INSCC - nS/pF //10-100 in steps of 10
    double PnsK =	1.3;		     // Permeability of K+
    double PnsNa = 0.9;	       // Permeability of Na+
    double gnsNa =	0.5	;	     // Conductance ratio of Na+ in INSCC
    double gnsK	=0.1;		       // Conductance ratio of K+ in INSCC

    double Enscc = ((R*T/F)*log((PnsK*K_out+PnsNa*Na_out)/(PnsK*theCell->K_in+PnsNa*theCell->Na_in)));

    return (gnsNa) * gns * (V-Enscc);

}

double Calculate_NSC_K(Cell *theCell) {

    double V = theCell->v;
    double gns = 0.0123*scaling10;	 // Maximal conductance for INSCC - nS/pF //10-100 in steps of 10
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
    double G_Cl_b = 0.0;
    double E_Cl = -R*T/F*log(Cl_out/theCell->Cl_in);

    return G_Cl_b*(V - E_Cl);

}


double Calculate_I_Kv21(Cell *theCell, double dt) {

  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);

    double V = theCell->v;
    double p_Kbar, tau_p_K, p_K, alpha_x, beta_x;
    double g_K = 2.0;     //Maximun conductance of K2.1 channels - nS

    double E_K = R*T/F*log(K_out/theCell->K_in);
    alpha_x =0.0109*exp(V/11.3912);
    beta_x =0.0084*exp(V/241.3324);
    p_Kbar = 1/(1+(beta_x/alpha_x));
    tau_p_K = 1/(alpha_x+beta_x)+46.6311;

    theCell->p_K = theCell->p_K + ((p_Kbar - theCell->p_K)/tau_p_K )*dt;

    return 1*g_K*theCell->p_K*(V-E_K);
}


double Calculate_I_Kv15(Cell *theCell, double dt) {

  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);

  double V = theCell->v;
  double alpha_15, beta_15,p_Kbar15, tau_p_K15,p_K15;
  double g_K1 = 1.5*1;    //Maximun conductance of K1.5 channels - nS

  double E_K = R*T/F*log(K_out/theCell->K_in);


  alpha_15 =0.0666*exp(V/9.6398);
  beta_15 =0.0637*exp(V/123.4370 );
  p_Kbar15 = 1/(1+(beta_15/alpha_15));
  tau_p_K15 = 1/(alpha_15+beta_15) - 0.0015;


  theCell->p_K15 = theCell->p_K15 + ((p_Kbar15 - theCell->p_K15)/tau_p_K15 )*dt;

  return g_K1*theCell->p_K15*(V-E_K);
}


double Calculate_BK(Cell *theCell, double dt) {

  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);

    double V = theCell->v;
    double I_BK , J_BK ;
    double BK_T = 6;
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



double Calculate_J_SERCA(Cell *theCell,double dt) {
  double theta = drand48();
  double UN =  drand48();
  double ND = -2*log(UN);
  double randomN =  sqrt(ND)*cos(2*3.1416*theta);

    double V = theCell->v;
    double I_SERCAbar = 0.0000035*1.5;   //Maximun flux via PMCA - mM/ms
    double K_mSERCA = 219E-6;          //Half saturation concentration - mM

    return I_SERCAbar*pow(theCell->Ca_in,2)/(pow(theCell->Ca_in,2) + pow(K_mSERCA,2)) + 0.00005*randomN*0  ;
}

/*
  ///熱交換器
#include<math.h>
#include<cmath>
#include <iostream>
#include <fstream>
using namespace std;
#include "CNewtonRaphsonMethod.h"
#include"CFluidParameter.h"
#include"PropertyLithiumBromide_ver05.h"
#include"PropertyWater_ver05.h"

PropertyLithiumBromide libr;
PropertyWater wat;

int main()
{

	int Z, G, l, i, j, N, M, O, n, p, t, o, z;
	///Z = 1;/////Tube mesh
	N = 99;/////gamma
	M = 499;/////parameter
	//O = 1;/////Parameter Analysis
	//G = 1;/////Mass flow rate
	//double u[100][100], v[100][100], Ti[100][100], Xi[100][100], Cw[100][100], Et[100][100], Ef[100][100], Ec[100][100], Ed[100][100];
	double K,DAL,w1,w2,wg,h1,h2,Ts,gamma,D1,D2,Ds,As,x1, x12, x2, xi, L1, L2, A1, A2, m1, m2, mi, myul, myug, myus, rhol, rhog, rhos, Sgen, dp_s, CA_b_s, st_s, E_Reb_sur, E_Reb_no_sur, SEC1, SEC2, Delta_st, CA_b, g, st, L, rho, myu, beta, hmin, hminp, Re, Rew, PI, f, gl, k, Reb, qf, mf, dp, Tw, d, hv, ro, ri, cp, Tsat, Te, m, P, Pw, hwi, hwo, U, Gw, Twl, aave, bave, Xave, Tave, uave, a, htcw, htc, mtc, Tin, Xin, B, S, C, Mw, hw0, sw0, Cwi;
	
	//double Gv[100], A[100], b[100], X[100], CA[100], h[100], dr[100], r[100], R[100], dtx[100], dty[100], dxx[100], dxy[100], dux[100], duy[100], dvx[100], dvy[100], Twi[100], htcl[100], mtcl[100];

	ofstream fout("GridTube.csv");
	if (!fout) {
		cout << "ファイルをオープンできませんでした" << endl;
		return 1;
	}
	else
		cout << "ファイルをオープンしました" << endl;

	for(j=1;j<=M;j++){



		///////////////////initialisation risolutore////////////////////
				
		Tsat = 120.41;			//[^C] (decide according to the operative temperature) 
		//fout << Tsat << ";";
		P = wat.t_p(Tsat);	//[kPa] (fixed to be at saturation)
		//fout << P << ";";
		Ts=Tsat+273.15;
		xi=0.1;
		//fout << xi << ";";
		x1=xi;
		//x1=0.9;
		h1=0.0;
		//fout << h1 << ";";
		h2=0.4;
		//fout << h2 << ";";
		L1=1.0;
		//fout << L1 << ";";
		L2=1.0;
		//fout << L2 << ";";
		D1=0.02;//double(j+1);
		//fout << D1 << ";";
		D2=0.02;
		//fout << D2 << ";";
		Ds=D2/D1;
		//fout << Ds << ";";
		A1=3.14*D1*D1/4.0;
		A2=3.14*D2*D2/4.0;
		As=A2/A1;
		mi=0.05;
		//fout << mi << ";";
		S=1.00;
		//fout << S << ";";
		DAL=pow(Ds,19.0/7.0)*pow(L1/L2,4.0/7.0);
		//fout << DAL << ";";

		myul=wat.sat_myul(P);
		myug=wat.sat_myuv(P);
		myus=wat.sat_myul(P)/wat.sat_myuv(P);
		rhol=wat.sat_roul(P);
		rhog=wat.sat_rouv(P);
		rhos=wat.sat_roul(P)/wat.sat_rouv(P);
		fout << myus << ";";
		fout << rhos << ";";

		x1=0.002*j;
		//fout << x1 << ",";
		//fout << 0 << endl;

		//gamma=(1.0/(pow(L1/L2,4.0/7.0)*pow(Ds,5.0/7.0)*As+1.0));//*(0.7+0.1*double(j-1));
		//fout << gamma << ",";
				//fout << 0 << endl;				
		for(i=1;i<=N;i++){
		m1=0.01*i*mi;
		//m2=mi-m1;
		gamma=m1/mi;
		//gamma=1.0/(pow(L1/L2,4.0/7.0)*pow(Ds,5.0/7.0)*As+1.0);
		//gamma=0.36;
		
		x12=x1;
		Sgen=1.0/Ts*mi*mi*(xi/rhog+(1.0-xi)/rhol)*8.0*3.14*gamma*L1*myul/A1/A1/rhol*((1.0+x1*(rhos-1.0))/(1.0+x1*(myus-1.0)));
		Sgen=1.0/Ts*mi*mi*(xi/rhog+(1.0-xi)/rhol)*8.0*3.14*gamma*1.00001*L1*myul/A1/A1/rhol*((1.0+x1*(rhos-1.0))/(1.0+x1*(myus-1.0)));
		K=0.3164*pow(myul,1.0/4.0)/2.0/rhol;
		wg=9.81*rhol/pow(mi,7.0/4.0);
		w1=(L1+h1)*pow(D1,-5.0/4.0)*pow(A1,-7.0/4.0);
		w2=(L2+h2)*pow(D2,-5.0/4.0)*pow(A2,-7.0/4.0);

		//SEC1=1.0*(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)-K*w2*pow(1.0-gamma,7.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)-wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0));
		SEC1=1.0*(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))-(K*w2*pow(1.0-gamma,7.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)+wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0));
		//fout << SEC1 << endl;	
		fout << SEC1 << ",";
		Sgen=mi/Ts*(gamma/rhol*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+(1.0-gamma)/rhol*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0))*(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))*pow(mi,7.0/4.0);
        //fout << Sgen << ",";                       
		////////////////////////////////////////////////////////////////


	/////////////////////////flowrate loop////////////////////////////////////////////////////////////////////////////////////////////////////////////

	

			m = 0.01;			///[kg/s] interations (define the value of total iterations G)
		
		
		////////////////////Experiment conditions////////////////////

		Pw = 101.325;			//cooling water pressure [kPa]

		Xin = 0.62;						//Solution inlet concentration
		Te = libr.sc_T_XTsat(Xin, Tsat);
		Tin = Te;
		Tave = Tin;
		Xave = Xin;
		
		/////////////////////properties calculation//////////////////
		st_s = (libr.sc_st_XT(Xave, Tave))/3.0;       // prova with surfactnat 
		st = libr.sc_st_XT(Xave, Tave);			//////[N/m]   no surfactant
		rho = libr.sc_rho_XT(Xave, Tave);			////[kg/m3]
		myu = libr.sc_visc_XT(Xave, Tave);			///[Pas]	
		gl = 9.81;								//[m/s2]
		//k = libr.sc_thc_XT(Xave, Tave);			////[W/mK]
		//d = libr.sc_d_XT(Xave, Tave);				//////[m2/s]
		//cp = libr.sc_cp_XT(Xave, Tave)*1000.0;		////[J/kgK]
		//a = k / rho / cp;
		//Tsat=wat.p_t( P );
		//Mw = 0.018015;							/////Molar weight [kg/mol]
		//hw0 = -241830;							/////standard water molar enthalpy [J/mol]
		//sw0 = 188.84;								//////standard water molar entropy [J/molK]
		//Cwi = (rho) / Mw * (1.0 - Xin);
		//Re = 4.0*m / L / myu;
		//fout << Re << endl;
		dp_s = (rho*pow(st_s, 3.0)) / (pow(myu, 4.0)*gl);	//dimensionless parameter group for minimum stable thickness calculation from Maron (1982) galileo number
		dp = (rho*pow(st, 3.0)) / (pow(myu, 4.0)*gl);


				
				//Delta_st = 0.00144;  //angolo 89 
				Delta_st = 0.03915;//angolo 29.7
				CA_b = acos(Delta_st / st); //new contact ange calculation 2020/1/10 - no surfactant
				CA_b_s = acos(Delta_st / st_s);// new contact ange calculation 2020/1/10 -  surfactant

			

	//////////////////////////////// calculation/////////////////////////////////////////////////////////////////////////////////////////////////////////////
					
						E_Reb_sur = 23.0;
						E_Reb_no_sur = 83.0;
						//E_Reb_sur = SEC1 * 0.054*(pow(dp, 0.2*SEC2))*pow((log(CA_b) + 3.45), 3);
						//E_Reb_no_sur = SEC1 * 0.054*(pow(dp, 0.2*SEC2))*pow((log(CA_b) + 3.45), 3);
						
							////////////////////NewtonRaphson Method//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

							
                                //fout << DAL << ",";
								//fout << gamma << ",";
								//fout << x1 << ",";
								//fout << x12 << "endl";
								//Sgen=1.0/Ts*mi*mi*(xi/rhog+(1.0-xi)/rhol)*8.0*3.14*gamma*L1*myul/A1/A1/rhol*((1.0+x1*(rhos-1.0))/(1.0+x1*(myus-1.0)));
								//Sgen=mi/Ts*(gamma/rhol*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+(1.0-gamma)/rhol*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0))*(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0));
                                //fout << Sgen << endl;
								}
                                fout << 0 << endl;
								}
                            
}

*/
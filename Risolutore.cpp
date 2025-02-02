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
	M = 1;/////parameter
	//O = 1;/////Parameter Analysis
	//G = 1;/////Mass flow rate
	//double u[100][100], v[100][100], Ti[100][100], Xi[100][100], Cw[100][100], Et[100][100], Ef[100][100], Ec[100][100], Ed[100][100];
	double eps1,eps2,slw,Sl1,Sl2,Sg1,Sg2,Di,K,DAL,w1,w2,wg,h1,h2,Ts,gamma,DP,D1,D2,Ds,As,x1, x12, x2, xi, L1, L2, A1, A2, m1, m2, mi, myul, myug, myus, rhol, rhog, rhos, Sgen, dp_s, CA_b_s, st_s, E_Reb_sur, E_Reb_no_sur, SEC1, SEC2, Delta_st, CA_b, g, st, L, rho, myu, beta, hmin, hminp, Re, Rew, PI, f, gl, k, Reb, qf, mf, dp, Tw, d, hv, ro, ri, cp, Tsat, Te, m, P, Pw, hwi, hwo, U, Gw, Twl, aave, bave, Xave, Tave, uave, a, htcw, htc, mtc, Tin, Xin, B, S, C, Mw, hw0, sw0, Cwi;
	
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
				
		//Tsat = 15.0;			//[^C] (decide according to the operative temperature) 
		Tsat = 120.41;
		fout << Tsat << ";";
		P = wat.t_p(Tsat);	//[kPa] (fixed to be at saturation)
		fout << P << ";";
		Ts=Tsat+273.15;
		xi=0.1;
		fout << xi << ";";
		x1=xi;
		//x1=0.9;
		Di=0.02;
		h1=0.0;
		fout << h1 << ",";
		h2=0.0;
		fout << h2 << ",";
		L1=3.0;
		//L1=0.10*(1.0+3.0*double(j-1));
		fout << L1 << ",";
		L2=3.0;
		//if(j>10){
		//L2=2.0*(0.3+0.3*double(j-1));
		//}
		//else{
		//L2=0.3+0.3*double(j-1);
		//}
		//L2=0.10*(1.0+3.0*double(j-1));
		fout << L2 << ",";
		D1=0.02;//double(j+1);
		fout << D1 << ",";
		D2=0.02;
		fout << D2 << ",";
		Ds=D2/D1;
		fout << Ds << ",";
		A1=3.14*D1*D1/4.0;
		A2=3.14*D2*D2/4.0;
		As=A2/A1;
		//mi=0.880/double(j+10);
		mi=0.059;
		fout << mi << ",";
		
		DAL=pow(Ds,19.0/7.0)*(L1+h1)/(L2+h2);
		fout << DAL << ",";

		eps1=0.001;
		
		eps2=0.001;

		myul=wat.sat_myul(P);
		//myul=0.00015291; //R1234yf
		//myul=0.00012936; //R454c
		//myul=0.00011366; //R32
		//myul=0.000089198; //R410A
		//myul=0.000057314; //CO2
		myug=wat.sat_myuv(P);
		//myug=0.000011425; //R1234yf
		//myug=0.000012383; //R454c
		//myug=0.000012817; //R32
		//myug=0.00001571; //R410A
		//myug=0.000019716; //CO2
		myus=myul/myug;
		rhol=wat.sat_roul(P);
		//rhol=1091.9;//R1234yf
		//rhol=1042.4;//R454c
		//rhol=961.01;//R32
		//rhol=1058.6;//R410A
		//rhol=710.5;//CO2
		rhog=wat.sat_rouv(P);
		//rhog=37.925;//R1234yf
		//rhog=48.547;//R454c
		//rhog=47.339;//R32
		//rhog=64.874;//R410A
		//rhog=242.73;//CO2
		rhos=rhol/rhog;
		fout << myus << ",";
		fout << rhos << ",";
		slw=wat.sat_laml(P);
		//S=1.0;
		S=pow(rhos,1.0/3.0);
		fout << S << ",";
		fout << 0 << endl;


		for(i=1;i<=N;i++){
		m1=0.01*i*mi;
		//m1=mi/2.0+0.005*i*mi;
		m2=mi-m1;
		gamma=m1/mi;
		//gamma=1.0/(pow(L1/L2,4.0/7.0)*pow(Ds,5.0/7.0)*As+1.0)/1.15;
		//gamma=0.36;

		x1=xi;
		//x1=0.99;
		//x1=0.0001;
		x12=x1;
		Sgen=1.0/Ts*mi*mi*(xi/rhog+(1.0-xi)/rhol)*8.0*3.14*gamma*L1*myul/A1/A1/rhol*((1.0+x1*(rhos-1.0))/(1.0+x1*(myus-1.0)));
		//Sgen=1.0/Ts*mi*mi*(xi/rhog+(1.0-xi)/rhol)*8.0*3.14*gamma*1.00001*L1*myul/A1/A1/rhol*((1.0+x1*(rhos-1.0))/(1.0+x1*(myus-1.0)));
		K=0.3164*pow(myul,1.0/4.0)/2.0/rhol;
		wg=9.81*rhol/pow(mi,7.0/4.0);
		w1=(L1+h1)*pow(D1,-5.0/4.0)*pow(A1,-7.0/4.0);
		w2=(L2+h2)*pow(D2,-5.0/4.0)*pow(A2,-7.0/4.0);

		//SEC1=1.0*(K*w1*pow(gamma,7.0/4.0)*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)-K*w2*pow(1.0-gamma,7.0/4.0)*pow(1.0+(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)-wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0));
		//SEC2=pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*myul;

		//SEC1=K*w1*pow(gamma,7.0/4.0)*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)*pow(mi,7.0/4.0);
		//SEC2=wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)*pow(mi,7.0/4.0);
		SEC1=K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)*pow(mi,7.0/4.0);
		SEC2=wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)*pow(mi,7.0/4.0);

		CA_b=pow(mi,7.0/4.0)*(K*w2*pow(1.0-gamma,7.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0));
		CA_b_s=(wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0))*pow(mi,7.0/4.0);
		////////////////////////////////////////////////////////////////
		x2=(xi/(1.0-gamma)-x1/(1.0/gamma-1.0));

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
		//dp_s = (rho*pow(st_s, 3.0)) / (pow(myu, 4.0)*gl);	//dimensionless parameter group for minimum stable thickness calculation from Maron (1982) galileo number
		//dp = (rho*pow(st, 3.0)) / (pow(myu, 4.0)*gl);


				
				//Delta_st = 0.00144;  //angolo 89 
				//Delta_st = 0.03915;//angolo 29.7
				//CA_b = acos(Delta_st / st); //new contact ange calculation 2020/1/10 - no surfactant
				//CA_b_s = acos(Delta_st / st_s);// new contact ange calculation 2020/1/10 -  surfactant

			

	//////////////////////////////// calculation/////////////////////////////////////////////////////////////////////////////////////////////////////////////
					
						//E_Reb_sur = 23.0;
						//E_Reb_no_sur = 83.0;
						//E_Reb_sur = SEC1 * 0.054*(pow(dp, 0.2*SEC2))*pow((log(CA_b) + 3.45), 3);
						//E_Reb_no_sur = SEC1 * 0.054*(pow(dp, 0.2*SEC2))*pow((log(CA_b) + 3.45), 3);
						
							////////////////////NewtonRaphson Method//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

							CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
							mnm.setup(2000, 1 + 1e-12, 1e-6);					//セットアップ　※引数はニュートン法の変数以上にする
							

							p = 0;
							//inizialization ______Semi Empirical Constants______ SEC1,SEC2

							//SEC1 = 1.0; 
							//SEC2 = 1.0;

								mnm.setValue(p, x1);
								p++;
								//mnm.setValue(p, x12);
								//p++;
								//mnm.setValue(p, gamma);
								//p++;
							

							mnm.setAcc(0.1);				//加速度勾配の入力　※なくても可

							mnm.initial();					//計算を開始する前に必ず初期化してください
							for (mnm.main_loop_init(); mnm.main_loop_check(); mnm.main_loop_reinit()) {		// おまじない
								for (mnm.sub_loop_init(); mnm.sub_loop_check(); mnm.sub_loop_reinit()) {	// おまじない
									p = 0;
									//値の設定
									
										x1 = mnm.getValue(p);
										p++;
										//x12 = mnm.getValue(p);
										//p++;
										//gamma = mnm.getValue(p);
										//p++;
									

									p = 0;

								//	mnm.setError(p, E_Reb_sur, SEC1 * 0.054*(pow(dp_s, 0.2*SEC2))*pow((log(CA_b_s) + 3.45), 3));  //old one surfactant
									//p++;
									//mnm.setError(p, E_Reb_no_sur, SEC1 * 0.054*(pow(dp, 0.2*SEC2))*pow((log(CA_b) + 3.45), 3));  //old one no surfactant
									//p++;

										//mnm.setError(p, 0.0, x1*x1+x1*1.0/(rhos-1.0)*1.0/(myus-1.0)*(rhos-1.0-(1.0/gamma-1.0)*(myus-1.0)-1.0/gamma*(myus-1.0-rhos+1.0)*1.0/((L2/L1*A1*A1/A2/A2*1.0/(1.0/gamma-1.0))-1.0))+1.0/(rhos-1.0)*1.0/(myus-1.0)*(1.0/gamma*xi*(myus-1.0-rhos+1.0)*1.0/((L2/L1*A1*A1/A2/A2*1.0/(1.0/gamma-1.0))-1.0)-(1.0/gamma-1.0+1.0/gamma*xi*(rhos-1.0))));
										//mnm.setError(p, 0.0, x1*x1*(1.0-L1*A2*A2/L2/A1/A1/(1.0/gamma-1.0))+x1*(-L1*A2*A2/L2/A1/A1/(rhos-1.0)/(1.0/gamma-1.0)+L1*A2*A2/L2/A1/A1/(myus-1.0)+L1*A2*A2/L2/A1/A1*xi/(1.0-gamma)+1.0/(myus-1.0)-(1.0/gamma-1.0)*1.0/(rhos-1.0)-1.0/gamma*xi)+L1*A2*A2/L2/A1/A1*1.0/((rhos-1.0)*(myus-1.0))*(1.0+xi*(myus-1.0)/(1.0-gamma))+1.0/((rhos-1.0)*(myus-1.0))*(1.0/gamma-1.0+1.0/gamma*xi*(rhos-1.0)));
										//mnm.setError(p, 1.0*(K*w1*pow(gamma,7.0/4.0)*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)), 1.0*(K*w2*pow(1.0-gamma,7.0/4.0)*pow(1.0+(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)+wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)));
										//p++;


									//////////homogeneous +slip///////////
					//mnm.setError(p, 0.0, K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)-K*w2*pow(1.0-gamma,7.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/(gamma)-1.0))-1.0)-wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/(gamma)-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/(gamma)-1.0))-1.0));
										//p++;
										///homogeneous + viscosity///
										//mnm.setError(p, (pow(mi,7.0/4.0)*K*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+rhol*9.81*h1*((1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))))+1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0)))/rhos))*1000.0, (pow(mi,7.0/4.0)*K*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*w2*pow(1.0-gamma,7.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/(gamma)-1.0))-1.0)+rhol*9.81*h2*((1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*S/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 1.0))))+1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*S/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 1.0)))/rhos))*1000.0);
										//p++;

						//mnm.setError(p, 1.0*(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)), 1.0*(K*w2*pow(1.0-gamma,7.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)+wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)));
						//p++;
									////mnm.setError(p, 1.0*(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)+mi*mi/pow(mi,7.0/4.0)*(16.0/(3.14*3.14*Di*Di*Di*Di*pow((rhol*((1.0-1.0/S)*xi-1.0)/((1.0-rhos/S)*xi-1.0)),2.0))*rhol*(((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)-(((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)))-gamma*gamma/(A1*A1*(rhol*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)))+(1.0-2.0*gamma+gamma*gamma)/(A2*A2*rhol*(((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0))))), 1.0*(K*w2*pow(1.0-gamma,7.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)+wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)));
									////p++;

									//////////homogeneous +slip+ viscosity///////////
									//*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)
									//*pow(1.0+(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*(1.0/myus-1.0),1.0/4.0)
								//mnm.setError(p, 1.0*(K*w1*pow(gamma,7.0/4.0)*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)), 1.0*(K*w2*pow(1.0-gamma,7.0/4.0)*pow(1.0+(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)+wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)));
								//p++;	
									
									

									///////// separated model +zivi ///////////
										////32.0*gamma*L1*myug*(myus*x1*((pow(rhos,2.0/3.0)-1.0)*x1+1.0)+rhos*(1-x1)*(1-x1))/(rhol*D1*D1*pow(rhos,2.0/3.0)*x1)
										////32.0*(1.0-gamma)*L2*myug*(myus*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*((pow(rhos,2.0/3.0)-1.0)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))+1.0)+rhos*(1-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*(1-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))))/(rhol*D2*D2*pow(rhos,2.0/3.0)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))
											
								//mnm.setError(p, 4.0/D1/D1/3.14*32.0*gamma*(L1+h1)*myug*(myus*x1*((pow(rhos,2.0/3.0)-1.0)*x1+1.0)+rhos*(1-x1)*(1-x1))/(rhol*D1*D1*pow(rhos,2.0/3.0)*x1)+gl*rhol/mi*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0), 4.0/D2/D2/3.14*32.0*(1.0-gamma)*(L2+h2)*myug*(myus*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*((pow(rhos,2.0/3.0)-1.0)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))+1.0)+rhos*(1-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*(1-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))))/(rhol*D2*D2*pow(rhos,2.0/3.0)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))+gl*rhol*h2/mi*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0));
								//p++;
										
										//Zivi//
										//1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0)))
										//512.0*gamma*myul/(3.14*rhol*pow(D1,4.0))*((1.0-x1)/(pow(1.0-1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))),3.0))+x1*rhos/(myus*pow(1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))),3.0)))
										//gl*rhol*((1.0-1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))))+1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0)))/rhos)


									//mi*512.0*gamma*myul/(3.14*rhol*pow(D1,4.0))*((1.0-x1)/(pow(1.0-1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))),3.0))+x1*rhos/(myus*pow(1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))),3.0)))+gl*h1*rhol*((1.0-1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))))+1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0)))/rhos)
										//(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))
									//1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0)))
									//mi*512.0*(1.0-gamma)*myul/(3.14*rhol*pow(D2,4.0))*((1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/(pow(1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))),3.0))+(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*rhos/(myus*pow(1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))),3.0)))+gl*rhol*h2*((1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))))+1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0)))/rhos)
										

			//mnm.setError(p,mi*512.0/4.0*gamma*myul*(L1+h1)/(3.14*rhol*pow(D1,4.0))*((1.0-x1)/(pow(1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))),3.0))+x1*rhos/(myus*pow(1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))),3.0)))+gl*h1*rhol*((1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))))+1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0)))/rhos) , mi*512.0/4.0*(1.0-gamma)*myul*(L2+h2)/(3.14*rhol*pow(D2,4.0))*((1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/(pow(1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*S/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 1.0))),3.0))+(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*rhos/(myus*pow(1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*S/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 1.0))),3.0)))+gl*rhol*h2*((1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*S/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 1.0))))+1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*S/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 1.0)))/rhos) );
			//p++;
										////mnm.setError(p, 0.0, K*w1*pow(gamma*1.00001,7.0/4.0)*((1.0-rhos/S)*x12-1.0)/((1.0-1.0/S)*x12-1.0)+wg*h1*((1.0-1.0/S)*x12-1.0)/((1.0-rhos/S)*x12-1.0)-K*w2*pow(1.0-gamma*1.00001,7.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma*1.00001)-x12/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma*1.00001)-x12/(1.0/(gamma*1.00001)-1.0))-1.0)-wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x12/(1.0/(gamma*1.00001)-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma*1.00001)-x12/(1.0/(gamma*1.00001)-1.0))-1.0));
										////p++;

										////mnm.setError(p, 0.0, ((mi/Ts*(gamma*1.00001/rhol*((1.0-rhos/S)*x12-1.0)/((1.0-1.0/S)*x12-1.0)+(1.0-gamma*1.00001)/rhol*((1.0-rhos/S)*(xi/(1.0-gamma*1.00001)-x12/(1.0/(gamma*1.00001)-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma*1.00001)-x12/(1.0/(gamma*1.00001)-1.0))-1.0))*(K*w1*pow(gamma*1.00001,7.0/4.0)*((1.0-rhos/S)*x12-1.0)/((1.0-1.0/S)*x12-1.0)+wg*h1*((1.0-1.0/S)*x12-1.0)/((1.0-rhos/S)*x12-1.0)))-(mi/Ts*(gamma/rhol*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+(1.0-gamma)/rhol*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0))*(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))))/(gamma*0.00001));
										////p++;
										
								//////////completely separated Zivi//////////// 
			
				
				//mnm.setError(p, eps1*(L1+h1)*(pow(D2,6.0))*gamma*gamma/(pow(D1,6.0)*eps2*(L2+h2))*((pow(rhos,2.0/3.0)-1.0)*x1+1.0)*((pow(rhos,2.0/3.0)-1.0)*x1+1.0) , (1.0-gamma)*(1.0-gamma)*((pow(rhos,2.0/3.0)-1.0)*(xi-gamma*x1)/(1.0-gamma)+1.0)*((pow(rhos,2.0/3.0)-1.0)*(xi-gamma*x1)/(1.0-gamma)+1.0) );
				//p++;
						//with static head
mnm.setError(p, 32.0*eps1*(L1+h1)*mi*mi*gamma*gamma/(pow(D1,6.0)*3.14*3.14*rhog*pow(rhos,5.0/3.0))*(pow(rhos,2.0/3.0)+1.0)*pow(((pow(rhos,2.0/3.0)-1.0)*x1+1.0),2.0)+gl*h1*rhol*((1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))))+1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0)))/rhos) , 32.0*eps2*(L2+h2)*mi*mi*(1.0-gamma)*(1.0-gamma)/(pow(D2,6.0)*3.14*3.14*rhog*pow(rhos,5.0/3.0))*(pow(rhos,2.0/3.0)+1.0)*pow(((pow(rhos,2.0/3.0)-1.0)*((xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))+1.0),2.0)+gl*rhol*h2*((1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*S/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 1.0))))+1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))*S/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 1.0)))/rhos));
p++;
				
				
									// mnm.prt();				//エラー表示
									mnm.prt_sum();			//エラーの合計を表示

								}
								
							}

								//////////completely separated Zivi////////////
							//K=eps1*(L1+h1)*(pow(D2,6.0))/(pow(D1,6.0)*eps2*(L2+h2));
                           // if (K == 1){
							//	x1=xi;
							//}
							///else{
							//x1=-((eps1*(L1+h1)*(pow(D2,6.0))/(pow(D1,6.0)*eps2*(L2+h2))-1.0)*gamma-(pow(eps1*(L1+h1)*(pow(D2,6.0))/(pow(D1,6.0)*eps2*(L2+h2)),0.5)-1.0)*(pow(rhos,2.0/3.0)-1.0)*xi-pow(eps1*(L1+h1)*(pow(D2,6.0))/(pow(D1,6.0)*eps2*(L2+h2)),0.5)+1)/((eps1*(L1+h1)*(pow(D2,6.0))/(pow(D1,6.0)*eps2*(L2+h2)-1.0)*(pow(rhos,2.0/3.0)-1.0)*gamma));	
							//}	
								//////////
							fout << DAL << ",";
								fout << gamma << ",";
								fout << x1 << ",";
								//fout << x12 << "endl";
								//Sgen=1.0/Ts*mi*mi*(xi/rhog+(1.0-xi)/rhol)*8.0*3.14*gamma*L1*myul/A1/A1/rhol*((1.0+x1*(rhos-1.0))/(1.0+x1*(myus-1.0)));
								//Sgen=mi/Ts*(gamma/rhol*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+(1.0-gamma)/rhol*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0))*(K*w1*pow(gamma,7.0/4.0)*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0));
								//SEC1=(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)));
								//fout << SEC1 << ",";
								
								///////////homogeneous +slip//////
						//Sgen=mi/Ts*(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)))*(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))*pow(mi,7.0/4.0);
								//////////homogeneous +viscosity////
								//Sgen=mi/Ts*(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)))*(K*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))))+1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0)))/rhos))*pow(mi,7.0/4.0);
                        
						//SEC1=(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)));        
						//DP=(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))*pow(mi,7.0/4.0);
								//DP=(K*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))))+1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0)))/rhos))*pow(mi,7.0/4.0);
								

								///////////homogeneous +slip + viscosity//////
								//*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)
								//Sgen=mi/Ts*(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)))*(K*w1*pow(gamma,7.0/4.0)*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))*pow(mi,7.0/4.0);
								//SEC1=(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)));        
								//DP=(K*w1*pow(gamma,7.0/4.0)*pow(1.0+x1*(1.0/myus-1.0),1.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))*pow(mi,7.0/4.0);
									
								//////////separated +zivi////
								/////Sgen=mi/Ts*mi*(xi/rhog+(1.0-xi)/rhol)*(32.0*gamma*L1*myug*(myus*x1*((pow(rhos,2.0/3.0)-1.0)*x1+1.0)+rhos*(1-x1)*(1-x1))/(rhol*D1*D1*pow(rhos,2.0/3.0)*x1));
                        //Sgen=mi/Ts*(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)))*(4.0/D1/D1/3.14*32.0*gamma*(L1+h1)*myug*(myus*x1*((pow(rhos,2.0/3.0)-1.0)*x1+1.0)+rhos*(1-x1)*(1-x1))/(rhol*D1*D1*pow(rhos,2.0/3.0)*x1)+gl/mi*rhol*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0));
                       // DP=4.0/D1/D1/3.14*32.0*gamma*(L1+h1)*myug*(myus*x1*((pow(rhos,2.0/3.0)-1.0)*x1+1.0)+rhos*(1-x1)*(1-x1))/(rhol*D1*D1*pow(rhos,2.0/3.0)*x1)+gl*rhol/mi*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0);
						//SEC1=(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)));						
							       
						// //Sl1=mi*mi*gamma*(1.0-x1)/rhol/Ts*32.0*myul/rhol*L1*gamma/D1/D1*(x1*(pow(rhos,2.0/3.0)-1.0)+1.0);
								////Sg1=mi*mi*gamma*x1/rhog/Ts*32.0*myug/rhog*L1*gamma/D1/D1*(1.0/(pow(rhos,2.0/3.0))+x1*(1.0-1.0/(pow(rhos,2.0/3.0))));
								////Sl2=mi*mi*(1.0-gamma)*(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/rhol/Ts*32.0*myul/rhol*L2*gamma/D2/D2*((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*(pow(rhos,2.0/3.0)-1.0)+1.0);
								////Sg2=mi*mi*(1.0-gamma)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))/rhog/Ts*32.0*myug/rhog*L2*gamma/D2/D2*(1.0/(pow(rhos,2.0/3.0))+(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*(1.0-1.0/(pow(rhos,2.0/3.0))));
								////Sgen=Sl1+Sl2+Sg1+Sg2;

								//Zivi//
								///Homo//
		//Sgen=mi/Ts*(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)))*(mi*512.0/4.0*gamma*myul*(L1+h1)/(3.14*rhol*pow(D1,4.0))*((1.0-x1)/(pow(1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))),3.0))+x1*rhos/(myus*pow(1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))),3.0)))+gl*h1*rhol*((1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))))+1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0)))/rhos));
		//DP=(mi*512.0/4.0*gamma*myul*(L1+h1)/(3.14*rhol*pow(D1,4.0))*((1.0-x1)/(pow(1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))),3.0))+x1*rhos/(myus*pow(1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))),3.0)))+gl*h1*rhol*((1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))))+1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0)))/rhos));
		//SEC1=(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)));						
								////Fully separated//
					//Sgen=mi*mi/Ts*512.0/4.0*(gamma*gamma*(L1+h1)/(3.14*pow(D1,4.0))*((myul/rhol/rhol*(pow(1.0-x1,2.0))/(pow(1.0-1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))),3.0)))+myug/rhog*rhog*((pow(x1,2.0))/(pow(1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))),3.0))))+(1.0-gamma)*(1.0-gamma)*((L2+h2)/(3.14*pow(D2,4.0))*((myul/rhol/rhol*(pow(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)),2.0))/(pow(1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))),3.0)))+myug/rhog*rhog*((pow((xi/(1.0-gamma)-x1/(1.0/gamma-1.0)),2.0))/(pow(1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))),3.0))))))+mi/Ts*(gl*h1*rhol*((1.0-1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))))+1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0)))/rhos)*(gamma/(rhol*(1.0-1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0))))+rhog*1.0/(1.0+(1.0-x1)/(x1*pow(rhos, 2.0/3.0)))))+gl*rhol*h2*((1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))))+1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0)))/rhos)*(gamma/(rhol*(1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))))+rhog*1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))))));
								
						
						
						////Fully separated Zivi almost analyitical//
						//Sgen=128.0*mi*mi*mi*gamma*gamma*gamma*eps1*(L1+h1)/(Ts*3.14*3.14*pow(D1,6.0)*rhog*rhog*rhos*rhos)*(1.0+x1*(pow(rhos,2.0/3.0)-1.0))*(x1*((pow(rhos,4.0/3.0)-pow(rhos,2.0/3.0)+1.0)*x1+pow(rhos,2.0/3.0)-2.0)+1.0)+128.0*mi*mi*mi*(1.0-gamma)*(1.0-gamma)*(1.0-gamma)*eps2*(L2+h2)/(Ts*3.14*3.14*pow(D2,6.0)*rhog*rhog*rhos*rhos)*(1.0+(xi-gamma*x1)/(1.0-gamma)*(pow(rhos,2.0/3.0)-1.0))*((xi-gamma*x1)/(1.0-gamma)*((pow(rhos,4.0/3.0)-pow(rhos,2.0/3.0)+1.0)*(xi-gamma*x1)/(1.0-gamma)+pow(rhos,2.0/3.0)-2.0)+1.0);
                      //  DP=32.0*eps1*(L1+h1)*mi*mi*gamma*gamma/(pow(D1,6.0)*3.14*3.14*rhog*pow(rhos,5.0/3.0))*(pow(rhos,2.0/3.0)+1.0)*pow(((pow(rhos,2.0/3.0)-1.0)*x1+1.0),2.0);
						
						//with static head
Sgen=mi/Ts*(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)))*(32.0*eps1*(L1+h1)*mi*mi*gamma*gamma/(pow(D1,6.0)*3.14*3.14*rhog*pow(rhos,5.0/3.0))*(pow(rhos,2.0/3.0)+1.0)*pow(((pow(rhos,2.0/3.0)-1.0)*x1+1.0),2.0)+gl*h1*rhol*((1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))))+1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0)))/rhos));
DP=32.0*eps1*(L1+h1)*mi*mi*gamma*gamma/(pow(D1,6.0)*3.14*3.14*rhog*pow(rhos,5.0/3.0))*(pow(rhos,2.0/3.0)+1.0)*pow(((pow(rhos,2.0/3.0)-1.0)*x1+1.0),2.0)+gl*h1*rhol*((1.0-1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0))))+1.0/(1.0+(1.0-x1)*S/(x1*pow(rhos, 1.0)))/rhos);
								//(L2+h2)*((myul/rhol/rhol*(pow(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)),2.0))/(pow(1.0-1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))),3.0)))+myug/rhog*rhog*((pow((xi/(1.0-gamma)-x1/(1.0/gamma-1.0)),2.0))/(pow(1.0/(1.0+(1.0-(xi/(1.0-gamma)-x1/(1.0/gamma-1.0)))/((xi/(1.0-gamma)-x1/(1.0/gamma-1.0))*pow(rhos, 2.0/3.0))),3.0))))
								//momentum//
								//Sgen=mi/Ts*(gamma/rhol*(((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))+(1.0-gamma)/rhol*(((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)))*((K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))*pow(mi,7.0/4.0)+mi*mi*(16.0/(3.14*3.14*Di*Di*Di*Di*pow((rhol*((1.0-1.0/S)*xi-1.0)/((1.0-rhos/S)*xi-1.0)),2.0))*rhol*(((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))-gamma*gamma/(A1*A1*(rhol*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)))));
                                
								//SEC2=mi*mi*(1.0/(3.14*3.14*Di*Di*Di*Di*(rhol*((1.0-1.0/S)*xi-1.0)/((1.0-rhos/S)*xi-1.0)))*rhol*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)-gamma*gamma/(A1*A1*(rhol*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0)))+1.0/(3.14*3.14*Di*Di*Di*Di*(rhol*((1.0-1.0/S)*xi-1.0)/((1.0-rhos/S)*xi-1.0)))*rhol*(((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0))-(1.0-2.0*gamma+gamma*gamma)/(A2*A2*rhol*(((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0))));
								//fout << SEC2 << ",";
								//Sgen=mi/Ts*(gamma/rhol*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0))*(K*w1*pow(gamma,7.0/4.0)*((1.0-rhos/S)*x1-1.0)/((1.0-1.0/S)*x1-1.0)+wg*h1*((1.0-1.0/S)*x1-1.0)/((1.0-rhos/S)*x1-1.0))*pow(mi,7.0/4.0)+mi/Ts*((1.0-gamma)/rhol*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0))*((K*w2*pow(1.0-gamma,7.0/4.0)*((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)+wg*h2*((1.0-1.0/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)/((1.0-rhos/S)*(xi/(1.0-gamma)-x1/(1.0/gamma-1.0))-1.0)))*pow(mi,7.0/4.0);
								//fout << SEC1 << ",";
								fout << DP << ",";
								fout << Sgen << endl;
								}
                                fout << 0 << endl;
								}
                            
}


*/
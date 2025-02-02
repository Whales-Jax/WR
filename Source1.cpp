/*
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

#include "CFluidParameter.h"
#include "CNewtonRaphsonMethod.h"
#include "DXProperty_ver06.h"
#include "parameterN.h"

int main()
{
	ofstream complex("HeatExchanger(results).csv");				//Results printout

	DXProperty ref;
	ref.reibai = "R32.FLD", "AIR.PPF";
	ref.joutai = "DEF";
	ref.ref_state[0] = "DEF";
	ref.ref_state[1] = "IIR";
	ref.LoadDLL("refprop.DLL");
	ref.setmulti();

	double dt;																//time step [s]

	//===Refrigerant westmost boundary Properties at previous time step===//(from/to other device)
	double G_RM_in_dt;														//refrigerant inlet mass flowrate [kg/s]
	double P_RM_in_dt;														//refrigerant inlet pressure [kPa]
	double h_RM_in_dt;														//refrigerant inlet specific enthalpy [kJ/kg]
	double T_RM_in_dt;														//refrigerant inlet temperature [degC]
	double u_RM_in_dt;														//refrigerant inlet specific internal energy [kJ/kg]
	double rho_RM_in_dt;													//refrigerant inlet density [kg/m^3]
	double myu_RM_in_dt;													//refrigerant inlet viscosity [Pa-s]
	double k_RM_in_dt;														//refrigerant inlet thermal conductivity [W/m-k]
	double cp_RM_in_dt;														//refrigerant inlet specific heat [J/kg-K]
	double x_RM_in_dt;														//refrigerant inlet quality

	//===Air frontal boundary Properties at previous time step===//(from environment)
	double G_AM_in_dt;														//air mass flowrate [kg/s]
	double P_AM_in_dt;														//air inlet pressure [kPa]
	double T_AM_in_dt;														//air inlet Temperature [degC]

	//===Refrigerant westmost boundary Properties current time step===//(from/to other device)
	double G_RM_in;															//refrigerant inlet mass flowrate [kg/s]
	double P_RM_in;															//refrigerant inlet pressure [kPa]
	double h_RM_in;															//refrigerant inlet specific enthalpy [kJ/kg]
	double T_RM_in;															//refrigerant inlet temperature [degC]
	double u_RM_in;															//refrigerant inlet specific internal energy [kJ/kg]
	double rho_RM_in;														//refrigerant inlet density [kg/m^3]
	double myu_RM_in;														//refrigerant inlet viscosity [Pa-s]
	double k_RM_in;															//refrigerant inlet thermal conductivity [W/m-k]
	double cp_RM_in;														//refrigerant inlet specific heat [J/kg-K]
	double x_RM_in;															//refrigerant inlet quality


	//===Air frontal boundary Properties current time step===//(from environment)
	double G_AM_in;															//air inlet mass flowrate [kg/s]
	double P_AM_in;															//air inlet pressure [kPa]
	double T_AM_in;															//air inlet Temperature [degC]
	double v_AM_in;															//air inlet Temperature [degC]

	//===Refrigerant eastmost boundary Properties current time step===//(to/from other device)
	double G_RM_out;														//refrigerant mass flowrate
	double P_RM_out;														//refrigerant inlet pressure
	double h_RM_out;														//refrigerant inlet specific enthalpy

	//===Heat Exchanger Geometry===//
	const double tube_D_o = 0.01;											//tube outer diameter [m]
	const double tube_D_i = 0.009;											//tube inner diameter [m]
	const double tube_L = 1.0;												//tube length [m]
	const double tube_T = 0.0005;											//tube thickness [m]
	const double tube_adjH = 0.03;											//horizontal distance between adjacent tubes [m]
	const double tube_adjV = 0.03;											//vertical distance between adjacent tubes [m]
	const double fin_S = 0.003;												//fin spacing [m]
	const double fin_T = 0.00013;											//fin thickness [m]
	const double PI = (6 * (asin(0.5)));									//pi
	const double beta = 0.0;												//tube inclination [rad]
	const double ktube = 205.0 / 1000.0;									//aluminum thermal conductivity [kW/m-K]

	//Other refrigerant-related parameters
	//double eta_fin[50], Re_r[50], Prandtl_r[50];

	//Other air-related parameters
	//double Re_a[50], Re_a_max[50], Prandtl_a[50];
	double D_H_air, Vmax_a, Sd, Vmesh;
	double jcf, jcf_t;														//Colburn factor
	double fin_ff, fin_ft, fin_dP_f, fin_dP_t, fin_dP[50], fin_f;			//air-side pressure drop variables

	//Other parameters for Newton-Raphson
	int i, j, k, p;															//Counter for Newton-Rahpson variables
	double alfa_r[50], alfa_a[50], KA[50], Q[50], f[50];
	double Rf, Rc, DP, DP_i, DP_o, DT, DH, Qtot, Gtot, Atot;
	double H_res, G_res;
	double SA_tot, SA_finned, Ac_fin;

	//Variables for the TTCM

	int Nrow, Ncol, Nmesh, Ntubes;
	int xrow, xcol, yrow, ycol;
	int checker[50][50];													//Used for checking the U-bend connectivity

	int size;																//Number of elements of the array that simplifies the TTCM

	int c_pdrop;															//Counter of the "longest" possible path
	double pdrop;															//Initialized pressure drop per tube

	SA_finned = (int)(tube_L / fin_S) * 2 * tube_adjH * tube_adjV;
	Ac_fin = (tube_L - (((int)(tube_L / fin_S))*fin_T)) * tube_adjV;
	SA_tot = SA_finned + ((tube_L - (((int)(tube_L / fin_S))* fin_T)) * PI * tube_D_o);

	Nmesh = 10;
	Ncol = 1;
	Nrow = 1;
	Ntubes = Nrow * Ncol;

	parameter::refr **R = new parameter::refr*[Nmesh];
	parameter::air **A= new parameter::air*[Nmesh];
	parameter::prop **P = new parameter::prop*[Nmesh];
	for (int i = 0; i < Nmesh; ++i) {
		parameter::refr* refri = new parameter::refr[Ncol];
		R[i] = refri;
		parameter::air* airp = new parameter::air[Ncol];
		A[i] = airp;
		parameter::prop* prope = new parameter::prop[Ncol];
		P[i] = prope;
	}
	//===Input Conditions===//
	G_RM_in = 0.015;														//refrigerant mass flowrate
	P_RM_in = 1500.0;														//refrigerant inlet pressure
	h_RM_in = 250.0;														//refrigerant inlet specific enthalpy
	G_AM_in = 0.27;															//refrigerant mass flowrate
	P_AM_in = 101.325;														//refrigerant inlet pressure
	T_AM_in = 35.0;															//refrigerant inlet specific enthalpy

	clock_t start_time, end_time;
	start_time = clock();

	Vmesh = PI * tube_D_o*tube_D_o / 4.0*tube_L / double(Nmesh);					//refrigerant mesh volume [m^3]

	ref.change_pure(1);
	ref.state_ph(P_RM_in, h_RM_in);
	T_RM_in = ref.Rc.T;														//refrigerant inlet temperature

	ref.change_pure(2);
	ref.state_tp(T_AM_in, P_AM_in);
	double rho_AM_in = ref.Rc.rho;											//air inlet density [kg/m^3]

	v_AM_in = G_AM_in / ((tube_L - (int)(tube_L / fin_S) * fin_T) * tube_adjV *double(Nrow)* rho_AM_in);		//air inlet velocity

	Sd = pow(((pow(tube_adjH, 0.5)) + (pow((tube_adjV / 2), 0.5))), 0.5);

	if (Sd < ((tube_adjV + tube_D_o) / 2))
	{
		Vmax_a = (tube_adjV / (2 * (Sd - tube_D_o))) * v_AM_in;
	}

	else
	{
		Vmax_a = (tube_adjV / (tube_adjV - tube_D_o)) * v_AM_in;
	}

	//===Initialization of Local and boundary Variables at previous time step===//
	//??//


	//===Initialization of local and boundary Variables at current time step===//

	R[0][1].G_R_B = G_RM_in;
	R[0][1].h_R_B = h_RM_in;
	R[0][1].P_R_B = P_RM_in;

	for (int j = 1; j <= Ncol; j++)
	{
		for (int i = 1; i <= Nmesh; i++)
		{
			R[i][j].G_R_B = R[i][j].G_R_B_dt;
			R[i][j].h_R_B = R[i][j].h_R_B_dt;
			R[i][j].P_R_B = R[i][j].P_R_B_dt;
			R[i][j].h_R_P = R[i][j].h_R_P_dt;
			R[i][j].P_R_P = R[i][j].P_R_P_dt;
			R[i][j].T_R_P = R[i][j].T_R_P_dt;
			R[i][j].u_R_P = R[i][j].u_R_P_dt;
		}
	}


	for (int i = 0; i <= Nmesh; i++)
	{
		A[i][0].G_A_B = G_AM_in / double(Nmesh) / double(Nrow);
		A[i][0].T_A_B = T_AM_in;
		A[i][0].P_A_B = P_AM_in;

		for (int j = 1; j <= Ncol; j++)
		{
			//Air-side local and boundary parameters
			A[i][j].T_A_B = A[i][j].T_A_B_dt;
			A[i][j].P_A_B = A[i][j].P_A_B_dt;
			A[i][j].T_A_P = A[i][j].T_A_P_dt;
			A[i][j].P_A_P = A[i][j].P_A_P_dt;
			A[i][j].G_A_B = A[i][j].G_A_B_dt;
		}
	}



	//===Newton-Rahpson Calculation===//
	CNewtonRaphsonMethod mnm;
	mnm.setup(2000, 1 + 1e-12, 1e-4);

	p = 0;
	for (int j = 1; j <= Ncol; j++)
	{
		for (int i = 1; i <= Nmesh; i++)
		{
			mnm.setValue(p, R[i][j].G_R_B);
			p++;
			mnm.setValue(p, R[i][j].h_R_P);
			p++;
			mnm.setValue(p, R[i][j].P_R_P);
			p++;

		}
	}

	mnm.setAcc(1.0);
	mnm.initial();

	for (mnm.main_loop_init(); mnm.main_loop_check(); mnm.main_loop_reinit())
	{
		for (mnm.sub_loop_init(); mnm.sub_loop_check(); mnm.sub_loop_reinit())
		{
			p = 0;

			for (int j = 1; j <= Ncol; j++)
			{
				for (int i = 1; i <= Nmesh; i++)
				{
					R[i][j].G_R_B = mnm.getValue(p);
					p++;
					R[i][j].h_R_P = mnm.getValue(p);
					p++;
					R[i][j].P_R_P = mnm.getValue(p);
					p++;

				}
			}
			p = 0;



			//Fundamental equations, mesh loop

			for (int j = 1; j <= Ncol; j++)
			{
				for (int i = 1; i <= Nmesh; i++)
				{

					//===Condition of windup integration===//
					if (R[i][j].G_R_B >= 0.0)
					{

						ref.change_pure(1);
						ref.state_ph(R[i][j].P_R_P, R[i][j].h_R_P);
						R[i][j].x_R_P = ref.Rc.x;
						R[i][j].T_R_P = ref.Rc.T;
						R[i][j].u_R_P = ref.Rc.u;

						//===Single phase properties===//
						if (R[i][j].x_R_P < 0.0 || R[i][j].x_R_P > 1.0)
						{

							R[i][j].myu_R_P = ref.Rc.visc / 1000000.0;
							R[i][j].k_R_P = ref.Rc.thc;
							R[i][j].cp_R_P = ref.Rc.cp * 1000.0;
							R[i][j].rho_R_P = ref.Rc.rho;

						}

						//===Two phase properties===//
						else {

							ref.change_pure(1);
							ref.sat_p(R[i][j].P_R_P);
							R[i][j].myu_Rl_P = ref.Rl.visc / 1000000.0;					//in Pa-s
							R[i][j].myu_Rv_P = ref.Rv.visc / 1000000.0;
							R[i][j].cp_Rl_P = ref.Rl.cp * 1000.0;							//in J/kg-K
							R[i][j].cp_Rv_P = ref.Rv.cp * 1000.0;
							R[i][j].k_Rl_P = ref.Rl.thc;									//in W/m-k
							R[i][j].k_Rv_P = ref.Rv.thc;
							R[i][j].rho_Rl_P = ref.Rl.rho;								//in kg/m^3
							R[i][j].rho_Rv_P = ref.Rv.rho;
							R[i][j].h_Rl_P = ref.Rl.h;									//in kJ/kg
							R[i][j].h_Rv_P = ref.Rv.h;
							R[i][j].st_Rlv_P = ref.Rl.st;									//in N/m

							R[i][j].myu_R_P = (R[i][j].x_R_P * R[i][j].myu_Rv_P) + ((1.0 - R[i][j].x_R_P) * R[i][j].myu_Rl_P);
							R[i][j].VF_R_P = 1.0 / (1.0 + (((1.0 - R[i][j].x_R_P) / R[i][j].x_R_P) * (R[i][j].rho_Rv_P / R[i][j].rho_Rl_P)));
							R[i][j].rho_R_P = (R[i][j].VF_R_P * R[i][j].rho_Rv_P) + ((1.0 - R[i][j].VF_R_P) * R[i][j].rho_Rl_P);

						}

						R[i][j].Re_R_P = R[i][j].G_R_B * 4.0 / PI / tube_D_i / R[i][j].myu_R_P;
						R[i][j].Pr_R_P = R[i][j].myu_R_P * R[i][j].cp_R_P / R[i][j].k_R_P;

						///===Air-side properties===/// (...A_P[i]=...A_B[i][j] inlet air state for calculation of local air properties)
						ref.change_pure(2);
						ref.state_tp(A[i][j-1].T_A_B, A[i][j - 1].P_A_B);

						A[i][j].rho_A_P = ref.Rc.rho;
						A[i][j].myu_A_P = ref.Rc.visc / 1000000.0;
						A[i][j].k_A_P = ref.Rc.thc;
						A[i][j].cp_A_P = ref.Rc.cp * 1000.0;
						A[i][j].T_A_P = A[i][j-1].T_A_B;
						A[i][j].P_A_P = A[i][j - 1].P_A_B;

						///===Air-side continuity===///
						A[i][j].G_A_B = A[i][j-1].G_A_B;

						A[i][j].v_A_B = A[i][j].G_A_B / ((tube_L / double(Nmesh) - (int)(tube_L / fin_S) * fin_T) * tube_adjV * A[i][j].rho_A_P);		//air inlet velocity
						A[i][j].Re_A_P = A[i][j].rho_A_P * A[i][j-1].v_A_B * tube_D_o / A[i][j].myu_A_P;
						A[i][j].Remax_A_P = A[i][j].rho_A_P * Vmax_a * tube_D_o / A[i][j].myu_A_P;
						A[i][j].Pr_A_P = A[i][j].myu_A_P * A[i][j].cp_A_P / A[i][j].k_A_P;

						///===simplified overall heat transfer coefficient===///
						P[i][j].K_P = 2.0 * PI * tube_D_o * tube_L / double(Nmesh);

						///===heat transfer parameters air-side===///

						fin_ft = (4.0 / PI) * (0.25 + (0.118 * pow(A[i][j].Remax_A_P, -0.16) / (pow((tube_adjV / tube_D_o) - 1.0, 1.08)))) * ((tube_adjV / tube_D_o) - 1);
						fin_ff = 1.455 * pow(A[i][j].Remax_A_P, -0.656) * pow((tube_adjV / tube_adjH), -0.347) * pow((fin_S / tube_D_o), -0.134) * pow((tube_adjV / tube_D_o), 1.23);
						fin_f = (fin_ff * SA_finned * SA_tot) + (fin_ft * (1.0 - (SA_finned / SA_tot)) * (1.0 - (fin_T / SA_finned)));


						///===Pressure drop===///
						A[i][j].P_A_B = A[i][j - 1].P_A_B - (fin_f * SA_tot * pow((A[i][j].rho_A_P * Vmax_a), 2) / (Ac_fin * 2 * A[i][j].rho_A_P)) / 1000;

						///===Air outlet temperature===///
						A[i][j].T_A_B = ((P[i][j].K_P - A[i][j].G_A_B * A[i][j].cp_A_P) * A[i][j-1].T_A_B - P[i][j].K_P * R[i][j].T_R_P) / (P[i][j].K_P + A[i][j].G_A_B * A[i][j].cp_A_P);

						///===Overall heat transfer===///
						P[i][j].q_P = P[i][j].K_P * ((A[i][j-1].T_A_B + A[i][j].T_A_B) / 2.0 - R[i][j].T_R_P);


						///===Continuity error===///
						R[i][j].M_R_P = Vmesh * R[i][j].rho_R_P;

						mnm.setError(p, (R[i][j].M_R_P - R[i][j].M_R_P_dt) / dt + R[i-1][j].G_R_B - R[i][j].G_R_B, 0.0);
						p++;

						///===Energy conservation error===/// (h_R_B[i] = h_R_P[i])
						mnm.setError(p, (R[i][j].M_R_P * R[i][j].u_R_P - R[i][j].M_R_P_dt * R[i][j].u_R_P_dt) / dt + R[i - 1][j].G_R_B * R[i - 1][j].h_R_P - R[i][j].G_R_B * R[i][j].h_R_P + P[i][j].q_P, 0.0);
						p++;

						if (i == 1) {

							///===Pressure drop error===/// Westmost boundary (P_R_B[i] = P_R_P[i])
							mnm.setError(p, R[i - 1][j].P_R_B - R[i][j].P_R_P - 128.0 * tube_L / double(Nmesh) / 2.0 * R[i][j].G_R_B * R[i][j].myu_R_P / (PI * R[i][j].rho_R_P * pow(tube_D_i, 4.0)), 0.0);
							p++;

						}

						else {

							///Eastmost boundary///
							if (i == Nmesh) {

								mnm.setError(p, R[i][j].P_R_P - R[i][j].P_R_B - 128.0 * tube_L / double(Nmesh) / 2.0 * R[i][j].G_R_B * R[i][j].myu_R_P / (PI * R[i][j].rho_R_P * pow(tube_D_i, 4.0)), 0.0);
								p++;

							}

							else {

								mnm.setError(p, R[i - 1][j].P_R_P - R[i][j].P_R_P - 128.0 * tube_L / double(Nmesh) * R[i][j].G_R_B * R[i][j].myu_R_P / (PI * R[i][j].rho_R_P * pow(tube_D_i, 4.0)), 0.0);
								p++;

							}
						}


					}

					else {

						ref.change_pure(1);
						ref.state_ph(R[i][j].P_R_P, R[i][j].h_R_P);
						R[i][j].x_R_P = ref.Rc.x;
						R[i][j].T_R_P = ref.Rc.T;
						R[i][j].u_R_P = ref.Rc.u;

						if (R[i][j].x_R_P < 0.0 || R[i][j].x_R_P > 1.0)
						{
							R[i][j].myu_R_P = ref.Rc.visc / 1000000.0;
							R[i][j].k_R_P = ref.Rc.thc;
							R[i][j].cp_R_P = ref.Rc.cp * 1000.0;
							R[i][j].rho_R_P = ref.Rc.rho;

						}

						else {
							ref.change_pure(1);
							ref.sat_p(R[i][j].P_R_P);
							R[i][j].myu_Rl_P = ref.Rl.visc / 1000000.0;					//in Pa-s
							R[i][j].myu_Rv_P = ref.Rv.visc / 1000000.0;
							R[i][j].cp_Rl_P = ref.Rl.cp * 1000.0;							//in J/kg-K
							R[i][j].cp_Rv_P = ref.Rv.cp * 1000.0;
							R[i][j].k_Rl_P = ref.Rl.thc;									//in W/m-k
							R[i][j].k_Rv_P = ref.Rv.thc;
							R[i][j].rho_Rl_P = ref.Rl.rho;								//in kg/m^3
							R[i][j].rho_Rv_P = ref.Rv.rho;
							R[i][j].h_Rl_P = ref.Rl.h;									//in kJ/kg
							R[i][j].h_Rv_P = ref.Rv.h;
							R[i][j].st_Rlv_P = ref.Rl.st;								//in N/m
							R[i][j].myu_R_P = (R[i][j].x_R_P * R[i][j].myu_Rv_P) + ((1.0 - R[i][j].x_R_P) * R[i][j].myu_Rl_P);
							R[i][j].VF_R_P = 1.0 / (1.0 + (((1.0 - R[i][j].x_R_P) / R[i][j].x_R_P) * (R[i][j].rho_Rv_P / R[i][j].rho_Rl_P)));
							R[i][j].rho_R_P = (R[i][j].VF_R_P * R[i][j].rho_Rv_P) + ((1.0 - R[i][j].VF_R_P) * R[i][j].rho_Rl_P);
						}

						R[i][j].Re_R_P = R[i][j].G_R_B * 4.0 / PI / tube_D_i / R[i][j].myu_R_P;
						R[i][j].Pr_R_P = R[i][j].myu_R_P * R[i][j].cp_R_P / R[i][j].k_R_P;

						///===Air-side properties===/// (...A_P[i]=...A_B[i][j] inlet air state for calculation of local air properties)
						ref.change_pure(2);
						ref.state_tp(A[i][j-1].T_A_B, A[i][j - 1].P_A_B);

						A[i][j].rho_A_P = ref.Rc.rho;
						A[i][j].myu_A_P = ref.Rc.visc / 1000000.0;
						A[i][j].k_A_P = ref.Rc.thc;
						A[i][j].cp_A_P = ref.Rc.cp * 1000.0;
						A[i][j].T_A_P = A[i][j-1].T_A_B;
						A[i][j].P_A_P = A[i][j - 1].P_A_B;

						A[i][j].G_A_B = A[i][j-1].G_A_B;

						A[i][j].v_A_B = A[i][j].G_A_B / ((tube_L / double(Nmesh) - (int)(tube_L / fin_S) * fin_T) * tube_adjV * A[i][j].rho_A_P);		//air inlet velocity
						A[i][j].Re_A_P = A[i][j].rho_A_P * A[i][j-1].v_A_B * tube_D_o / A[i][j].myu_A_P;
						A[i][j].Remax_A_P = A[i][j].rho_A_P * Vmax_a * tube_D_o / A[i][j].myu_A_P;
						A[i][j].Pr_A_P = A[i][j].myu_A_P * A[i][j].cp_A_P / A[i][j].k_A_P;

						///===simplified overall heat transfer coefficient===///
						P[i][j].K_P = 2.0 * PI * tube_D_o * tube_L / double(Nmesh);

						///===heat transfer air-side===///

						fin_ft = (4.0 / PI) * (0.25 + (0.118 * pow(A[i][j].Remax_A_P, -0.16) / (pow((tube_adjV / tube_D_o) - 1.0, 1.08)))) * ((tube_adjV / tube_D_o) - 1);
						fin_ff = 1.455 * pow(A[i][j].Remax_A_P, -0.656) * pow((tube_adjV / tube_adjH), -0.347) * pow((fin_S / tube_D_o), -0.134) * pow((tube_adjV / tube_D_o), 1.23);
						fin_f = (fin_ff * SA_finned * SA_tot) + (fin_ft * (1.0 - (SA_finned / SA_tot)) * (1.0 - (fin_T / SA_finned)));



						A[i][j].P_A_B = A[i][j - 1].P_A_B - (fin_f * SA_tot * pow((A[i][j].rho_A_P * Vmax_a), 2.0) / (Ac_fin * 2 * A[i][j].rho_A_P)) / 1000.0;

						A[i][j].T_A_B = ((P[i][j].K_P - A[i][j].G_A_B * A[i][j].cp_A_P) * A[i][j-1].T_A_B - P[i][j].K_P * R[i][j].T_R_P) / (P[i][j].K_P + A[i][j].G_A_B * A[i][j].cp_A_P);

						P[i][j].q_P = P[i][j].K_P * ((A[i][j-1].T_A_B + A[i][j].T_A_B) / 2.0 - R[i][j].T_R_P);

						///===Continuity error===///
						R[i][j].M_R_P = Vmesh * R[i][j].rho_R_P;

						mnm.setError(p, (R[i][j].M_R_P - R[i][j].M_R_P_dt) / dt + R[i][j].G_R_B - R[i+1][j].G_R_B, 0.0);
						p++;

						///===Energy conservation error===/// (h_R_B[i] = h_R_P[i])
						mnm.setError(p, (R[i][j].M_R_P * R[i][j].u_R_P - R[i][j].M_R_P_dt * R[i][j].u_R_P_dt) / dt + R[i][j].G_R_B * R[i][j].h_R_P - R[i+1][j].G_R_B * R[i+1][j].h_R_P + P[i][j].q_P, 0.0);
						p++;

						if (i == 1) {

							///===Pressure drop error===/// Westmost boundary (P_R_B[i] = P_R_P[i])
							mnm.setError(p, R[i - 1][j].P_R_B - R[i][j].P_R_P + 128.0 * tube_L / double(Nmesh) / 2.0 * R[i][j].G_R_B * R[i][j].myu_R_P / (PI * R[i][j].rho_R_P * pow(tube_D_i, 4.0)), 0.0);
							p++;

						}

						else {

							///Eastmost boundary///
							if (i == Nmesh) {

								mnm.setError(p, R[i][j].P_R_P - R[i][j].P_R_B + 128.0 * tube_L / double(Nmesh) / 2.0 * R[i][j].G_R_B * R[i][j].myu_R_P / (PI * R[i][j].rho_R_P * pow(tube_D_i, 4.0)), 0.0);
								p++;

							}

							else {

								mnm.setError(p, R[i-1][j].P_R_P - R[i][j].P_R_P + 128.0 * tube_L / double(Nmesh) * R[i][j].G_R_B * R[i][j].myu_R_P / (PI * R[i][j].rho_R_P * pow(tube_D_i, 4.0)), 0.0);
								p++;

							}
						}



					}
				}
			}



			mnm.prt_sum();
		}


		///===Pressure and enthalpy windup===/// 
		for (int j = 1; j <= Ncol; j++)
		{
			for (int i = 1; i <= Nmesh; i++)
			{
				if (R[i][j].G_R_B >= 0.0)
				{

					R[i][j].h_R_B = R[i][j].h_R_P;
					R[i][j].P_R_B = R[i][j].P_R_P;

				}

				else {

					R[i][j].h_R_B = R[i+1][j].h_R_P;
					R[i][j].P_R_B = R[i+1][j].P_R_P;

				}
				A[i][j].T_A_B = A[i][j].T_A_P;
				A[i][j].P_A_B = A[i][j].P_A_P;

			}
		}

	}
	//End of Newton-Rahpson loop//

	//===Calculation of Final Outputs//
	DP = 0.0;
	Qtot = 0.0;
	Gtot = 0.0;
	DP_i = 0.0;
	DP_o = 0.0;


	cout << "\nRESULTS";
	for (int j = 1; j <= Ncol; j++)
	{
		for (int i = 1; i <= Nmesh; i++)
		{
			Qtot += P[i][j].q_P;

			cout << "\n\nAt tube_" << i << "\trefrigerant's\tP_ = " << R[i][j].P_R_P << "\t\th_ = " << R[i][j].h_R_P << "\t\tx_ = " << R[i][j].x_R_P << "\t\tT_ = " << R[i][j].T_R_P << "\t\tG_ = " << R[i][j].G_R_B;
			cout << "\n\t\tair's\t\tP_ = " << A[i][j].P_A_P << "\t\tT_ = " << A[i][j].T_A_B;
		}
		cout << endl;

		Atot = (tube_D_o * PI * (tube_L - (int)(tube_L / fin_S) * fin_T) + (int)(tube_L / fin_S) * 2.0 * (tube_adjV * tube_adjH - PI / 4.0 * tube_D_o * tube_D_o))* double(Ntubes);
		cout << "\nQtot = " << Qtot << " kW\tDP = " << R[0][1].P_R_B - R[Nmesh][1].P_R_B << " kPa ";

		end_time = clock();
		double run_time = (double)(end_time - start_time) / 1000;			//[sec]
		cout << "\n\nTotal run time: " << run_time << " s";

		//Writing results on .csv file
		cout << "\n\nWriting to .csv file...";

		complex << "Qtot = " << Qtot << "," << " DP = " << DP_i + DP_o << "," << " Gtot = " << Gtot << "," << " Atot = " << Atot << ",";


		//===Initialization of local and boundary Variables at current time step===//

		G_RM_in_dt = G_RM_in;
		h_RM_in_dt = h_RM_in;
		P_RM_in_dt = P_RM_in;

		for (int j = 1; j <= Ncol; j++)
		{
			for (int i = 1; i <= Nmesh; i++)
			{
				R[i][j].G_R_B_dt = R[i][j].G_R_B;
				R[i][j].h_R_B_dt = R[i][j].h_R_B;
				R[i][j].P_R_B_dt = R[i][j].P_R_B;
				R[i][j].h_R_P_dt = R[i][j].h_R_P;
				R[i][j].P_R_P_dt = R[i][j].P_R_P;
				R[i][j].T_R_P_dt = R[i][j].T_R_P;
				R[i][j].u_R_P_dt = R[i][j].u_R_P;
				R[i][j].M_R_P_dt = R[i][j].M_R_P;
				P[i][j].q_P_dt = P[i][j].q_P;
			}
		}


		for (int i = 0; i <= Nmesh; i++)
		{
			G_AM_in_dt = G_AM_in;
			T_AM_in_dt = T_AM_in;
			P_AM_in_dt = P_AM_in;

			for (int j = 1; j <= Ncol; j++)
			{
				//Air-side local and boundary parameters
				A[i][j].T_A_B_dt = A[i][j].T_A_B;
				A[i][j].P_A_B_dt = A[i][j].P_A_B;
				A[i][j].T_A_P_dt = A[i][j].T_A_P;
				A[i][j].P_A_P_dt = A[i][j].P_A_P;
				A[i][j].G_A_B_dt = A[i][j].G_A_B;
			}
		}





		///1
		int i = 0;
		for (int c = 0; c < Ncol; c++)
		{
			for (int r = 0; r < Nrow; r++)
			{
				++i;
				//T_in_r_mtrx[r][c] = T_in_r[i];
			}
		}

		complex << "\nRefrigerant Inlet Temperature\n";
		for (int r = 0; r < Nrow; r++)
		{
			for (int c = 0; c < Ncol; c++)
			{
				//complex << T_in_r_mtrx[r][c] << ",";
			}
			complex << endl;
		}

		///2
		i = 0;
		for (int c = 0; c < Ncol; c++)
		{
			for (int r = 0; r < Nrow; r++)
			{
				++i;
				//T_out_r_mtrx[r][c] = T_out_r[i];
			}
		}

		complex << "\nRefrigerant Outlet Temperature\n";
		for (int r = 0; r < Nrow; r++)
		{
			for (int c = 0; c < Ncol; c++)
			{
				//complex << T_out_r_mtrx[r][c] << ",";
			}
			complex << endl;
		}

		///3
		i = 0;
		for (int c = 0; c < Ncol; c++)
		{
			for (int r = 0; r < Nrow; r++)
			{
				++i;
				//P_in_r_mtrx[r][c] = P_in_r[i];
			}
		}

		complex << "\nRefrigerant Inlet Pressure\n";
		for (int r = 0; r < Nrow; r++)
		{
			for (int c = 0; c < Ncol; c++)
			{
				//complex << P_in_r_mtrx[r][c] << ",";
			}
			complex << endl;
		}

		///4
		i = 0;
		for (int c = 0; c < Ncol; c++)
		{
			for (int r = 0; r < Nrow; r++)
			{
				++i;
				//P_out_r_mtrx[r][c] = P_out_r[i];
			}
		}

		complex << "\nRefrigerant Outlet Pressure\n";
		for (int r = 0; r < Nrow; r++)
		{
			for (int c = 0; c < Ncol; c++)
			{
				//complex << P_out_r_mtrx[r][c] << ",";
			}
			complex << endl;
		}

		///5
		i = 0;
		for (int c = 0; c < Ncol; c++)
		{
			for (int r = 0; r < Nrow; r++)
			{
				++i;
				//T_in_air_mtrx[r][c] = T_in_air[i];
			}
		}

		complex << "\nAir Inlet Temperature\n";
		for (int r = 0; r < Nrow; r++)
		{
			for (int c = 0; c < Ncol; c++)
			{
				//complex << T_in_air_mtrx[r][c] << ",";
			}
			complex << endl;
		}

		///6
		i = 0;
		for (int c = 0; c < Ncol; c++)
		{
			for (int r = 0; r < Nrow; r++)
			{
				++i;
				//T_out_air_mtrx[r][c] = T_out_air[i];
			}
		}

		complex << "\nAir Outlet Temperature\n";
		for (int r = 0; r < Nrow; r++)
		{
			for (int c = 0; c < Ncol; c++)
			{
				//complex << T_out_air_mtrx[r][c] << ",";
			}
			complex << endl;
		}

		///7
		i = 0;
		for (int c = 0; c < Ncol; c++)
		{
			for (int r = 0; r < Nrow; r++)
			{
				++i;
				//h_in_r_mtrx[r][c] = h_in_r[i];
			}
		}

		complex << "\nRefrigerant Inlet Specific Enthalpy\n";
		for (int r = 0; r < Nrow; r++)
		{
			for (int c = 0; c < Ncol; c++)
			{
				//complex << h_in_r_mtrx[r][c] << ",";
			}
			complex << endl;
		}

		///8
		i = 0;
		for (int c = 0; c < Ncol; c++)
		{
			for (int r = 0; r < Nrow; r++)
			{
				++i;
				//h_out_r_mtrx[r][c] = h_out_r[i];
			}
		}

		complex << "\nRefrigerant Outlet Specific Enthalpy\n";
		for (int r = 0; r < Nrow; r++)
		{
			for (int c = 0; c < Ncol; c++)
			{
				//complex << h_out_r_mtrx[r][c] << ",";
			}
			complex << endl;
		}

		///9
		i = 0;
		for (int c = 0; c < Ncol; c++)
		{
			for (int r = 0; r < Nrow; r++)
			{
				++i;
				//G_r_mtrx[r][c] = G_r[i];
			}
		}

		complex << "\nRefrigerant Flowrate\n";
		for (int r = 0; r < Nrow; r++)
		{
			for (int c = 0; c < Ncol; c++)
			{
				//complex << G_r_mtrx[r][c] << ",";
			}
			complex << endl;
		}

		complex << "\nDisplaying the Tube-Tube Connectivity Matrix\n";
		complex << ",";
		for (int i = 1; i <= Ntubes; i++)
		{
			complex << "t(" << i << ")" << ",";
		}

		for (int i = 1; i <= Ntubes; i++)
		{
			complex << "\nj(" << i << ")" << ",";
			for (int j = 1; j <= Ntubes; j++)
			{
				//complex << tube_mtrx[i][j] << ",";
			}
		}

		cout << "\n\n";
		system("pause");
		return 0;
	}

	}*/
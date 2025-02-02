/*
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

#include "CFluidParameter.h"
#include "CNewtonRaphsonMethod.h"
#include "DXProperty_ver06.h"

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
	
	//===Refrigerant Local Properties at previous time step===//
	double T_R_P_dt[50][50];												//refrigerant local temperature [degC]
	double P_R_P_dt[50][50];												//refrigerant local pressure [kPa]
	double h_R_P_dt[50][50];												//refrigerant local specific enthalpy [kJ/kg]
	double u_R_P_dt[50][50];												//refrigerant local specific internal energy [kJ/kg]
	double M_R_P_dt[50][50];
//	double rho_R_P_dt[50][50];												//refrigerant local density [kg/m^3]
//	double myu_R_P_dt[50][50];												//refrigerant local viscosity [Pa-s]
//	double k_R_P_dt[50][50];												//refrigerant local thermal conductivity [W/m-k]
//	double cp_R_P_dt[50][50];												//refrigerant local specific heat [J/kg-K]
//	double x_R_P_dt[50][50];												//refrigerant local quality

//	double Re_R_P_dt[50][50];													//
//	double Pr_R_P_dt[50][50];													//

	
	//===Refrigerant Boundary Properties at previous time step===//
	double G_R_B_dt[50][50];												//refrigerant boundary mass flowrate [kg/s]
	double P_R_B_dt[50][50];												//refrigerant boundary pressure [kPa]
	double h_R_B_dt[50][50];												//refrigerant boundary specific enthalpy [kJ/kg]
	double T_R_B_dt[50][50];												//refrigerant boundary temperature [degC]
//	double rho_R_B_dt[50][50];												//refrigerant boundary density [kg/m^3]
//	double myu_R_B_dt[50][50];												//refrigerant boundary viscosity [Pa-s]
//	double k_R_B_dt[50][50];												//refrigerant boundary thermal conductivity [W/m-k]
//	double cp_R_B_dt[50][50];												//refrigerant boundary specific heat [J/kg-K]
//	double x_R_B_dt[50][50];												//refrigerant boundary quality

//	double Re_R_B_dt[50][50];													//
//	double Pr_R_B_dt[50][50];													//

	//===Air Local Properties at previous time step===//
	double P_A_P_dt[50][50];												//air local ipressure [kPa]
	double h_A_P_dt[50][50];												//air local specific enthalpy [kJ/kg]
	double T_A_P_dt[50][50];												//air local temperature [degC]
	double u_A_P_dt[50][50];												//air local specific internal energy [kJ/kg]
	double v_A_P_dt[50][50];												//air local velocity [m/s]
//	double rho_A_P_dt[50][50];												//air local density [kg/m^3]
//	double myu_A_P_dt[50][50];												//air local viscosity [Pa-s]
//	double k_A_P_dt[50][50];												//air local thermal conductivity [W/m-k]
//	double cp_A_P_dt[50][50];												//air local specific heat [J/kg-K]
//	double x_A_P_dt[50][50];												//air local absolute humidity

//	double Re_A_P_dt[50][50];												//
//	double Pr_A_P_dt[50][50];												//


	//===Air Boundary Properties at previous time step===//
	double G_A_B_dt[50][50];												//air boundary mass flowrate [kg/s]
	double h_A_B_dt[50][50];												//air boundary specific enthalpy [kJ/kg]
	double P_A_B_dt[50][50];												//air boundary pressure [kPa]
	double T_A_B_dt[50][50];												//air boundary temperature [degC]
	double v_A_B_dt[50][50];												//air boundary velocity [m/s]

//	double rho_A_B_dt[50][50];												//air boundary density [kg/m^3]
//	double myu_A_B_dt[50][50];												//air boundary viscosity [Pa-s]
//	double k_A_B_dt[50][50];												//air boundary thermal conductivity [W/m-k]
//	double cp_A_B_dt[50][50];												//air boundary specific heat [J/kg-K]
//	double x_A_B_dt[50][50];												//air boundary absolute humidity
	
//	double Re_A_B_dt[50][50];												//
//	double Pr_A_B_dt[50][50];												//

	//===Transfer properties at previous time step===//
	double q_P_dt[50][50];													//local heat transfer rate [kW]
	double K_P_dt[50][50];													//local overall heat transfer coefficeint [kW/m^2-K]
	
	
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
		
	//===Refrigerant Local Properties current time step===//
	double T_R_P[50][50];													//refrigerant local temperature [degC]
	double P_R_P[50][50];													//refrigerant local pressure [kPa]
	double h_R_P[50][50];													//refrigerant local specific enthalpy [kJ/kg]
	double u_R_P[50][50];													//refrigerant local specific internal energy [kJ/kg]
	double M_R_P[50][50];
	double rho_R_P[50][50];													//refrigerant local density [kg/m^3]
	double myu_R_P[50][50];													//refrigerant local viscosity [Pa-s]
	double k_R_P[50][50];													//refrigerant local thermal conductivity [W/m-k]
	double cp_R_P[50][50];													//refrigerant local specific heat [J/kg-K]
	double x_R_P[50][50];													//refrigerant local quality
	double VF_R_P[50][50];													//refrigerant local void fraction
	double myu_Rl_P[50][50];
	double k_Rl_P[50][50];
	double rho_Rl_P[50][50];
	double cp_Rl_P[50][50];
	double h_Rl_P[50][50];
	double myu_Rv_P[50][50];
	double k_Rv_P[50][50];
	double rho_Rv_P[50][50];
	double cp_Rv_P[50][50];
	double h_Rv_P[50][50];
	double st_Rlv_P[50][50];

	double Re_R_P[50][50];													//
	double Pr_R_P[50][50];													//

	//===Refrigerant Boundary Properties current time step===//
	double G_R_B[50][50];													//refrigerant boundary mass flowrate [kg/s]
	double P_R_B[50][50];													//refrigerant boundary pressure [kPa]
	double h_R_B[50][50];													//refrigerant boundary specific enthalpy [kJ/kg]
	double T_R_B[50][50];													//refrigerant boundary temperature [degC]

	double rho_R_B[50][50];													//refrigerant boundary density [kg/m^3]
	double myu_R_B[50][50];													//refrigerant boundary viscosity [Pa-s]
	double k_R_B[50][50];													//refrigerant boundary thermal conductivity [W/m-k]
	double cp_R_B[50][50];													//refrigerant boundary specific heat [J/kg-K]
	double x_R_B[50][50];													//refrigerant boundary quality

	double Re_R_B[50][50];													//
	double Pr_R_B[50][50];													//

	//===Air Local Properties current time step===//
	double P_A_P[50][50];													//air local ipressure [kPa]
	double h_A_P[50][50];													//air local specific enthalpy [kJ/kg]
	double T_A_P[50][50];													//air local temperature [degC]
	double u_A_P[50][50];													//air local specific internal energy [kJ/kg]
	double v_A_P[50][50];													//air local velocity [m/s]
	double rho_A_P[50][50];													//air local density [kg/m^3]
	double myu_A_P[50][50];													//air local viscosity [Pa-s]
	double k_A_P[50][50];													//air local thermal conductivity [W/m-k]
	double cp_A_P[50][50];													//air local specific heat [J/kg-K]
	double x_A_P[50][50];													//air local absolute humidity

	double Re_A_P[50][50];													//
	double Pr_A_P[50][50];													//
	double Remax_A_P[50][50];

	//===Air Boundary Properties current time step===//
	double G_A_B[50][50];													//air boundary mass flowrate [kg/s]
	double h_A_B[50][50];													//air boundary specific enthalpy [kJ/kg]
	double P_A_B[50][50];													//air boundary pressure [kPa]
	double T_A_B[50][50];													//air boundary temperature [degC]
	double v_A_B[50][50];													//air boundary velocity [m/s]
	double u_A_B[50][50];													//air boundary specific internal energy [kJ/kg]
	double rho_A_B[50][50];													//air boundary density [kg/m^3]
	double myu_A_B[50][50];													//air boundary viscosity [Pa-s]
	double k_A_B[50][50];													//air boundary thermal conductivity [W/m-k]
	double cp_A_B[50][50];													//air boundary specific heat [J/kg-K]
	double x_A_B[50][50];													//air boundary absolute humidity

	double Re_A_B[50][50];													//
	double Pr_A_B[50][50];													//
	
	//===Transfer properties at current time step===//
	double q_P[50][50];														//local heat transfer rate [kW]
	double K_P[50][50];														//local overall heat transfer coefficeint [kW/m^2-K]
	
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
	int i,j,k,p;															//Counter for Newton-Rahpson variables
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

	//===Input Conditions===//
	G_RM_in = 0.015;														//refrigerant mass flowrate
	P_RM_in = 1500.0;														//refrigerant inlet pressure
	h_RM_in = 250.0;														//refrigerant inlet specific enthalpy
	G_AM_in = 0.27;															//refrigerant mass flowrate
	P_AM_in = 101.325;														//refrigerant inlet pressure
	T_AM_in = 35.0;															//refrigerant inlet specific enthalpy

	clock_t start_time, end_time;
	start_time = clock();

	Vmesh=PI*tube_D_o*tube_D_o/4.0*tube_L/double(Nmesh);					//refrigerant mesh volume [m^3]

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
	
	G_R_B[0][1] = G_RM_in;
	h_R_B[0][1] = h_RM_in;
	P_R_B[0][1] = P_RM_in;
	
	for (int j = 1; j <= Ncol; j++)
	{
	for (int i = 1; i <= Nmesh; i++)
	{
		G_R_B[i][j] = G_R_B_dt[i][j];
		h_R_B[i][j] = h_R_B_dt[i][j];
		P_R_B[i][j] = P_R_B_dt[i][j];
		h_R_P[i][j] = h_R_P_dt[i][j];
		P_R_P[i][j] = P_R_P_dt[i][j];
		T_R_P[i][j] = T_R_P_dt[i][j];	
		u_R_P[i][j] = u_R_P_dt[i][j];	
	}
	}
	

	for (int i = 0; i <= Nmesh; i++)
	{
		G_A_B[i][0] = G_AM_in / double(Nmesh) / double(Nrow);
		T_A_B[i][0] = T_AM_in;
		P_A_B[i][0] = P_AM_in;

		for (int j = 1; j <= Ncol; j++)
	{
		//Air-side local and boundary parameters
		T_A_B[i][j] = T_A_B_dt[i][j];
		P_A_B[i][j] = P_A_B_dt[i][j];
		T_A_P[i][j] = T_A_P_dt[i][j];
		P_A_P[i][j] = P_A_P_dt[i][j];
		G_A_B[i][j] = G_A_B_dt[i][j];
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
		mnm.setValue(p, G_R_B[i][j]);
		p++;
		mnm.setValue(p, h_R_P[i][j]);
		p++;
		mnm.setValue(p, P_R_P[i][j]);
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
				G_R_B[i][j] = mnm.getValue(p);
				p++;
				h_R_P[i][j] = mnm.getValue(p);
				p++;
				P_R_P[i][j] = mnm.getValue(p);
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
		if (G_R_B[i][j] >= 0.0)					
	{
		
				ref.change_pure(1);
				ref.state_ph(P_R_P[i][j], h_R_P[i][j]);
				x_R_P[i][j] = ref.Rc.x;
				T_R_P[i][j] = ref.Rc.T;
				u_R_P[i][j] = ref.Rc.u;

		//===Single phase properties===//
		if (x_R_P[i][j] < 0.0 || x_R_P[i][j] > 1.0)
				{

				myu_R_P[i][j] = ref.Rc.visc / 1000000.0;
				k_R_P[i][j] = ref.Rc.thc;
				cp_R_P[i][j] = ref.Rc.cp * 1000.0;
				rho_R_P[i][j] = ref.Rc.rho;
		
		}
				
		//===Two phase properties===//
		else{
		
				ref.change_pure(1);
				ref.sat_p(P_R_P[i][j]);
				myu_Rl_P[i][j] = ref.Rl.visc / 1000000.0;					//in Pa-s
				myu_Rv_P[i][j] = ref.Rv.visc / 1000000.0;
				cp_Rl_P[i][j] = ref.Rl.cp * 1000.0;							//in J/kg-K
				cp_Rv_P[i][j] = ref.Rv.cp * 1000.0;
				k_Rl_P[i][j] = ref.Rl.thc;									//in W/m-k
				k_Rv_P[i][j] = ref.Rv.thc;
				rho_Rl_P[i][j] = ref.Rl.rho;								//in kg/m^3
				rho_Rv_P[i][j] = ref.Rv.rho;
				h_Rl_P[i][j] = ref.Rl.h;									//in kJ/kg
				h_Rv_P[i][j] = ref.Rv.h;
				st_Rlv_P[i][j] = ref.Rl.st;									//in N/m
		
				myu_R_P[i][j] = (x_R_P[i][j] * myu_Rv_P[i][j]) + ((1.0 - x_R_P[i][j]) * myu_Rl_P[i][j]);
				VF_R_P[i][j] = 1.0 / (1.0 + (((1.0 - x_R_P[i][j]) / x_R_P[i][j]) * (rho_Rv_P[i][j] / rho_Rl_P[i][j])));
				rho_R_P[i][j] = (VF_R_P[i][j] * rho_Rv_P[i][j]) + ((1.0 - VF_R_P[i][j]) * rho_Rl_P[i][j]);
	
		}

				Re_R_P[i][j] = G_R_B[i][j] * 4.0 / PI / tube_D_i / myu_R_P[i][j];
				Pr_R_P[i][j] = myu_R_P[i][j] * cp_R_P[i][j] / k_R_P[i][j];

			///===Air-side properties===/// (...A_P[i]=...A_B[i][j] inlet air state for calculation of local air properties)
				ref.change_pure(2);
				ref.state_tp(T_A_B[i][j-1], P_A_B[i][j-1]);

				rho_A_P[i][j] = ref.Rc.rho;
				myu_A_P[i][j] = ref.Rc.visc / 1000000.0;
				k_A_P[i][j] = ref.Rc.thc;
				cp_A_P[i][j] = ref.Rc.cp * 1000.0;
				T_A_P[i][j] = T_A_B[i][j-1];
				P_A_P[i][j] = P_A_B[i][j-1];

			///===Air-side continuity===///
			G_A_B[i][j] = G_A_B[i][j-1];

				v_A_B[i][j] = G_A_B[i][j] / ((tube_L / double(Nmesh) - (int)(tube_L / fin_S) * fin_T) * tube_adjV * rho_A_P[i][j]);		//air inlet velocity
				Re_A_P[i][j] = rho_A_P[i][j] * v_A_B[i][j-1] * tube_D_o / myu_A_P[i][j];
				Remax_A_P[i][j] = rho_A_P[i][j] * Vmax_a * tube_D_o / myu_A_P[i][j];
				Pr_A_P[i][j] = myu_A_P[i][j] * cp_A_P[i][j] / k_A_P[i][j];

			///===simplified overall heat transfer coefficient===///
			K_P[i][j] = 2.0 * PI * tube_D_o * tube_L / double(Nmesh);
		
				///===heat transfer parameters air-side===///

				fin_ft = (4.0 / PI) * (0.25 + (0.118 * pow(Remax_A_P[i][j], -0.16) / (pow((tube_adjV / tube_D_o) - 1.0, 1.08)))) * ((tube_adjV / tube_D_o) - 1);
				fin_ff = 1.455 * pow(Remax_A_P[i][j], -0.656) * pow((tube_adjV / tube_adjH), -0.347) * pow((fin_S / tube_D_o), -0.134) * pow((tube_adjV / tube_D_o), 1.23);
				fin_f = (fin_ff * SA_finned * SA_tot) + (fin_ft * (1.0 - (SA_finned / SA_tot)) * (1.0 - (fin_T / SA_finned)));
				
				
			///===Pressure drop===///
			P_A_B[i][j] = P_A_B[i][j-1] - (fin_f * SA_tot * pow((rho_A_P[i][j] * Vmax_a), 2) / (Ac_fin * 2 * rho_A_P[i][j])) / 1000;

			///===Air outlet temperature===///
			T_A_B[i][j] = ((K_P[i][j] - G_A_B[i][j] * cp_A_P[i][j]) * T_A_B[i][j-1] - K_P[i][j] * T_R_P[i][j])/(K_P[i][j] + G_A_B[i][j] * cp_A_P[i][j]);			
		
			///===Overall heat transfer===///
			q_P[i][j] = K_P[i][j] * ((T_A_B[i][j-1]+T_A_B[i][j])/2.0-T_R_P[i][j]);

		
		///===Continuity error===///
		M_R_P[i][j] = Vmesh * rho_R_P[i][j];

 		mnm.setError(p, (M_R_P[i][j] - M_R_P_dt[i][j]) / dt + G_R_B[i-1][j] - G_R_B[i][j], 0.0);
		p++;

		///===Energy conservation error===/// (h_R_B[i] = h_R_P[i])
		mnm.setError(p, (M_R_P[i][j] * u_R_P[i][j] - M_R_P_dt[i][j] * u_R_P_dt[i][j]) / dt + G_R_B[i-1][j] * h_R_P[i-1][j] - G_R_B[i][j] * h_R_P[i][j] + q_P[i][j], 0.0);
		p++;
		
		if (i == 1){

		///===Pressure drop error===/// Westmost boundary (P_R_B[i] = P_R_P[i])
		mnm.setError(p, P_R_B[i-1][j] - P_R_P[i][j] -128.0 * tube_L / double(Nmesh) /2.0 * G_R_B[i][j] * myu_R_P[i][j] / (PI * rho_R_P[i][j] * pow(tube_D_i, 4.0)), 0.0);
		p++;

		}

		else{

			///Eastmost boundary///
			if(i == Nmesh){

			mnm.setError(p, P_R_P[i][j] - P_R_B[i][j] -128.0 * tube_L / double(Nmesh) /2.0 * G_R_B[i][j] * myu_R_P[i][j] / (PI * rho_R_P[i][j] * pow(tube_D_i, 4.0)), 0.0);
			p++;

			}

			else{

			mnm.setError(p, P_R_P[i-1][j] - P_R_P[i][j] -128.0 * tube_L / double(Nmesh) * G_R_B[i][j] * myu_R_P[i][j] / (PI * rho_R_P[i][j] * pow(tube_D_i, 4.0)), 0.0);
			p++;

			}
		}
		
		
		}

		else{

			ref.change_pure(1);
		ref.state_ph(P_R_P[i][j], h_R_P[i][j]);
		x_R_P[i][j] = ref.Rc.x;
		T_R_P[i][j] = ref.Rc.T;
		u_R_P[i][j] = ref.Rc.u;

		if (x_R_P[i][j] < 0.0 || x_R_P[i][j] > 1.0)
				{
		myu_R_P[i][j] = ref.Rc.visc / 1000000.0;
		k_R_P[i][j] = ref.Rc.thc;
		cp_R_P[i][j] = ref.Rc.cp * 1000.0;
		rho_R_P[i][j] = ref.Rc.rho;
		
		}

		else{
		ref.change_pure(1);
				ref.sat_p(P_R_P[i][j]);
				myu_Rl_P[i][j] = ref.Rl.visc / 1000000.0;					//in Pa-s
				myu_Rv_P[i][j] = ref.Rv.visc / 1000000.0;
				cp_Rl_P[i][j] = ref.Rl.cp * 1000.0;							//in J/kg-K
				cp_Rv_P[i][j] = ref.Rv.cp * 1000.0;
				k_Rl_P[i][j] = ref.Rl.thc;									//in W/m-k
				k_Rv_P[i][j] = ref.Rv.thc;
				rho_Rl_P[i][j] = ref.Rl.rho;								//in kg/m^3
				rho_Rv_P[i][j] = ref.Rv.rho;
				h_Rl_P[i][j] = ref.Rl.h;									//in kJ/kg
				h_Rv_P[i][j] = ref.Rv.h;
				st_Rlv_P[i][j] = ref.Rl.st;								//in N/m
		myu_R_P[i][j] = (x_R_P[i][j] * myu_Rv_P[i][j]) + ((1.0 - x_R_P[i][j]) * myu_Rl_P[i][j]);
		VF_R_P[i][j] = 1.0 / (1.0 + (((1.0 - x_R_P[i][j]) / x_R_P[i][j]) * (rho_Rv_P[i][j] / rho_Rl_P[i][j])));
		rho_R_P[i][j] = (VF_R_P[i][j] * rho_Rv_P[i][j]) + ((1.0 - VF_R_P[i][j]) * rho_Rl_P[i][j]);
		}

		Re_R_P[i][j] = G_R_B[i][j] * 4.0 / PI / tube_D_i / myu_R_P[i][j];
					Pr_R_P[i][j] = myu_R_P[i][j] * cp_R_P[i][j] / k_R_P[i][j];

		///===Air-side properties===/// (...A_P[i]=...A_B[i][j] inlet air state for calculation of local air properties)
				ref.change_pure(2);
				ref.state_tp(T_A_B[i][j-1], P_A_B[i][j-1]);

				rho_A_P[i][j] = ref.Rc.rho;
				myu_A_P[i][j] = ref.Rc.visc / 1000000.0;
				k_A_P[i][j] = ref.Rc.thc;
				cp_A_P[i][j] = ref.Rc.cp * 1000.0;
				T_A_P[i][j] = T_A_B[i][j-1];
				P_A_P[i][j] = P_A_B[i][j-1];

				G_A_B[i][j] = G_A_B[i][j-1];

				v_A_B[i][j] = G_A_B[i][j] / ((tube_L / double(Nmesh) - (int)(tube_L / fin_S) * fin_T) * tube_adjV * rho_A_P[i][j]);		//air inlet velocity
				Re_A_P[i][j] = rho_A_P[i][j] * v_A_B[i][j-1] * tube_D_o / myu_A_P[i][j];
				Remax_A_P[i][j] = rho_A_P[i][j] * Vmax_a * tube_D_o / myu_A_P[i][j];
				Pr_A_P[i][j] = myu_A_P[i][j] * cp_A_P[i][j] / k_A_P[i][j];

			///===simplified overall heat transfer coefficient===///
			K_P[i][j] = 2.0 * PI * tube_D_o * tube_L / double(Nmesh);
		
				///===heat transfer air-side===///

				fin_ft = (4.0 / PI) * (0.25 + (0.118 * pow(Remax_A_P[i][j], -0.16) / (pow((tube_adjV / tube_D_o) - 1.0, 1.08)))) * ((tube_adjV / tube_D_o) - 1);
				fin_ff = 1.455 * pow(Remax_A_P[i][j], -0.656) * pow((tube_adjV / tube_adjH), -0.347) * pow((fin_S / tube_D_o), -0.134) * pow((tube_adjV / tube_D_o), 1.23);
				fin_f = (fin_ff * SA_finned * SA_tot) + (fin_ft * (1.0 - (SA_finned / SA_tot)) * (1.0 - (fin_T / SA_finned)));
				
				

			P_A_B[i][j] = P_A_B[i][j-1] - (fin_f * SA_tot * pow((rho_A_P[i][j] * Vmax_a), 2.0) / (Ac_fin * 2 * rho_A_P[i][j])) / 1000.0;

			T_A_B[i][j] = ((K_P[i][j] - G_A_B[i][j] * cp_A_P[i][j]) * T_A_B[i][j-1] - K_P[i][j] * T_R_P[i][j])/(K_P[i][j] + G_A_B[i][j] * cp_A_P[i][j]);			
		
			q_P[i][j] = K_P[i][j] * ((T_A_B[i][j-1]+T_A_B[i][j])/2.0-T_R_P[i][j]);

		///===Continuity error===///
		M_R_P[i][j] = Vmesh * rho_R_P[i][j];

 		mnm.setError(p, (M_R_P[i][j] - M_R_P_dt[i][j]) / dt + G_R_B[i][j] - G_R_B[i+1][j], 0.0);
		p++;

		///===Energy conservation error===/// (h_R_B[i] = h_R_P[i])
		mnm.setError(p, (M_R_P[i][j] * u_R_P[i][j] - M_R_P_dt[i][j] * u_R_P_dt[i][j]) / dt + G_R_B[i][j] * h_R_P[i][j] - G_R_B[i+1][j] * h_R_P[i+1][j] + q_P[i][j], 0.0);
		p++;
		
		if (i == 1){

		///===Pressure drop error===/// Westmost boundary (P_R_B[i] = P_R_P[i])
		mnm.setError(p, P_R_B[i-1][j] - P_R_P[i][j] + 128.0 * tube_L / double(Nmesh) /2.0 * G_R_B[i][j] * myu_R_P[i][j] / (PI * rho_R_P[i][j] * pow(tube_D_i, 4.0)), 0.0);
		p++;
		
		}

		else{

			///Eastmost boundary///
			if(i == Nmesh){

			mnm.setError(p, P_R_P[i][j] - P_R_B[i][j] + 128.0 * tube_L / double(Nmesh) /2.0 * G_R_B[i][j] * myu_R_P[i][j] / (PI * rho_R_P[i][j] * pow(tube_D_i, 4.0)), 0.0);
			p++;

			}

			else{

			mnm.setError(p, P_R_P[i-1][j] - P_R_P[i][j] + 128.0 * tube_L / double(Nmesh) * G_R_B[i][j] * myu_R_P[i][j] / (PI * rho_R_P[i][j] * pow(tube_D_i, 4.0)), 0.0);
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
		if (G_R_B[i][j] >= 0.0)
	{

		h_R_B[i][j] = h_R_P[i][j];
		P_R_B[i][j] = P_R_P[i][j];

		}

		else{
			
		h_R_B[i][j] = h_R_P[i+1][j];
		P_R_B[i][j] = P_R_P[i+1][j];

		}
		T_A_B[i][j] = T_A_P[i][j];
		P_A_B[i][j] = P_A_P[i][j];

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
		Qtot += q_P[i][j];

		cout << "\n\nAt tube_" << i << "\trefrigerant's\tP_ = " << P_R_P[i][j] <<  "\t\th_ = " << h_R_P[i][j] <<  "\t\tx_ = " << x_R_P[i][j] << "\t\tT_ = " << T_R_P[i][j] << "\t\tG_ = " << G_R_B[i][j];
		cout << "\n\t\tair's\t\tP_ = " << P_A_P[i][j] << "\t\tT_ = " << T_A_B[i][j];
	}
	cout << endl;

	Atot = (tube_D_o * PI * (tube_L - (int)(tube_L / fin_S) * fin_T) + (int)(tube_L / fin_S) * 2.0 * (tube_adjV * tube_adjH - PI / 4.0 * tube_D_o * tube_D_o))* double(Ntubes);
	cout << "\nQtot = " << Qtot << " kW\tDP = " << P_R_B[0][1]  - P_R_B[Nmesh][1]<< " kPa ";

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
		G_R_B_dt[i][j] = G_R_B[i][j];
		h_R_B_dt[i][j] = h_R_B[i][j];
		P_R_B_dt[i][j] = P_R_B[i][j];
		h_R_P_dt[i][j] = h_R_P[i][j];
		P_R_P_dt[i][j] = P_R_P[i][j];
		T_R_P_dt[i][j] = T_R_P[i][j];	
		u_R_P_dt[i][j] = u_R_P[i][j];
		M_R_P_dt[i][j] = M_R_P[i][j];
		q_P_dt[i][j] = q_P[i][j];
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
		T_A_B_dt[i][j] = T_A_B[i][j];
		P_A_B_dt[i][j] = P_A_B[i][j];
		T_A_P_dt[i][j] = T_A_P[i][j];
		P_A_P_dt[i][j] = P_A_P[i][j];
		G_A_B_dt[i][j] = G_A_B[i][j];
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

*/
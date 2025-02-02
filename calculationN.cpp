/*#include <math.h>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

#include "CNewtonRaphsonMethod.h"
#include "DXProperty_ver06.h"
#include "matrix.h"
#include "parameter.h"

double run_time[1];

double solution(int ref_fluid, int type, int ref_cor, int Ntubes, int Nrow, int Ncol)
{
	//Variables for Refrigerant
	double T_r_inlet, x_r_inlet, myu_r_inlet, rho_r_inlet;
	double myu_rl, myu_rv, rho_rl, rho_rv, cp_rl, cp_rv, k_rl, k_rv, h_rl, h_rv, st_rl, st_rv, VF;
	double myu_r, k_r, cp_r, rho_r, Re_r, Pr_r;
	double DP;
	double Ffl, MW, Pcrit;

	//Variables for Air
	double SA_finned, Ac_fin, SA_tot;
	double rho_a_inlet, myu_a_inlet, k_a_inlet, cp_a_inlet, h_a_inlet, V_a_inlet;
	double Sd, Vmax_a;
	double Re_a, Re_a_max, Pr_a;
	double fin_ft, fin_ff, fin_f, fin_DP, jcf, jcf_t;	
	double rho_a, myu_a, k_a, cp_a;

	//Other variables
	const double PI = 6 * (asin(0.5));
	const double eta_fin = 0.77, Rc = 0.0, Rf = 0.0;
	double htc, alfa_lo, mass_flux, Co, Fr, Xtt, Pr;
	int p;										//counter for Newton-Raphson variables
	double G_res, H_res;
	
	DXProperty refr;
	DXProperty refa;

	if (ref_fluid == 1) { refr.reibai = "R11.FLD"; Ffl = 1.30, MW = 137.4; Pcrit = 4394.0; }
	else if (ref_fluid == 2) { refr.reibai = "R12.FLD"; Ffl = 1.50; MW = 120.91; Pcrit = 4136.1; }
	else if (ref_fluid == 3) { refr.reibai = "R22.FLD"; Ffl = 2.20; MW = 86.47; Pcrit = 4990.0; }
	else if (ref_fluid == 4) { refr.reibai = "R113.FLD"; Ffl = 1.30; MW = 187.37; Pcrit = 3392.2; }
	else if (ref_fluid == 5) { refr.reibai = "R114.FLD"; Ffl = 1.24; MW = 170.92; Pcrit = 3257.0; }
	else if (ref_fluid == 6) { refr.reibai = "R134A.FLD"; Ffl = 1.63; MW = 102.03; Pcrit = 4059.28; }
	else if (ref_fluid == 7) { refr.reibai = "WATER.FLD"; Ffl = 1.10; MW = 18.02; Pcrit = 22064.00; }

	refa.reibai = "AIR.PPF";
	refr.joutai = "DFT";
	refa.joutai = "IIR";
	refr.LoadDLL("refprop.DLL");
	refa.LoadDLL("refprop2.DLL");
	refr.setup();
	refa.setup();

	SA_finned = (int)(tube_L / fin_S) * 2 * tube_adjH * tube_adjV;
	Ac_fin = (tube_L - (((int)(tube_L / fin_S)) * fin_T)) * tube_adjV;
	SA_tot = SA_finned + ((tube_L - (((int)(tube_L / fin_S))* fin_T)) * PI * tube_D_o);
	
	refr.state_ph(P_r_inlet, h_r_inlet);
	T_r_inlet = refr.Rc.T;
	x_r_inlet = refr.Rc.x;

	refa.state_tp(T_a_inlet, P_a_inlet);
	rho_a_inlet = refa.Rc.rho;
	myu_a_inlet = refa.Rc.visc / 1000000;
	k_a_inlet = refa.Rc.thc;
	cp_a_inlet = refa.Rc.cp * 1000;
	h_a_inlet = refa.Rc.h;
	
	V_a_inlet = G_a_total / ((tube_L - (int)(tube_L / fin_S) * fin_T) * tube_adjV * Nrow * rho_a_inlet);
	Sd = pow(((pow(tube_adjH, 0.5)) + (pow((tube_adjV / 2), 0.5))), 0.5);

	if (Sd < ((tube_adjV + tube_D_o) / 2))
	{
		Vmax_a = (tube_adjV / (2 * (Sd - tube_D_o))) * V_a_inlet;
	}

	else
	{
		Vmax_a = (tube_adjV / (tube_adjV - tube_D_o)) * V_a_inlet;
	}
	
	//===Local Variables Initialization===//
	for (int i = 1; i <= Ntubes; i++)
	{
		//Air-side parameters
		T_in_a[BFS_path[i]][0] = T_a_inlet;
		T_out_a[BFS_path[i]][0] = T_a_inlet;
		P_in_a[BFS_path[i]][0] = P_a_inlet;
		G_a[BFS_path[i]][0] = G_a_total / Nrow;

		//Refrigerant-side parameters
		h_in_r[BFS_path[i]][0] = h_r_inlet;
		h_out_r[BFS_path[i]][0] = h_r_inlet;
		T_in_r[BFS_path[i]][0] = T_r_inlet;
		T_out_r[BFS_path[i]][0] = T_r_inlet;

		if (row_neg[BFS_path[i]] == 0)
		{
			G_r[BFS_path[i]][0] = G_r_total / Ninlet[0];
			P_in_r[BFS_path[i]][0] = P_r_inlet;
			x_in_r[BFS_path[i]][0] = x_r_inlet;

			if (x_r_inlet < 0 || x_r_inlet > 1)
			{
				refr.state_ph(P_r_inlet, h_r_inlet);
				myu_r_inlet = refr.Rc.visc / 1000000;
				rho_r_inlet = refr.Rc.rho;
			}

			else
			{
				refr.sat_p(P_r_inlet);
				myu_rl = refr.Rl.visc / 1000000;
				myu_rv = refr.Rv.visc / 1000000;
				rho_rl = refr.Rl.rho;
				rho_rv = refr.Rv.rho;
				
				VF = 1 / (1 + (((1 - x_r_inlet) / x_r_inlet) * (rho_rv / rho_rl)));

				myu_r_inlet = (x_r_inlet * myu_rv) + ((1 - x_r_inlet) * myu_rl);				
				rho_r_inlet = (VF * rho_rv) + ((1 - VF) * rho_rl);
			}

			DP = 128 * tube_L * G_r[BFS_path[i]][0] * myu_r_inlet / (PI * rho_r_inlet * pow(tube_D_i, 4.0) * 1000);		//in kPa
			P_out_r[BFS_path[i]][0] = P_in_r[BFS_path[i]][0] - DP;

			refr.state_ph(P_out_r[BFS_path[i]][0], h_r_inlet);
			x_out_r[BFS_path[i]][0] = refr.Rc.x;
		}

		else
		{
			for (int j = 1; j <= Ntubes; j++)
			{
				if (TTCM[BFS_path[i]][j] == -1)
				{
					if (row_neg[BFS_path[i]] > 1)
					{
						G_r[BFS_path[i]][0] = G_r[BFS_path[i]][0] + G_r[j][0];

						double temp_P_in = P_out_r[j][0];

						for (int k = 1; k <= Ntubes; k++)		//Search for the lowest outlet pressure among the merging flows
						{
							if (TTCM[BFS_path[i]][k] == -1 && k != j)
							{
								if (P_out_r[k][0] < temp_P_in && P_out_r[k][0] > 0)
								{
									temp_P_in = P_out_r[k][0];
								}
							}
						}

						for (int k = 1; k <= Ntubes; k++)
						{
							if (TTCM[BFS_path[i]][k] == -1) { P_out_r[k][0] = temp_P_in; }
						}

						P_in_r[BFS_path[i]][0] = temp_P_in;
						P_out_r[BFS_path[i]][0] = P_in_r[BFS_path[i]][0] - DP;

						refr.state_ph(P_in_r[BFS_path[i]][0], h_r_inlet);
						x_in_r[BFS_path[i]][0] = refr.Rc.x;

						refr.state_ph(P_out_r[BFS_path[i]][0], h_r_inlet);
						x_out_r[BFS_path[i]][0] = refr.Rc.x;
					}

					else if (row_pos[j] > 1)
					{
						G_r[BFS_path[i]][0] = G_r[j][0] / row_pos[j];
						P_in_r[BFS_path[i]][0] = P_out_r[j][0];
						P_out_r[BFS_path[i]][0] = P_in_r[BFS_path[i]][0] - DP;

						refr.state_ph(P_in_r[BFS_path[i]][0], h_r_inlet);
						x_in_r[BFS_path[i]][0] = refr.Rc.x;

						refr.state_ph(P_out_r[BFS_path[i]][0], h_r_inlet);
						x_out_r[BFS_path[i]][0] = refr.Rc.x;
					}

					else
					{
						G_r[BFS_path[i]][0] = G_r[j][0];
						P_in_r[BFS_path[i]][0] = P_out_r[j][0];
						P_out_r[BFS_path[i]][0] = P_in_r[BFS_path[i]][0] - DP;

						refr.state_ph(P_in_r[BFS_path[i]][0], h_r_inlet);
						x_in_r[BFS_path[i]][0] = refr.Rc.x;

						refr.state_ph(P_out_r[BFS_path[i]][0], h_r_inlet);
						x_out_r[BFS_path[i]][0] = refr.Rc.x;
					}
				}
			}
		}
	
		Re_a = rho_a_inlet * V_a_inlet * tube_D_o / myu_a_inlet;
		Re_a_max = rho_a_inlet * Vmax_a * tube_D_o / myu_a_inlet;
		Pr_a = myu_a_inlet * cp_a_inlet / k_a_inlet;

		//Air-side pressure drop
		fin_ft = (4 / PI) * (0.25 + (0.118 * pow(Re_a_max, -0.16) / (pow((tube_adjV / tube_D_o) - 1, 1.08)))) * ((tube_adjV / tube_D_o) - 1);
		fin_ff = 1.455 * pow(Re_a_max, -0.656) * pow((tube_adjV / tube_adjH), -0.347) * pow((fin_S / tube_D_o), -0.134) * pow((tube_adjV / tube_D_o), 1.23);
		fin_f = (fin_ff * SA_finned * SA_tot) + (fin_ft * (1 - (SA_finned / SA_tot)) * (1 - (fin_T / SA_finned)));
		fin_DP = (fin_f * SA_tot * pow((rho_a_inlet * Vmax_a), 2) / (Ac_fin * 2 * rho_a_inlet)) / 1000;
		P_out_a[BFS_path[i]][0] = P_in_a[BFS_path[i]][0] - fin_DP;

		//Air-side heat transfer (Kim-Youn-Webb)
		jcf = 0.163 * pow(Re_a_max, -0.369) * pow((tube_adjV / tube_adjH), 0.106) * pow((fin_S / tube_D_o), 0.0138) * pow((tube_adjV / tube_D_o), 0.13);

		if (Nrow >= 3) { jcf_t = jcf; }

		else
		{
			jcf_t = 1.043 * pow(0.163 * pow(Re_a_max, -0.14) * pow((tube_adjV / tube_adjH), -0.564) * pow((fin_S / tube_D_o), -0.123) * pow((tube_adjV / tube_D_o), 1.17), 3 - Nrow) * jcf;
		}

		alfa_a[BFS_path[i]][0] = jcf_t * (rho_a_inlet * Vmax_a * cp_a_inlet) * pow(Pr_a, 1.5) / 1000;

		//Refrigerant-side heat transfer
		htc = (1.0 / (tube_k * PI * tube_L)) + (Rc / (PI * tube_D_o * tube_L)) + (Rf / (PI * tube_D_o * tube_L)) + (1.0 / (alfa_a[BFS_path[i]][0] * (tube_D_o * PI * (tube_L - (int)(tube_L / fin_S) * fin_T) + (int)(tube_L / fin_S) * 2.0 * (tube_adjV * tube_adjH - PI / 4.0 * tube_D_o * tube_D_o)) * eta_fin));

		refr.sat_p(P_r_inlet);
		myu_rl = refr.Rl.visc / 1000000;				//in Pa-s
		myu_rv = refr.Rv.visc / 1000000;
		cp_rl = refr.Rl.cp * 1000;						//in J/kg-K
		k_rl = refr.Rl.thc;								//in W/m-K
		rho_rl = refr.Rl.rho;							//in kg/m^3
		rho_rv = refr.Rv.rho;

		if (type == 1)
		{
			if (ref_cor == 1)		//Kandlikar
			{				
				mass_flux = 4 * G_r[BFS_path[i]][0] / (PI * pow(tube_D_i, 2));
				alfa_lo = 0.023 * k_rl * (pow(((mass_flux * tube_D_i * (1 - x_in_r[BFS_path[i]][0])) / myu_rl), 0.8)) * (pow((cp_rl * myu_rl / k_rl), 0.4)) / (tube_D_i * 1000);			//in kJ/kg

				Co = (pow(((1 - x_in_r[BFS_path[i]][0]) / x_in_r[BFS_path[i]][0]), 0.8)) * (pow((rho_rv / rho_rl), 0.5));
				Fr = pow(mass_flux, 2) / (9.81 * tube_D_i * (pow(rho_rl, 2)));
				double C;

				if (Fr > 0.04) { C = 0; }
				else { C = 0.3; }

				double conv = 1.1360 * pow(Co, -0.9) * pow((25 * Fr), C) * alfa_lo;
				double boil = 0.0683 * pow(Co, -0.2) * pow((25 * Fr), C) * alfa_lo;

				if (conv > boil) { alfa_r[BFS_path[i]][0] = conv; }
				else { alfa_r[BFS_path[i]][0] = boil; }

				Q[BFS_path[i]][0] = pow((htc + pow((alfa_r[BFS_path[i]][0] * PI * tube_D_i * tube_L), -1)), -1) * (T_in_a[BFS_path[i]][0] - T_in_r[BFS_path[i]][0]);
				Q_conv[BFS_path[i]][0] = Q[BFS_path[i]][0];
				Q_boil[BFS_path[i]][0] = Q[BFS_path[i]][0];
			}

			else if (ref_cor == 2)	//Gungor-Winterton
			{
				Xtt = (pow(((1 - x_in_r[BFS_path[i]][0]) / x_in_r[BFS_path[i]][0]), 0.9)) * (pow((rho_rv / rho_rl), 0.5)) * (pow((myu_rl / myu_rv), 0.1));
				mass_flux = 4 * G_r[BFS_path[i]][0] / (PI * pow(tube_D_i, 2));
				Fr = pow(mass_flux, 2) / (9.81 * tube_D_i * pow(rho_rl, 2));																																						// in g/mol
				alfa_lo = 0.023 * k_rl * (pow(((mass_flux * tube_D_i * (1 - x_in_r[BFS_path[i]][0])) / myu_rl), 0.8)) * (pow((cp_rl * myu_rl / k_rl), 0.4)) / (tube_D_i * 1000);			//in kJ/kg

				if (Fr < 0.05) { alfa_r[BFS_path[i]][0] = (alfa_lo * (1 + (1.37 * (pow((1 / Xtt), 0.86))))) * (pow(Fr, (0.1 - (2 * Fr)))) / 1000; }
				else { alfa_r[BFS_path[i]][0] = alfa_lo * (1 + (1.37 * (pow((1 / Xtt), 0.86)))) / 1000; }
				
				Q[BFS_path[i]][0] = pow((htc + pow((alfa_r[BFS_path[i]][0] * PI * tube_D_i * tube_L), -1)), -1) * (T_in_a[BFS_path[i]][0] - T_in_r[BFS_path[i]][0]);
				Q_conv[BFS_path[i]][0] = Q[BFS_path[i]][0];
				Q_boil[BFS_path[i]][0] = Q[BFS_path[i]][0];			
			}

			else if (ref_cor == 3)	//Liu-Winterton
			{
				Pr = myu_rl * cp_rl / k_rl;
				Xtt = (pow(((1 - x_in_r[BFS_path[i]][0]) / x_in_r[BFS_path[i]][0]), 0.9)) * (pow((rho_rv / rho_rl), 0.5)) * (pow((myu_rl / myu_rv), 0.1));
				mass_flux = 4 * G_r[BFS_path[i]][0] / (PI * pow(tube_D_i, 2));
				Fr = pow(mass_flux, 2) / (9.81 * tube_D_i * pow(rho_rl, 2));
				alfa_lo = 0.023 * k_rl * (pow(((mass_flux * tube_D_i * (1 - x_in_r[BFS_path[i]][0])) / myu_rl), 0.8)) * (pow((cp_rl * myu_rl / k_rl), 0.4)) / (tube_D_i * 1000);			//in kJ/kg
				double F = pow((1 + (x_in_r[BFS_path[i]][0] * Pr * ((rho_rl / rho_rv) - 1))), 0.35);
				
				if (Fr < 0.05) { alfa_r[BFS_path[i]][0] = (alfa_lo * (pow(Fr, (0.1 - (2 * Fr)))) * (pow(F, 0.35))) / 1000; }
				else { alfa_r[BFS_path[i]][0] = (alfa_lo * (pow(F, 0.35))) / 1000; }

				Q[BFS_path[i]][0] = pow((htc + pow((alfa_r[BFS_path[i]][0] * PI * tube_D_i * tube_L), -1)), -1) * (T_in_a[BFS_path[i]][0] - T_in_r[BFS_path[i]][0]);
				Q_conv[BFS_path[i]][0] = Q[BFS_path[i]][0];
				Q_boil[BFS_path[i]][0] = Q[BFS_path[i]][0];
			}

			else if (ref_cor == 4)	//Shah
			{
				mass_flux = 4 * G_r[BFS_path[i]][0] / (PI * pow(tube_D_i, 2));
				alfa_lo = 0.023 * k_rl * (pow(((mass_flux * tube_D_i * (1 - x_in_r[BFS_path[i]][0])) / myu_rl), 0.8)) * (pow((cp_rl * myu_rl / k_rl), 0.4)) / (tube_D_i * 1000);			//in kJ/kg

				Co = (pow(((1 - x_in_r[BFS_path[i]][0]) / x_in_r[BFS_path[i]][0]), 0.8)) * (pow((rho_rv / rho_rl), 0.5));
				Fr = pow(mass_flux, 2) / (9.81 * tube_D_i * (pow(rho_rl, 2)));
				double N;

				if (Fr > 0.04) { N = Co; }
				else { N = 0.38 * (pow(Fr, -0.3)) * Co; }

				alfa_r[BFS_path[i]][0] = alfa_lo * 1.8 * pow(N, -0.8);

				Q[BFS_path[i]][0] = pow((htc + pow((alfa_r[BFS_path[i]][0] * PI * tube_D_i * tube_L), -1)), -1) * (T_in_a[BFS_path[i]][0] - T_in_r[BFS_path[i]][0]);
				Q_conv[BFS_path[i]][0] = Q[BFS_path[i]][0];
				Q_boil[BFS_path[i]][0] = Q[BFS_path[i]][0];
			}
		}

		else
		{
			//in progress
			refr.state_ph(P_in_r[BFS_path[i]][0], h_in_r[BFS_path[i]][0]);
			myu_r = refr.Rc.visc / 1000000;
			k_r = refr.Rc.thc;
			cp_r = refr.Rc.cp * 1000;
			rho_r = refr.Rc.rho;

			Re_r = (G_r[BFS_path[i]][0] * 4.0) / (PI * tube_D_i * myu_r);
			Pr_r = myu_r * cp_r / k_r;
			alfa_r[BFS_path[i]][0] = 0.023 * k_r * pow(Re_r, 0.8) * pow(Pr_r, 0.3) / (tube_D_i * 1000);

			Q[BFS_path[i]][0] = (T_in_a[BFS_path[i]][0] - T_in_r[BFS_path[i]][0]) * (pow(((pow((alfa_r[BFS_path[i]][0] * PI * tube_D_i * tube_L), -1)) + htc), -1));
			Q_conv[BFS_path[i]][0] = Q[BFS_path[i]][0];
			Q_boil[BFS_path[i]][0] = Q[BFS_path[i]][0];
		}

	}

	if (Noutlet[0] > 1)				//All outlet tubes have the same initialized outlet pressure
	{
		for (int i = 1; i <= Ntubes; i++)
		{
			if (row_pos[i] == 0)
			{
				double temp_P_out = P_out_r[i][0];

				for (int j = 1; j <= Ntubes; j++)
				{
					if (row_pos[j] == 0 && j != i)
					{
						if (P_out_r[j][0] < temp_P_out) { temp_P_out = P_out_r[j][0]; }
					}
				}

				for (int j = 1; j <= Ntubes; j++)
				{
					if (row_pos[j] == 0) { P_out_r[j][0] = temp_P_out; }
				}
			}
		}
	}

	//==Printout of Initialized Variables==//
	for (int i = 1; i <= Ntubes; i++)
	{
		cout << "\nAt tube " << i << "\trefrigerant's\tG = " << G_r[i][0] << "\tP_in = " << P_in_r[i][0] << "\tP_out = " << P_out_r[i][0] << "\t\tT_in = " << T_in_r[i][0] << "\t\tT_out = " << T_out_r[i][0];
		cout << "\n\t\t\t\th_in = " << h_in_r[i][0] << "\th_out = " << h_out_r[i][0] << "\tx_in = " << x_in_r[i][0] << "\tx_out = " << x_out_r[i][0];
		cout << "\n\t\tair's\t\tG = " << G_a[i][0] << "\tP_in = " << P_in_a[i][0] << "\tP_out = " << P_out_a[i][0] << "\t\tT_in = " << T_in_a[i][0] << "\t\tT_out = " << T_out_a[i][0];
		cout << "\n";
	}

	clock_t start_time, end_time;
	start_time = clock();

	//===Newton-Raphson Calculation===//
	CNewtonRaphsonMethod mnm;
	mnm.setup(2000, 1 + 1e-12, 1e-4);

	p = 0;

	for (int a = 1; a <= Ntubes; a++)
	{
		int i = BFS_path[a];
		mnm.setValue(p, G_r[i][0]);
		p++;
		mnm.setValue(p, P_out_r[i][0]);
		p++;
		mnm.setValue(p, Q_conv[i][0]);
		p++;
		mnm.setValue(p, Q_boil[i][0]);
		p++;
	}

	mnm.setAcc(1);
	mnm.initial();

	for (mnm.main_loop_init(); mnm.main_loop_check(); mnm.main_loop_reinit())
	{
		for (mnm.sub_loop_init(); mnm.sub_loop_check(); mnm.sub_loop_reinit())
		{
			p = 0;

			for (int a = 1; a <= Ntubes; a++)
			{
				int i = BFS_path[a];
				G_r[i][0] = mnm.getValue(p);
				p++;
				P_out_r[i][0] = mnm.getValue(p);
				p++;
				Q_conv[i][0] = mnm.getValue(p);
				p++;
				Q_boil[i][0] = mnm.getValue(p);
				p++;
			}

			p = 0;

			//Fundamental equations, tube loop
			for (int a = 1; a <= Ntubes; a++)
			{
				int i = BFS_path[a];

				//Air parameters
				refa.state_tp(T_in_a[i][0], P_in_a[i][0]);

				rho_a = refa.Rc.rho;
				myu_a = refa.Rc.visc / 1000000;
				k_a = refa.Rc.thc;
				cp_a = refa.Rc.cp * 1000;
				h_in_a[i][0] = refa.Rc.h;

				Re_a = rho_a * V_a_inlet * tube_D_o / myu_a;
				Re_a_max = rho_a * Vmax_a * tube_D_o / myu_a;
				Pr_a = myu_a * cp_a / k_a;

				//Refrigerant parameters
				refr.sat_p(P_in_r[i][0]);
				myu_rl = refr.Rl.visc / 1000000;				//in Pa-s
				myu_rv = refr.Rv.visc / 1000000;
				cp_rl = refr.Rl.cp * 1000;						//in J/kg-K
				cp_rv = refr.Rv.cp * 1000;
				k_rl = refr.Rl.thc;								//in W/m-k
				k_rv = refr.Rv.thc;
				rho_rl = refr.Rl.rho;							//in kg/m^3
				rho_rv = refr.Rv.rho;
				h_rl = refr.Rl.h;								//in kJ/kg
				h_rv = refr.Rv.h;
				st_rl = refr.Rl.st;								//in N/m
				st_rv = refr.Rv.st;

				refr.state_ph(P_in_r[i][0], h_in_r[i][0]);
				x_in_r[i][0] = refr.Rc.x;

				if (x_in_r[i][0] < 0 || x_in_r[i][0] > 1)
				{
					refr.state_ph(P_in_r[i][0], h_in_r[i][0]);
					myu_r = refr.Rc.visc / 1000000;
					k_r = refr.Rc.thc;
					cp_r = refr.Rc.cp * 1000;
					rho_r = refr.Rc.rho;

					Re_r = (G_r[i][0] * 4.0) / (PI * tube_D_i * myu_r);
					Pr_r = myu_r * cp_r / k_r;
				}

				else
				{
					myu_r = (x_in_r[i][0] * myu_rv) + ((1 - x_in_r[i][0]) * myu_rl);
					VF = 1 / (1 + (((1 - x_in_r[i][0]) / x_in_r[i][0]) * (rho_rv / rho_rl)));
					rho_r = (VF * rho_rv) + ((1 - VF) * rho_rl);
					
					Re_r = (G_r[i][0] * 4.0) / (PI * tube_D_i * myu_r);
				}

				//Continuity and momentum equations						//observes the dowstream (destination) tube(s)
				G_res = 0.0;

				if (row_pos[i] == 0)									//outlet tube
				{
					if (Noutlet[0] > 1)									//multiple outlet tubes
					{
						if (i == ref_out[0])							//reference outlet tube
						{
							G_res = G_r_total;

							for (int j = 1; j <= Ntubes; j++)
							{
								if (row_pos[j] == 0 && j != 1)
								{
									G_res = G_res - G_r[j][0];
								}
							}

							mnm.setError(p, G_r[i][0], G_res);
							p++;

							mnm.setError(p, -128.0 * tube_L * G_r[i][0] * myu_r / (PI * rho_r * pow(tube_D_i, 4.0)), (P_out_r[i][0] - P_in_r[i][0]) * 1000.0);
							p++;
						}

						else												//non-reference outlet tube(s)
						{
							mnm.setError(p, G_r[i][0], (P_out_r[i][0] - P_in_r[i][0]) * 1000.0 * PI * rho_r * pow(tube_D_i, 4.0) / (-128.0 * tube_L * myu_r));
							p++;

							mnm.setError(p, P_out_r[i][0], P_out_r[ref_out[0]][0]);		//should have the same outlet pressure with the reference outlet tube
							p++;
						}
					}

					else													//single outlet tube
					{
						mnm.setError(p, G_r[i][0], G_r_total);
						p++;

						mnm.setError(p, -128.0 * tube_L * G_r[i][0] * myu_r / (PI * rho_r * pow(tube_D_i, 4.0)), (P_out_r[i][0] - P_in_r[i][0]) * 1000.0);
						p++;
					}
				}

				else if (merge[i] > 0)										//dowstream tube is a merge
				{
					int end_merge = 0, ref_merge = 0;

					for (int j = 1; j <= Ntubes; j++)
					{
						if (TTCM[i][j] == 1 && row_neg[j] > 1)
						{
							end_merge = j;
							break;
						}
					}

					for (int j = 1; j <= Ntubes; j++)						//looks for the reference tube among the merging flows
					{
						if (TTCM[end_merge][j] == -1 && merge[j] == 2)
						{
							ref_merge = j;
							break;
						}
					}

					if (merge[i] == 2)										//reference tube
					{
						G_res = G_r[end_merge][0];

						for (int j = 1; j <= Ntubes; j++)
						{
							if (TTCM[end_merge][j] == -1 && j != i)
							{
								G_res = G_res - G_r[j][0];
							}
						}

						mnm.setError(p, G_r[i][0], G_res);
						p++;

						mnm.setError(p, -128.0 * tube_L * G_r[i][0] * myu_r / (PI * rho_r * pow(tube_D_i, 4.0)), (P_out_r[i][0] - P_in_r[i][0]) * 1000.0);
						p++;
					}

					else
					{
						mnm.setError(p, G_r[i][0], (P_out_r[i][0] - P_in_r[i][0]) * 1000.0 * PI * rho_r * pow(tube_D_i, 4.0) / (-128.0 * tube_L * myu_r));
						p++;

						mnm.setError(p, P_out_r[i][0], P_out_r[ref_merge][0]);		//Should have the same outlet pressure with the reference tube
						p++;
					}
				}

				else
				{
					for (int j = 1; j <= Ntubes; j++)					//observes the downstream (destination) flow(s)
					{
						if (TTCM[i][j] == 1)
						{
							G_res = G_res + G_r[j][0];
						}
					}

					mnm.setError(p, G_r[i][0], G_res);
					p++;

					mnm.setError(p, -128.0 * tube_L * G_r[i][0] * myu_r / (PI * rho_r * pow(tube_D_i, 4.0)), (P_out_r[i][0] - P_in_r[i][0]) * 1000.0);
					p++;
				}
				
				//Air-side heat transfer (Kim-Youn-Webb)
				jcf = 0.163 * pow(Re_a_max, -0.369) * pow((tube_adjV / tube_adjH), 0.106) * pow((fin_S / tube_D_o), 0.0138) * pow((tube_adjV / tube_D_o), 0.13);

				if (Nrow >= 3) { jcf_t = jcf; }

				else
				{
					jcf_t = 1.043 * pow(0.163 * pow(Re_a_max, -0.14) * pow((tube_adjV / tube_adjH), -0.564) * pow((fin_S / tube_D_o), -0.123) * pow((tube_adjV / tube_D_o), 1.17), 3 - Nrow) * jcf;
				}

				alfa_a[i][0] = jcf_t * (rho_a * Vmax_a * cp_a) * pow(Pr_a, 1.5) / 1000;

				//Air-side pressure drop
				fin_ft = (4 / PI) * (0.25 + (0.118 * pow(Re_a_max, -0.16) / (pow((tube_adjV / tube_D_o) - 1, 1.08)))) * ((tube_adjV / tube_D_o) - 1);
				fin_ff = 1.455 * pow(Re_a_max, -0.656) * pow((tube_adjV / tube_adjH), -0.347) * pow((fin_S / tube_D_o), -0.134) * pow((tube_adjV / tube_D_o), 1.23);
				fin_f = (fin_ff * SA_finned * SA_tot) + (fin_ft * (1 - (SA_finned / SA_tot)) * (1 - (fin_T / SA_finned)));
				fin_DP = (fin_f * SA_tot * pow((rho_a * Vmax_a), 2) / (Ac_fin * 2 * rho_a)) / 1000;
				P_out_a[i][0] = P_in_a[i][0] - fin_DP;
				
				//Refrigerant-side heat transfer
				htc = (1.0 / (tube_k * PI * tube_L)) + (Rc / (PI * tube_D_o * tube_L)) + (Rf / (PI * tube_D_o * tube_L)) + (1.0 / (alfa_a[i][0] * (tube_D_o * PI * (tube_L - (int)(tube_L / fin_S) * fin_T) + (int)(tube_L / fin_S) * 2.0 * (tube_adjV * tube_adjH - PI / 4.0 * tube_D_o * tube_D_o)) * eta_fin));
				
				if (x_in_r[i][0] < 0 || x_in_r[i][0] > 1)		//single-phase heat transfer
				{
					alfa_r[i][0] = 0.023 * k_r * pow(Re_r, 0.8) * pow(Pr_r, 0.3) / (tube_D_i * 1000);

					Q[i][0] = (T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow((alfa_r[i][0] * PI * tube_D_i * tube_L), -1)) + htc), -1));

					mnm.setError(p, Q_conv[i][0], 0); p++;
					mnm.setError(p, Q_boil[i][0], 0); p++;
				}

				else									//two-phase heat transfer
				{
					if (type == 1)						//evaporation
					{
						if (ref_cor == 1)				//Kandlikar
						{
							mass_flux = 4 * G_r[i][0] / (PI * pow(tube_D_i, 2));
							alfa_lo = 0.023 * k_rl * (pow(((mass_flux * tube_D_i * (1 - x_in_r[i][0])) / myu_rl), 0.8)) * (pow((cp_rl * myu_rl / k_rl), 0.4)) / (tube_D_i * 1000);			//in kJ/kg

							Co = (pow(((1 - x_in_r[i][0]) / x_in_r[i][0]), 0.8)) * (pow((rho_rv / rho_rl), 0.5));
							Fr = pow(mass_flux, 2) / (9.81 * tube_D_i * (pow(rho_rl, 2)));
							double C;

							if (Fr > 0.04) { C = 0; }
							else { C = 0.3; }

							mnm.setError(p, Q_conv[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow((((1.1360 * (pow(Co, -0.9)) * (pow((25 * pow(mass_flux, 2) / (9.81 * tube_D_i * pow(rho_rl, 2))), 0.3))) + (667.2 * Ffl * pow((Q_conv[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 0.7))) * alfa_lo * PI * tube_D_i * tube_L), -1)) + htc), -1))));
							p++;

							mnm.setError(p, Q_boil[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow((((0.6683 * (pow(Co, -0.2)) * (pow((25 * pow(mass_flux, 2) / (9.81 * tube_D_i * pow(rho_rl, 2))), 0.3))) + (1058.0 * Ffl * pow((Q_boil[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 0.7))) * alfa_lo * PI * tube_D_i * tube_L), -1)) + htc), -1))));
							p++;

							if (Q_conv[i][0] > Q_boil[i][0])
							{
								Q_boil[i][0] = Q_conv[i][0];
								Q[i][0] = Q_conv[i][0];
								alfa_r[i][0] = ((1.1360 * (pow(Co, -0.9)) * (pow((25 * pow(mass_flux, 2) / (9.81 * tube_D_i * pow(rho_rl, 2))), 0.3))) + (667.2 * Ffl * pow((Q[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 0.7))) * alfa_lo;
							}

							else
							{
								Q_conv[i][0] = Q_boil[i][0];
								Q[i][0] = Q_boil[i][0];
								alfa_r[i][0] = ((0.6683 * (pow(Co, -0.2)) * (pow((25 * pow(mass_flux, 2) / (9.81 * tube_D_i * pow(rho_rl, 2))), 0.3))) + (1058.0 * Ffl * pow((Q[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 0.7))) * alfa_lo;
							}							
						}

						else if (ref_cor == 2)			//Gungor-Winterton
						{
							myu_r = (x_in_r[i][0] * myu_rv) + ((1 - x_in_r[i][0]) * myu_rl);
							VF = 1 / (1 + (((1 - x_in_r[i][0]) / x_in_r[i][0]) * (rho_rv / rho_rl)));
							rho_r = (VF * rho_rv) + ((1 - VF) * rho_rl);

							Xtt = (pow(((1 - x_in_r[i][0]) / x_in_r[i][0]), 0.9)) * (pow((rho_rv / rho_rl), 0.5)) * (pow((myu_rl / myu_rv), 0.1));
							mass_flux = 4 * G_r[i][0] / (PI * pow(tube_D_i, 2));
							Fr = pow(mass_flux, 2) / (9.81 * tube_D_i * pow(rho_rl, 2));
							
							double Re_l = mass_flux * tube_D_i * (1 - x_in_r[i][0]) / myu_rl;
							double P_red = P_in_r[i][0] / Pcrit;																																				// in g/mol
							double alfa_lo = 0.023 * k_rl * (pow(((mass_flux * tube_D_i * (1 - x_in_r[i][0])) / myu_rl), 0.8)) * (pow((cp_rl * myu_rl / k_rl), 0.4)) / (tube_D_i * 1000);			//in kJ/kg

							if (Fr < 0.05)
							{
								mnm.setError(p, Q_conv[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow((((alfa_lo * (pow(Fr, (0.1 - (2 * Fr)))) * (1 + (24000 * (pow((Q_conv[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 1.16))) + (1.37 * (pow((1 / Xtt), 0.86))))) + ((55 * (pow(P_red, 0.12)) * (pow((-log10(P_red)), -0.55)) * (pow(MW, -0.5)) * (pow((Q_conv[i][0] / (PI * tube_D_i * tube_L)), 2 / 3))) * ((pow(Fr, 0.5)) / (1 + (0.0000015 * (pow(Re_l, 1.17)) * (pow((1 + (24000 * (pow((Q_conv[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 1.16))) + (1.37 * (pow((1 / Xtt), 0.86)))), 2))))))) * PI * tube_D_i * tube_L), -1)) + htc), -1))));
								p++;

								mnm.setError(p, Q_boil[i][0], Q_conv[i][0]);
								p++;

								Q[i][0] = Q_conv[i][0];
								alfa_r[i][0] = (alfa_lo * (pow(Fr, (0.1 - (2 * Fr)))) * (1 + (24000 * (pow((Q_conv[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 1.16))) + (1.37 * (pow((1 / Xtt), 0.86))))) + ((55 * (pow(P_red, 0.12)) * (pow((-log10(P_red)), -0.55)) * (pow(MW, -0.5)) * (pow((Q_conv[i][0] / (PI * tube_D_i * tube_L)), 2 / 3))) * ((pow(Fr, 0.5)) / (1 + (0.0000015 * (pow(Re_l, 1.17)) * (pow((1 + (24000 * (pow((Q_conv[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 1.16))) + (1.37 * (pow((1 / Xtt), 0.86)))), 2))))));
							}

							else
							{
								mnm.setError(p, Q_conv[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow((((alfa_lo * (1 + (24000 * (pow((Q_conv[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 1.16))) + (1.37 * (pow((1 / Xtt), 0.86))))) + ((55 * (pow(P_red, 0.12)) * (pow((-log10(P_red)), -0.55)) * (pow(MW, -0.5)) * (pow((Q_conv[i][0] / (PI * tube_D_i * tube_L)), 2 / 3))) * (1 / (1 + (0.0000015 * (pow(Re_l, 1.17)) * (pow((1 + (24000 * (pow((Q_conv[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 1.16))) + (1.37 * (pow((1 / Xtt), 0.86)))), 2))))))) * PI * tube_D_i * tube_L), -1)) + htc), -1))));
								p++;

								mnm.setError(p, Q_boil[i][0], Q_conv[i][0]);
								p++;

								Q[i][0] = Q_conv[i][0];
								alfa_r[i][0] = (alfa_lo * (1 + (24000 * (pow((Q_conv[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 1.16))) + (1.37 * (pow((1 / Xtt), 0.86))))) + ((55 * (pow(P_red, 0.12)) * (pow((-log10(P_red)), -0.55)) * (pow(MW, -0.5)) * (pow((Q_conv[i][0] / (PI * tube_D_i * tube_L)), 2 / 3))) * (1 / (1 + (0.0000015 * (pow(Re_l, 1.17)) * (pow((1 + (24000 * (pow((Q_conv[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_D_i * tube_L)), 1.16))) + (1.37 * (pow((1 / Xtt), 0.86)))), 2))))));
							}
						}

						else if (ref_cor == 3)			//Liu-Winterton
						{
							myu_r = (x_in_r[i][0] * myu_rv) + ((1 - x_in_r[i][0]) * myu_rl);
							VF = 1 / (1 + (((1 - x_in_r[i][0]) / x_in_r[i][0]) * (rho_rv / rho_rl)));
							rho_r = (VF * rho_rv) + ((1 - VF) * rho_rl);
							Pr_r = myu_rl * cp_rl / k_rl;

							mass_flux = 4 * G_r[i][0] / (PI * pow(tube_D_i, 2));
							double Re_l = mass_flux * tube_D_i * (1 - x_in_r[i][0]) / myu_rl;
							double P_red = P_in_r[i][0] / Pcrit;
																																											// in g/mol
							double F = pow((1 + (x_in_r[i][0] * Pr_r * ((rho_rl / rho_rv) - 1))), 0.35);
							double S = 1 / (1 + (0.055 * (pow(F, 0.1)) * (pow(Re_l, 0.16))));
							double Fr = pow(mass_flux, 2) / (9.81 * tube_D_i * pow(rho_rl, 2));
							double alfa_lo = 0.023 * k_rl * (pow(((mass_flux * tube_D_i * (1 - x_in_r[i][0])) / myu_rl), 0.8)) * (pow((cp_rl * myu_rl / k_rl), 0.4)) / (tube_D_i * 1000);			//in kJ/kg

							if (Fr < 0.05)
							{
								mnm.setError(p, Q_conv[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow(((pow(((pow((F * (pow(Fr, (0.1 - (2 * Fr)))) * alfa_lo), 2)) + (pow((S * (pow(Fr, 0.5)) * (55 * (pow(P_red, 0.12)) * (pow((-log10(P_red)), -0.55)) * (pow(MW, -0.5)) * (pow((Q_conv[i][0] / (PI* tube_D_i * tube_L)), (2 / 3))))), 2))), 0.5)) * PI * tube_D_i * tube_L), -1)) + htc), -1))));
								p++;

								mnm.setError(p, Q_boil[i][0], Q_conv[i][0]);
								p++;

								Q[i][0] = Q_conv[i][0];
								alfa_r[i][0] = pow(((pow((F * (pow(Fr, (0.1 - (2 * Fr)))) * alfa_lo), 2)) + (pow((S * (pow(Fr, 0.5)) * (55 * (pow(P_red, 0.12)) * (pow((-log10(P_red)), -0.55)) * (pow(MW, -0.5)) * (pow((Q_conv[i][0] / (PI* tube_D_i * tube_L)), (2 / 3))))), 2))), 0.5);
							}

							else
							{
								mnm.setError(p, Q_conv[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow(((pow(((pow((F * alfa_lo), 2)) + (pow((S * (55 * (pow(P_red, 0.12)) * (pow((-log10(P_red)), -0.55)) * (pow(MW, -0.5)) * (pow((Q_conv[i][0] / (PI* tube_D_i * tube_L)), (2 / 3))))), 2))), 0.5)) * PI * tube_D_i * tube_L), -1)) + htc), -1))));
								p++;

								mnm.setError(p, Q_boil[i][0], Q_conv[i][0]);
								p++;

								Q[i][0] = Q_conv[i][0];
								alfa_r[i][0] = pow(((pow((F * alfa_lo), 2)) + (pow((S * (55 * (pow(P_red, 0.12)) * (pow((-log10(P_red)), -0.55)) * (pow(MW, -0.5)) * (pow((Q_conv[i][0] / (PI* tube_D_i * tube_L)), (2 / 3))))), 2))), 0.5);
							}
						}

						else if (ref_cor == 4)			//Shah
						{
							mass_flux = 4 * G_r[i][0] / (PI * pow(tube_D_i, 2));
							alfa_lo = 0.023 * k_rl * (pow(((mass_flux * tube_D_i * (1 - x_in_r[i][0])) / myu_rl), 0.8)) * (pow((cp_rl * myu_rl / k_rl), 0.4)) / (tube_D_i * 1000);			//in kJ/kg

							Co = (pow(((1 - x_in_r[i][0]) / x_in_r[i][0]), 0.8)) * (pow((rho_rv / rho_rl), 0.5));
							Fr = pow(mass_flux, 2) / (9.81 * tube_D_i * (pow(rho_rl, 2)));
							double N, C, alfa_conv, alfa_boil;

							if (Fr > 0.04) { N = Co; }
							else { N = 0.38 * pow(Fr, -0.3) * Co; }

							double Bo = Q[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_L * tube_D_i);
							if (Bo < 0.0011) { C = 15.43; }
							else { C = 14.7; }

							mnm.setError(p, Q_conv[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow((1.8 * pow(N, -0.8) * alfa_lo * PI * tube_D_i * tube_L), -1)) + htc), -1))));
							p++;

							alfa_conv = 1.8 * pow(N, -0.8) * alfa_lo;

							if (N <= 0.1)
							{
								mnm.setError(p, Q_boil[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow((C * pow((Q_boil[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_L * tube_D_i)), 0.5) * exp(2.47 * pow(N, -0.15)) * alfa_lo * PI * tube_D_i * tube_L), -1)) + htc), -1))));
								p++;

								alfa_boil = C * pow((Q_boil[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_L * tube_D_i)), 0.5) * exp(2.47 * pow(N, -0.15)) * alfa_lo;
							}

							else if (N > 1)
							{
								if (Bo < 0.00003)
								{
									mnm.setError(p, Q_boil[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow(((1 + (46 * pow((Q_boil[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_L * tube_D_i)), 0.5))) *alfa_lo * PI * tube_D_i * tube_L), -1)) + htc), -1))));
									p++;

									alfa_boil = (1 + (46 * pow((Q_boil[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_L * tube_D_i)), 0.5))) *alfa_lo;
								}

								else
								{
									mnm.setError(p, Q_boil[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow(230 * (pow((Q_boil[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_L * tube_D_i)), 0.5) *alfa_lo * PI * tube_D_i * tube_L), -1)) + htc), -1))));
									p++;

									alfa_boil = 230 * pow((Q_boil[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_L * tube_D_i)), 0.5) *alfa_lo;
								}
							}

							else
							{
								mnm.setError(p, Q_boil[i][0], ((T_in_a[i][0] - T_in_r[i][0]) * (pow(((pow((C * pow((Q_boil[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_L * tube_D_i)), 0.5) * exp(2.74 * pow(N, -0.1)) * alfa_lo * PI * tube_D_i * tube_L), -1)) + htc), -1))));
								p++;

								alfa_boil = C * pow((Q_boil[i][0] / (mass_flux * (h_rv - h_rl) * PI * tube_L * tube_D_i)), 0.5) * exp(2.74 * pow(N, -0.1)) * alfa_lo;
							}

							if (Q_conv[i][0] > Q_boil[i][0])
							{
								Q_boil[i][0] = Q_conv[i][0];
								Q[i][0] = Q_conv[i][0];
								alfa_r[i][0] = alfa_conv;
							}

							else
							{
								Q_conv[i][0] = Q_boil[i][0];
								Q[i][0] = Q_boil[i][0];
								alfa_r[i][0] = alfa_boil;
							}
						}
						
						else {}

						//Outlet parameter values
						h_out_r[i][0] = h_in_r[i][0] + (Q[i][0] / G_r[i][0]);

						refr.state_ph(P_out_r[i][0], h_out_r[i][0]);
						x_out_r[i][0] = refr.Rc.x;
						T_out_r[i][0] = refr.Rc.T;

						h_out_a[i][0] = h_in_a[i][0] - (Q[i][0] / G_a[i][0]);

						refa.state_ph(P_out_a[i][0], h_out_a[i][0]);
						T_out_a[i][0] = refa.Rc.T;
					}

					else								//condensation
					{
						//in progress
						mnm.setError(p, Q_conv[i][0], 0); p++;

						mnm.setError(p, Q_boil[i][0], 0); p++;

						Q[i][0] = 0;

						//Outlet parameter values
						h_out_r[i][0] = h_in_r[i][0] - (Q[i][0] / G_r[i][0]);

						refr.state_ph(P_out_r[i][0], h_out_r[i][0]);
						x_out_r[i][0] = refr.Rc.x;
						T_out_r[i][0] = refr.Rc.T;

						h_out_a[i][0] = h_in_a[i][0] + (Q[i][0] / G_a[i][0]);

						refa.state_ph(P_out_a[i][0], h_out_a[i][0]);
						T_out_a[i][0] = refa.Rc.T;
					}
				}

				//Outlet parameter values single phase
				h_out_r[i][0] = h_in_r[i][0] + (Q[i][0] / G_r[i][0]);

				refr.state_ph(P_out_r[i][0], h_out_r[i][0]);
				x_out_r[i][0] = refr.Rc.x;
				T_out_r[i][0] = refr.Rc.T;

				h_out_a[i][0] = h_in_a[i][0] - (Q[i][0] / G_a[i][0]);

				refa.state_ph(P_out_a[i][0], h_out_a[i][0]);
				T_out_a[i][0] = refa.Rc.T;
			}

			mnm.prt_sum();
		}

		//Input to the next tube/column					//observes the downstream (destination) tube(s)
		for (int a = 1; a <= Ntubes; a++)
		{
			int i = BFS_path[a];

			//Air-side
			if (i <= Ntubes - Nrow)
			{
				T_in_a[i + Nrow][0] = T_out_a[i][0];
				P_in_a[i + Nrow][0] = P_out_a[i][0];
			}
			else {}

			//Refrigerant-side
			//Mass and energy (enthalpy) balance at the junction
			G_res = 0.0;
			H_res = 0.0;

			for (int j = 1; j <= Ntubes; j++)
			{
				G_res = G_r[i][0];
				H_res = G_res * h_out_r[i][0];

				if (TTCM[i][j] == 1)								//Downstream flow
				{
					if (row_neg[j] > 1)								//Connected tube has a merge
					{
						for (int k = 1; k <= Ntubes; k++)
						{
							if (TTCM[j][k] == -1 && k != i)
							{
								G_res = G_res + G_r[k][0];
								H_res = H_res + (G_r[k][0] * h_out_r[k][0]);
							}
						}
					}
				}

				//Assigning to connected tube
				for (int j = 1; j <= Ntubes; j++)
				{
					if (TTCM[i][j] == 1)
					{
						h_in_r[j][0] = H_res / G_res;
						P_in_r[j][0] = P_out_r[i][0];

						refr.state_ph(P_in_r[j][0], h_in_r[j][0]);
						T_in_r[j][0] = refr.Rc.T;
					}
				}
			}
		}
	}

	//End of Newton-Raphson loop//
	end_time = clock();
	run_time[0] = (double)(end_time - start_time) / 1000;			//[sec]

	return 0;
}   */
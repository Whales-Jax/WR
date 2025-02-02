/*
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

#include "CFluidParameter.h"
#include "CNewtonRaphsonMethod.h"
#include "DXProperty_ver06.h"

int main()
{
	//Result printout
	ofstream complex("results_TTM.csv");

	DXProperty refr;
	DXProperty refa;

	refr.reibai = "WATER.FLD";
	refa.reibai = "AIR.PPF";
	refr.joutai = "DFT";
	refa.joutai = "IIR";
	refr.LoadDLL("refprop.DLL");
	refa.LoadDLL("refprop2.DLL");
	refr.setup();
	refa.setup();

	//Refrigerant properties
	double T_in_r[50];						//refrigerant local temperature in [degC]
	double T_out_r[50];						//refrigerant local temperature out [deg C]
	double P_in_r[50];						//refrigerant local pressure in [kPa]
	double P_out_r[50];						//refrigerant local pressure out [kPa]
	double G_r[50];							//refrigerant local mass flowrate [kg/s]
	double h_in_r[50];						//refrigerant local specific enthalpy in [kJ/kg]
	double h_out_r[50];						//refrigerant local specific enthalpy out [kJ/kg]

	double T_in_r1;							//refrigerant inlet temperature [degC]
	double T_out_r1;						//refrigerant outlet temperature [degC]
	double P_in_r1;							//refrigerant inlet pressure [kPa]
	double P_out_r1;						//refrigerant outlet pressure [kPa]
	double G_r1;							//refrigerant total mass flowrate [kg/s]
	double h_in_r1;							//refrigerant inlet enthalpy [kJ/kg]
	double h_out_r1;						//refrigerant outlet enthalpy [kJ/kg]

	//Air properties
	double T_in_air[50];					//air temperature in [degC]
	double T_out_air[50];					//air temperature out [degC]
	double P_in_air[50];					//air pressure in [kPa]
	double P_out_air[50];					//air pressure out [kPa]
	double G_air[50];						//air mass flowrate [kg/s]
	double h_in_air[50];					//air specific enthalpy in [kJ/kg]
	double h_out_air[50];					//air specific enthalpy out [kJ/kg]
	double V_air[50];						//air velocity [m/s]

	double G_air1;							//air total mass flowrate [kJ/kg]
	double T_air1;							//air inlet temperature [degC]
	double V_air1;							//air inlet velocity [m/s]
	double P_air1 = 101.325;				//air inlet pressure [kPa]

	double eta_fin[50];
	double Re_r[50];
	double Re_a[50], Re_a_max[50];
	double Prandtl_r[50];
	double Prandtl_a[50];
	double alfa_r[50], alfa_a[50], KA[50], Q[50], f[50];

	double tube[50];


	int i, j, p, t, k, J, K, O, x, y;
	double Rf, Rc, myu_r, k_a, k_r, P_r, T_r, P_a, T_a, DP, DP_i, DP_o, DT, DH, Qtot, Gtot, Atot, A_f_a;
	double D_H_air, H_res, G_res, Vmax_a, Sd;
	double jcf, jcf_t;													//Colburn factor
	double fin_ff, fin_ft, fin_dP_f, fin_dP_t, fin_dP[50], fin_f;		//air-side pressure drop variables
	double SA_tot, SA_finned, Ac_fin;

	int sumcol[50], sumrow[50];
	int tube_mtrx[100][100];


	//Heat exchanger geometry
	const double tube_D_o = 0.01;								//tube outer diameter [m]
	const double tube_D_i = 0.009;								//tube inner diameter [m]
	const double tube_L = 1.0;									//tube length [m]
	const double tube_T = 0.0005;								//tube thickness [m]
	const double tube_adjH = 0.03;								//horizontal distance between adjacent tubes [m]
	const double tube_adjV = 0.03;								//vertical distance between adjacent tubes [m]
	const double fin_S = 0.003;									//fin spacing [m]
	const double fin_T = 0.00013;								//fin thickness [m]
	const double PI = (6 * (asin(0.5)));						//pi
	const double beta = 0.0;									//tube inclination [rad]
	const double ktube = 205.0 / 1000.0;						//aluminum thermal conductivity [kW/m-K]

	D_H_air = (fin_S * tube_adjV * 2) / (fin_S + tube_adjV);	//hydraulic diameter

	int Ntubes, Njunc, Nrow, Ncol;

	SA_finned = (int)(tube_L / fin_S) * 2 * tube_adjH * tube_adjV;
	Ac_fin = (tube_L - (((int)(tube_L / fin_S))*fin_T)) * tube_adjV;
	SA_tot = SA_finned + ((tube_L - (((int)(tube_L / fin_S))* fin_T)) * PI * tube_D_o);

	//Tube-tube connectivity matrix (TTM)//
	cout << "Enter the number of rows and columns, separated by space." << endl;
	cin >> Nrow >> Ncol;

	Ntubes = Nrow * Ncol;
	Njunc = Ntubes - 1;						//no splits, no merges

	//TTM initialization
	for (i = 1; i <= Ntubes; i++)
	{
		for (j = 1; j <= Ntubes; j++)
		{
			tube_mtrx[i][j] = 0;
		}
	}

	//Input of Tube Connectivity
	for (i = 1; i <= Njunc; i++)
	{
		cout << "\nEnter the tubes connected (Note: source destination)" << endl;
		cin >> x >> y;

		tube_mtrx[x][y] = 1;
		tube_mtrx[y][x] = -1;
	}

	cout << "\n";

	//Search for inlet and outlet tubes, and splitting and merging
	int c_pos[50], c_neg[50];

	for (i = 1; i <= Ntubes; i++)
	{
		c_pos[i] = 0;
		c_neg[i] = 0;

		for (j = 1; j <= Ntubes; j++)
		{
			if (tube_mtrx[i][j] == 1)
			{
				c_pos[i] = c_pos[i] + 1;
			}

			if (tube_mtrx[i][j] == -1)
			{
				c_neg[i] = c_neg[i] + 1;
			}
		}
				
		if (c_neg[i] == 0) cout << "The input tube is at tube [" << i << "]." << endl;
		if (c_pos[i] == 0) cout << "The output tube is at tube [" << i << "]." << endl;

		if (c_neg[i] > 1) cout << "Merging at tube [" << i << "]." << endl;
		if (c_pos[i] > 1) cout << "Splitting at tube [" << i << "]." << endl;
	}

	clock_t start_time, end_time;
	start_time = clock();

	//Input conditions
	P_in_r1 = 102.325;							//refrigerant inlet pressure
	P_out_r1 = 101.325;							//refrigerant outlet pressure
	h_in_r1 = 84.0;								//refrigerant inlet enthalpy
	T_air1 = 35.0;								//air inlet temperature

	refr.state_ph(P_in_r1, h_in_r1);
	T_in_r1 = refr.Rc.T;						//refrigerant inlet temperature

	G_air1 = 0.27;								//air mass flowrate
	G_r1 = 0.015;								//refrigerant mass flowrate

	refa.state_tp(T_air1, P_air1);
	double rho_a1 = refa.Rc.rho;
	V_air1 = G_air1 / ((tube_L - (int)(tube_L / fin_S) * fin_T) * tube_adjV *double(Nrow)* rho_a1);

	Sd = pow(((pow(tube_adjH, 0.5)) + (pow((tube_adjV / 2), 0.5))), 0.5);

	if (Sd < ((tube_adjV + tube_D_o) / 2))
	{
		Vmax_a = (tube_adjV / (2 * (Sd - tube_D_o))) * V_air1;
	}

	else
	{
		Vmax_a = (tube_adjV / (tube_adjV - tube_D_o)) * V_air1;
	}

	//Local variables initialization
	for (i = 1; i <= Ntubes; i++)
	{
		T_in_air[i] = T_air1;
		T_out_air[i] = T_air1;
		P_in_air[i] = P_air1;
		G_air[i] = G_air1 / double(Nrow);
		G_r[i] = G_r1;

		refr.state_ph(P_in_r1, h_in_r1);
		double myu_r = refr.Rc.visc;
		double rho_r = refr.Rc.rho;

		DP = -128 * tube_L*G_r[i] * myu_r / (PI*rho_r*pow(tube_D_i, 4.0) * 1000000);			//in Pa

		if (c_pos[i] == 0)
		{
			O = i;
		}

		else
		{
			if (c_neg[i] == 0)
			{
				P_in_r[i] = P_in_r1;
				P_out_r[i] = P_in_r1 + DP / 1000.0;
				K = i;

				for (j = 1; j <= Ntubes; j++)
				{
					for (J = 1; J <= Ntubes; J++)
					{
						if (tube_mtrx[K][J] == 1)
						{
							P_in_r[J] = P_out_r[K];
							P_out_r[J] = P_in_r[J] + DP / 1000.0;
							K = J;
							cout << "\ni = " << i << "\tj = " << j << "\tJ = " << J << "\tK = " << K;
							break;
						}
						else {}
					}
				}

			}
			else {}
		}

		h_in_r[i] = h_in_r1;
		h_out_r[i] = h_in_r1;
		T_out_r[i] = T_in_r1;
		T_in_r[i] = T_in_r1;
	}

	for (int j = 1; j <= Njunc; j++)
	{
		if (sumrow[j] > 0)
		{
			for (i = 1; i <= Ntubes; i++)
			{
				if (tube_mtrx[j][i] == 1)
				{
					G_r[i] = G_r1 / double(sumrow[j] + 1);

					for (k = j + 1; k <= Njunc; k++)
					{
						if (tube_mtrx[k][i] == -1)
						{
							for (t = 1; t <= Ntubes; t++)
							{
								if (tube_mtrx[k][t] == 1)
								{
									G_r[t] = G_r1 / double(sumrow[j] + 1);
								}
								else {}

							}
						}
						else {}
					}
				}
				else {}
			}
		}
		else {}
	}

	//Printout of the initialized variables
	for (i = 1; i <= Ntubes; i++)
	{
		cout << "\n\nAt tube_" << i << "\trefrigerant's\tP_in = " << P_in_r[i] << "\t\tP_out = " << P_out_r[i] << "\t\th_in = " << h_in_r[i] << "\t\tT_out = " << T_out_r[i] << "\t\tG_r = " << G_r[i];
		cout << "\n\t\tair's\t\tP_in = " << P_in_air[i] << "\t\tT_out = " << T_out_air[i] << "\t\tT_in = " << T_in_air[i];
	}
	cout << endl;
	
	//Newton-Raphson Method//
	CNewtonRaphsonMethod mnm;
	mnm.setup(2000, 1 + 1e-12, 1e-4);

	p = 0;

	for (i = 1; i <= Ntubes; i++)
	{
		mnm.setValue(p, G_r[i]);
		p++;
		mnm.setValue(p, P_out_r[i]);
		p++;
		mnm.setValue(p, T_out_r[i]);
		p++;
		mnm.setValue(p, T_out_air[i]);
		p++;
	}

	mnm.setAcc(0.2);
	mnm.initial();

	for (mnm.main_loop_init(); mnm.main_loop_check(); mnm.main_loop_reinit())
	{
		for (mnm.sub_loop_init(); mnm.sub_loop_check(); mnm.sub_loop_reinit())
		{
			p = 0;

			for (i = 1; i <= Ntubes; i++)
			{
				G_r[i] = mnm.getValue(p);
				p++;
				P_out_r[i] = mnm.getValue(p);
				p++;
				T_out_r[i] = mnm.getValue(p);
				p++;
				T_out_air[i] = mnm.getValue(p);
				p++;
			}

			p = 0;

			//Fundamental equations, tube loop
			for (i = 1; i <= Ntubes; i++)
			{
				//Fin efficiency
				eta_fin[i] = 0.77;

				//Fouling resistance factor
				Rf = 0.0;

				//Fin-contact resistance factor
				Rc = 0.0;

				//Refrigerant- and air-heat transfer coefficients
				refr.state_tp(T_in_r[i], P_in_r[i]);
				refa.state_tp(T_in_air[i], P_in_air[i]);

				double myu_r = refr.Rc.visc / 1000000;
				double rho_a = refa.Rc.rho;
				double myu_a = refa.Rc.visc / 1000000;


				Re_r[i] = G_r[i] * 4.0 / PI / tube_D_i / myu_r;
				Re_a[i] = rho_a * V_air1 * tube_D_o / myu_a;
				Re_a_max[i] = rho_a * Vmax_a * tube_D_o / myu_a;

				double k_a = refa.Rc.thc;
				double cp_a = refa.Rc.cp * 1000;
				double k_r = refr.Rc.thc;
				double cp_r = refr.Rc.cp * 1000;

				Prandtl_a[i] = myu_a * cp_a / k_a;
				Prandtl_r[i] = myu_r * cp_r / k_r;

				jcf = (rho_a*Vmax_a*cp_a) * pow(Prandtl_a[i], 1.5) * 0.163 * pow(Re_a_max[i], -0.369) * pow((tube_adjV / tube_adjH), 0.106) * pow((fin_S / tube_D_o), 0.0138) * pow((tube_adjV / tube_D_o), 0.13);

				//Water-side heat transfer coefficient from Dittus-Boelter	
				alfa_r[i] = 0.023 * k_r * pow(Re_r[i], 0.8) * pow(Prandtl_r[i], 0.3) / (1000 * tube_D_i);

				//Air-side heat transfer coefficient from Kim-Youn-Webb
				//alfa_a[i] = k_a / tube_D_o * 0.097 * pow(Re_a[i], 0.671) / 1000;

				if (Nrow >= 3)
				{
					jcf_t = jcf;
				}

				else
				{
					jcf_t = 1.043 * pow(0.163 * pow(Re_a_max[i], -0.14) * pow((tube_adjV / tube_adjH), -0.564) * pow((fin_S / tube_D_o), -0.123) * pow((tube_adjV / tube_D_o), 1.17), 3 - Nrow) / jcf;
				}

				alfa_a[i] = jcf_t / 1000;

				//Air-side pressure drop
				fin_ft = (4 / PI) * (0.25 + (0.118 * pow(Re_a_max[i], -0.16) / (pow((tube_adjV / tube_D_o) - 1, 1.08)))) * ((tube_adjV / tube_D_o) - 1);
				fin_ff = 1.455 * pow(Re_a_max[i], -0.656) * pow((tube_adjV / tube_adjH), -0.347) * pow((fin_S / tube_D_o), -0.134) * pow((tube_adjV / tube_D_o), 1.23);
				fin_f = (fin_ff * SA_finned * SA_tot) + (fin_ft * (1 - (SA_finned / SA_tot)) * (1 - (fin_T / SA_finned)));
				fin_dP[i] = (fin_f * SA_tot * pow((rho_a * Vmax_a), 2) / (Ac_fin * 2 * rho_a)) / 1000;
				P_out_air[i] = P_in_air[i] - fin_dP[i];

				//Thermal conductivity and local heat transfer
				KA[i] = pow(1.0 / (alfa_r[i] * PI * tube_D_i * tube_L) + (1.0) / (ktube * (PI * tube_L)) + Rc / (PI * tube_D_o * tube_L) + Rf / (PI * tube_D_o * tube_L) + 1.0 / (alfa_a[i] * (tube_D_o * PI * (tube_L - (int)(tube_L / fin_S) * fin_T) + (int)(tube_L / fin_S) * 2.0 * (tube_adjV * tube_adjH - PI / 4.0 * tube_D_o * tube_D_o)) * eta_fin[i]), -1.0);
				Q[i] = KA[i] * (T_in_air[i] - T_in_r[i]);

				//Equations at the outlet
				if (sumcol[i] == 100)
				{
					f[i] = 64.0 * Re_r[i];

					//Hagen-Poisuelle equation
					mnm.setError(p, G_r[i], G_r1);
					p++;

					//Hagen-Poisuelle equation
					double rho_r = refr.Rc.rho;
					mnm.setError(p, -128.0 * tube_L*G_r[i] * myu_r / (PI * rho_r * pow(tube_D_i, 4.0)), (P_out_r[i] - P_in_r[i]) * 1000.0);
					p++;

					//cout << "\nTube " << i << "\trho_r = " << rho_r;
				}

				else
				{
					//Residual mass flow rate at junction j
					G_res = 0.0;

					for (j = 1; j <= Ntubes; j++)
					{
						if (tube_mtrx[i][j] == -1)
						{
							G_res = G_r[j];

							if (c_pos[j] > 1)
							{
								for (int t = 1; t <= Ntubes; t++)
								{
									if (tube_mtrx[j][t] == 1)
									{
										if (t == i)
										{
										}

										else
										{
											G_res = G_res - G_r[t];
										}
									}

									else {}
								}
							}

							else {}
						}

						else
						{
							if (c_neg[i] == 0)
							{
								G_res = G_r1;
							}

							else {}
						}
					}

					mnm.setError(p, G_r[i], G_res);
					p++;

					f[i] = 64.0 * Re_r[i];

					//Hagen-Poiseulle equation
					double rho_r = refr.Rc.rho;
					mnm.setError(p, -128.0 * tube_L * G_r[i] * myu_r / (PI * rho_r * pow(tube_D_i, 4.0)), (P_out_r[i] - P_in_r[i]) * 1000.0);
					p++;

					h_out_r[i] = h_in_r[i] + Q[i] / G_r[i];

					refr.state_ph(P_out_r[i], h_in_r[i] + Q[i] / G_r[i]);
					double T_in_rx = refr.Rc.T;

					mnm.setError(p, T_out_r[i], T_in_rx);
					p++;

					refa.state_tp(T_in_air[i], P_in_air[i]);
					double cp_ax = refa.Rc.cp;

					DT = Q[i] / G_air[i] / cp_ax;

					mnm.setError(p, T_out_air[i], T_in_air[i] - DT);
					p++;

					//cout << "\nTube " << i << "\th_out = " << h_out_r[i] << "\tDT = " << DT;
				}
			}
			mnm.prt_sum();
		}

		//Input to next tube or next column
		for (i = 1; i <= Ntubes; i++)
		{

			//Output
			if (sumcol[i] == 100)
			{
				P_out_r[i] = P_out_r1;

				if (i <= Ntubes - Nrow)
				{
					T_in_air[i + Nrow] = T_out_air[i];
					P_in_air[i + Nrow] = P_out_air[i];
				}

				else {}
			}

			//Generic tube
			else
			{
				if (i <= Ntubes - Nrow)
				{
					T_in_air[i + Nrow] = T_out_air[i];
					P_in_air[i + Nrow] = P_out_air[i];
				}
				else {}

				//Mass and enthalpy balance at the junction
				if (c_neg[i] == 0)
				{
					H_res = G_r[i] * h_out_r[i];
					G_res = G_r[i];

					for (j = 1; j <= Ntubes; j++)
					{
						if (tube_mtrx[i][j] == 1 && c_neg[j] > 1)
						{
							for (int t = 1; t <= Ntubes; t++)
							{
								if (t == i) {}

								else
								{
									H_res = H_res + (G_r[t] * h_out_r[t]);
									G_res = G_res + G_r[t];
								}
							}
						}

						else {}
					}

					for (j = 1; j <= Ntubes; j++)
					{
						if (tube_mtrx[i][j] == 1)
						{
							h_in_r[j] = H_res / G_res;
							P_in_r[j] = P_out_r[i];

							refr.state_ph(P_in_r[j], h_in_r[j]);
							double T_rt = refr.Rc.T;

							T_in_r[j] = T_rt;
						}
						else {}
					}
				}

				else
				{
					H_res = G_r[i] * h_out_r[i];
					G_res = G_r[i];

					for (j = 1; j <= Ntubes; j++)
					{
						if (tube_mtrx[i][j] == 1 && c_neg[j] > 1)
						{
							for (int t = 1; t <= Ntubes; t++)
							{
								if (t == i) {}

								else
								{
									H_res = H_res + (G_r[t] * h_out_r[t]);
									G_res = G_res + G_r[t];
								}
							}
						}

						else {}
					}

					for (j = 1; j <= Ntubes; j++)
					{
						if (tube_mtrx[i][j] == 1)
						{
							h_in_r[j] = H_res / G_res;
							P_in_r[j] = P_out_r[i];

							refr.state_ph(P_in_r[j], h_in_r[j]);
							double T_rt = refr.Rc.T;

							T_in_r[j] = T_rt;
						}
						else {}
					}
				}
			}
		}
	}

	//End of Newton-Raphson loop//

	//Calculation of total outputs
	K = 0;
	DP = 0.0;
	Qtot = 0.0;
	Gtot = 0.0;
	DP_i = 0.0;
	DP_o = 0.0;

	for (i = 1; i <= Ntubes; i++)
	{
		Qtot += Q[i];

		if (c_pos[i] == 0)
		{
			DP_o = -P_out_r[i];
			Gtot += G_r[i];
		}

		else
		{
			if (c_neg[i] == 0)
			{
				DP_i = P_in_r[i];
			}
			else {}
		}

		cout << "\nTube [" << i << "]. G_r = " << G_r[i] << "\nTube [" << i << "]. P_out_r = " << P_out_r[i] << "\nTube [" << i << "]. T_out_r = " << T_out_r[i] << "\nTube [" << i << "]. P_out_air = " << P_out_air[i] << "\nTube [" << i << "]. T_out_air = " << T_out_air[i] << "\nTube [" << i << "]. P_in_r= " << P_in_r[i] << "\nTube [" << i << "]. T_in_r = " << T_in_r[i] << "\nTube [" << i << "]. P_in_air = " << P_in_air[i] << "\nTube [" << i << "]. T_in_air = " << T_in_air[i];
	}

	Atot = (tube_D_o * PI * (tube_L - (int)(tube_L / fin_S) * fin_T) + (int)(tube_L / fin_S) * 2.0 * (tube_adjV * tube_adjH - PI / 4.0 * tube_D_o * tube_D_o))* double(Ntubes);
	cout << "\nQtot = " << Qtot << " DP = " << DP_i + DP_o << " Gtot = " << Gtot << " Atot = " << Atot;

	end_time = clock();
	double run_time = (double)(end_time - start_time) / 1000; // [sec]

	//Write results on .csv file
	cout << "Total run time: " << run_time << "[sec]\n";

	cout << "Writing to .csv file...\n";

	complex << "total run time: " << run_time << "[sec]\n";
	complex << "Qtot = " << Qtot << "," << " DP = " << DP_i + DP_o << "," << " Gtot = " << Gtot << "," << " Atot = " << Atot << ",";

	double	T_in_r_mtrx[100][100],
		T_out_r_mtrx[100][100],
		P_in_r_mtrx[100][100],
		P_out_r_mtrx[100][100],
		T_in_air_mtrx[100][100],
		T_out_air_mtrx[100][100],
		h_in_r_mtrx[100][100],
		h_out_r_mtrx[100][100],
		G_r_mtrx[100][100];

	int r, c;

	///1
	i = 0;
	for (c = 0; c < Ncol; c++)
	{
		for (r = 0; r < Nrow; r++)
		{
			++i;
			T_in_r_mtrx[r][c] = T_in_r[i];
		}
	}

	complex << "\nInput Temperature\n";
	for (r = 0; r < Nrow; r++)
	{
		for (c = 0; c < Ncol; c++)
		{
			complex << T_in_r_mtrx[r][c] << ",";
		}
		complex << endl;
	}

	///2
	i = 0;
	for (c = 0; c < Ncol; c++)
	{
		for (r = 0; r < Nrow; r++)
		{
			++i;
			T_out_r_mtrx[r][c] = T_out_r[i];
		}
	}

	complex << "\nOutput Temperature\n";
	for (r = 0; r < Nrow; r++)
	{
		for (c = 0; c < Ncol; c++)
		{
			complex << T_out_r_mtrx[r][c] << ",";
		}
		complex << endl;
	}

	///3
	i = 0;
	for (c = 0; c < Ncol; c++)
	{
		for (r = 0; r < Nrow; r++)
		{
			++i;
			P_in_r_mtrx[r][c] = P_in_r[i];
		}
	}

	complex << "\nInput Pressure\n";
	for (r = 0; r < Nrow; r++)
	{
		for (c = 0; c < Ncol; c++)
		{
			complex << P_in_r_mtrx[r][c] << ",";
		}
		complex << endl;
	}

	///4
	i = 0;
	for (c = 0; c < Ncol; c++)
	{
		for (r = 0; r < Nrow; r++)
		{
			++i;
			P_out_r_mtrx[r][c] = P_out_r[i];
		}
	}

	complex << "\nOutput Pressure\n";
	for (r = 0; r < Nrow; r++)
	{
		for (c = 0; c < Ncol; c++)
		{
			complex << P_out_r_mtrx[r][c] << ",";
		}
		complex << endl;
	}

	///5
	i = 0;
	for (c = 0; c < Ncol; c++)
	{
		for (r = 0; r < Nrow; r++)
		{
			++i;
			T_in_air_mtrx[r][c] = T_in_air[i];
		}
	}

	complex << "\nAir Input Temperature\n";
	for (r = 0; r < Nrow; r++)
	{
		for (c = 0; c < Ncol; c++)
		{
			complex << T_in_air_mtrx[r][c] << ",";
		}
		complex << endl;
	}

	///6
	i = 0;
	for (c = 0; c < Ncol; c++)
	{
		for (r = 0; r < Nrow; r++)
		{
			++i;
			T_out_air_mtrx[r][c] = T_out_air[i];
		}
	}

	complex << "\nAir Output Temperature\n";
	for (r = 0; r < Nrow; r++)
	{
		for (c = 0; c < Ncol; c++)
		{
			complex << T_out_air_mtrx[r][c] << ",";
		}
		complex << endl;
	}

	///7
	i = 0;
	for (c = 0; c < Ncol; c++)
	{
		for (r = 0; r < Nrow; r++)
		{
			++i;
			h_in_r_mtrx[r][c] = h_in_r[i];
		}
	}

	complex << "\nInput Enthalpy\n";
	for (r = 0; r < Nrow; r++)
	{
		for (c = 0; c < Ncol; c++)
		{
			complex << h_in_r_mtrx[r][c] << ",";
		}
		complex << endl;
	}

	///8
	i = 0;
	for (c = 0; c < Ncol; c++)
	{
		for (r = 0; r < Nrow; r++)
		{
			++i;
			h_out_r_mtrx[r][c] = h_out_r[i];
		}
	}

	complex << "\nOutput Enthalpy\n";
	for (r = 0; r < Nrow; r++)
	{
		for (c = 0; c < Ncol; c++)
		{
			complex << h_out_r_mtrx[r][c] << ",";
		}
		complex << endl;
	}

	///9
	i = 0;
	for (c = 0; c < Ncol; c++)
	{
		for (r = 0; r < Nrow; r++)
		{
			++i;
			G_r_mtrx[r][c] = G_r[i];
		}
	}

	complex << "\nRefrigerant flowrate\n";
	for (r = 0; r < Nrow; r++)
	{
		for (c = 0; c < Ncol; c++) {
			complex << G_r_mtrx[r][c] << ",";
		}
		complex << endl;
	}

	complex << "\nDisplaying the JUNCTION-TUBE MATRIX\n";
	complex << ",";
	for (i = 1; i <= Ntubes; i++)
	{
		complex << "t(" << i << ")" << ",";
	}

	for (i = 1; i <= Njunc; i++)
	{
		complex << "\nj(" << i << ")" << ",";
		for (j = 1; j <= Ntubes; j++)
		{
			complex << tube_mtrx[i][j] << ",";
		}
	}
	
	system("pause");
	return 0;
}*/
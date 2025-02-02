/*
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

#include "CFluidParameter.h"
#include "CNewtonRaphsonMethod.h"
#include "DXProperty_ver06.h"

int mass[50];														//Array for the path finding where tube numbers are inputted once explored
void mass_initial(int K, int Ntubes);

int main()
{
	ofstream complex("Basic Simulator (results).csv");				//Results printout

	DXProperty ref;
	ref.reibai = "R32.FLD", "AIR.PPF";
	ref.joutai = "DEF";
	ref.ref_state[0] = "DEF";
	ref.ref_state[1] = "IIR";
	ref.LoadDLL("refprop.DLL");
	ref.setmulti();

	//===Refrigerant Properties===//
	double T_in_r[50];													//refrigerant local inlet temperature [degC]
	double T_out_r[50];													//refrigerant local outlet temperature [deg C]
	double P_in_r[50];													//refrigerant local inlet pressure [kPa]
	double P_out_r[50];													//refrigerant local outlet pressure [kPa]
	double G_r[50];														//refrigerant local mass flowrate [kg/s]
	double h_in_r[50];													//refrigerant local inlet specific enthalpy [kJ/kg]
	double h_out_r[50];													//refrigerant local outlet specific enthalpy [kJ/kg]
	double x_in_r[50];													//refrigerant local inlet quality
	double x_out_r[50];													//refrigerant local outlet quality

	double T_in_r1;														//refrigerant inlet temperature [degC]
	double T_out_r1;													//refrigerant outlet temperature [degC]
	double P_in_r1;														//refrigerant inlet pressure [kPa]
	double P_out_r1;													//refrigerant outlet pressure [kPa]
	double G_r1;														//refrigerant total mass flowrate [kg/s]
	double h_in_r1;														//refrigerant inlet enthalpy [kJ/kg]
	double h_out_r1;													//refrigerant outlet enthalpy [kJ/kg]

	//===Air Properties===//
	double T_in_air[50];												//air local inlet temperature [degC]
	double T_out_air[50];												//air local outlet temperature out [degC]
	double P_in_air[50];												//air local inlet pressure [kPa]
	double P_out_air[50];												//air local outlet pressure [kPa]
	double G_air[50];													//air mass flowrate [kg/s]
	double h_in_air[50];												//air local inlet specific enthalpy in [kJ/kg]
	double h_out_air[50];												//air local outlet specific enthalpy [kJ/kg]
	double V_air[50];													//air velocity [m/s]

	double G_air1;														//air total mass flowrate [kJ/kg]
	double T_air1;														//air inlet temperature [degC]
	double V_air1;														//air inlet velocity [m/s]
	double P_air1 = 101.325;											//air inlet pressure [kPa]

	//===Heat Exchanger Geometry===//
	const double tube_D_o = 0.01;										//tube outer diameter [m]
	const double tube_D_i = 0.009;										//tube inner diameter [m]
	const double tube_L = 1.0;											//tube length [m]
	const double tube_T = 0.0005;										//tube thickness [m]
	const double tube_adjH = 0.03;										//horizontal distance between adjacent tubes [m]
	const double tube_adjV = 0.03;										//vertical distance between adjacent tubes [m]
	const double fin_S = 0.003;											//fin spacing [m]
	const double fin_T = 0.00013;										//fin thickness [m]
	const double PI = (6 * (asin(0.5)));								//pi
	const double beta = 0.0;											//tube inclination [rad]
	const double ktube = 205.0 / 1000.0;								//aluminum thermal conductivity [kW/m-K]

	//Other refrigerant-related parameters
	double eta_fin[50], Re_r[50], Prandtl_r[50];
	double myu_r1, k_r1, rho_r1;

	//Other air-related parameters
	double Re_a[50], Re_a_max[50], Prandtl_a[50];
	double k_a, D_H_air, Vmax_a, Sd;
	double jcf, jcf_t;													//Colburn factor
	double fin_ff, fin_ft, fin_dP_f, fin_dP_t, fin_dP[50], fin_f;		//air-side pressure drop variables


	//Other parameters for Newton-Raphson
	int p;																//Counter for Newton-Rahpson variables
	int ref_out = 0;													//Assigned reference outlet tube for circuitry that has multiple outlets
	double alfa_r[50], alfa_a[50], KA[50], Q[50], f[50];				//f[50] is not yet used in this code
	double Rf, Rc, DP, DP_i, DP_o, DT, DH, Qtot, Gtot, Atot;
	double H_res, G_res;
	double SA_tot, SA_finned, Ac_fin;
	
	//Variables for the TTCM
	int tube_mtrx[50][50], c_pos[50], c_neg[50];
	int Nrow, Ncol, Ntubes, Ninlet, Noutlet, Nloop, Npair, Ncircuit;
	int xrow, xcol, yrow, ycol;
	int temp_merge, temp_out;
	int checker[50][50];												//Used for checking the U-bend connectivity
	int mass_done[50];													//Stores tube numbers when the observed tube has been completely explored (all of the connected tubes of the observed tube are inputted in the dummy array)
	int loop[25];														//Stores the tube number of the start of each circuit loop based from the mass_done array
	int simple[500];													//Simplified array of the TTCM
	int random[50];														//Stores the randomized element of the simplified array
	int size;															//Number of elements of the array that simplifies the TTCM
	int c_simple;														//Counter for the conversion of simplified array into TTCM
	int merge[50];														//Used for determining tubes with merges and the connected tubes to them
	int ref_merge[75];													//Used for establishing reference and connected merges
	int count;															//Counter for tube numbering illustration
	int c_mass, c_m = mass[1];											//Counter for mass initialization
	int c_pdrop;														//Counter of the "longest" possible path
	double pdrop;														//Initialized pressure drop per tube
	bool feasible, check_loop;											//Argument for connectivity checking
	bool start, assign, find;											//Arguments for the Breadth Depth Search (path finding)

	SA_finned = (int)(tube_L / fin_S) * 2 * tube_adjH * tube_adjV;
	Ac_fin = (tube_L - (((int)(tube_L / fin_S))*fin_T)) * tube_adjV;
	SA_tot = SA_finned + ((tube_L - (((int)(tube_L / fin_S))* fin_T)) * PI * tube_D_o);

	clock_t start_time, end_time;
	start_time = clock();
	
	cout << "Enter the number of rows and columns, separated by a space." << endl;
	cin >> Nrow >> Ncol;

	Ntubes = Nrow * Ncol;

	do
	{
		Ninlet = 0, Noutlet = 0, Nloop = 0, Ncircuit = 1;
		count = 1;
		temp_merge = 0, temp_out = 0;
		
		feasible = true, check_loop = true;
		start = true, assign = false, find = false;

		//===Tube-Tube Connectivity Matrix (TTCM)===//
		//TTCM Initialization
		for (int i = 1; i <= Ntubes; i++)
		{
			mass[i] = 0;
			mass_done[i] = 0;
			loop[i] = 0;

			for (int j = 1; j <= Ntubes; j++)
			{
				tube_mtrx[i][j] = 0;
				checker[i][j] = 0;
			}
		}
		
		
		//===User-Defined Connectivity===//
		//Tube Numbering Illustration									//Shows the tube numbering to the user
		for (int j = 0; j < Ncol; j++)
		{
			for (int i = 0; i < Nrow; i++)
			{
				tube_mtrx[i][j] = count;
				count++;
			}
		}

		cout << "\nTube numbering is shown below" << endl;

		for (int i = 0; i < Nrow; i++)
		{
			for (int j = 0; j < Ncol; j++)
			{
				cout << "\t" << tube_mtrx[i][j];
			}

			cout << "\n";
		}

		//Input of Tube Connectivity						
		int x, y;
		
		for (int i = 1; i <= (Ntubes - 1); i++)
		{
			cout << "\nEnter the tube connection (Note: source <space> destination)" << endl;
			cin >> x >> y;

			tube_mtrx[x][y] = 1;
			tube_mtrx[y][x] = -1;
		}

		Npair = 0; 

		for (int i = 1; i <= Ntubes; i++)
		{
			for (int j = 1; j <= Ntubes; j++)
			{
				if (tube_mtrx[i][j] == 1) {Npair++;}
			}
		}
		
		/*
		//Randomized TTCM
		size = ((Ntubes*Ntubes) - Ntubes) / 2;
		c_simple = 1;

		for (int i = 0; i <= Ntubes; i++)
		{
			simple[i] = 0;
		}

		do
		{
			Npair = 0;

			for (int i = 1; i <= (Ntubes - 1); i++)					//Random assignment of -1, 0, or 1 to a random element of the simplified array
			{
				bool check;
				int n;

				srand((unsigned)time(NULL));

				do
				{
					check = true;
					n = 1 + rand() % size;

					for (int j = 1; j < i; j++)
					{
						if (random[j] == n)							//Checks if the randomized element has been previously chosen
						{
							check = false;
							break;
						}
					}
				} while (!check);

				random[i] = n;
				simple[n] = -1 + (rand() % 3);						//Assigns a value of either -1, 0, or 1
				Npair = Npair + abs(simple[n]);
			}
		} while (Npair < int(Ntubes / 2));							//Avoids too much 0s assigned

		for (int i = 1; i <= Ntubes - 1; i++)
		{
			for (int j = i + 1; j <= Ntubes; j++)
			{
				tube_mtrx[i][j] = simple[c_simple];
				tube_mtrx[j][i] = -simple[c_simple];
				c_simple++;
			}
		}
		

		for (int i = 1; i <= Ntubes; i++)
		{
			c_pos[i] = 0;											//Counts the number of (+1) in a row of the TTCM
			c_neg[i] = 0;											//Counts the number of (-1) in a row of the TTCM

			for (int j = 1; j <= Ntubes; j++)
			{
				if (tube_mtrx[i][j] == 1) {c_pos[i]++;}
				if (tube_mtrx[i][j] == -1) {c_neg[i]++;}
			}

			if (c_neg[i] == 0) {Ninlet++;}
			if (c_pos[i] == 0) {Noutlet++;}
		}

		//===Constraints checking===//
		//Non-connected tube//
		for (int i = 1; i <= Ntubes; i++)
		{
			if (c_pos[i] == 0 && c_neg[i] == 0)
			{
				feasible = false;
				cout << "\nINFEASIBLE. Single pass at tube " << i << "." << endl;
			}
		}

		//Long U-bend length//
		for (int i = 1; i <= Ntubes; i++)
		{
			for (int j = i + 1; j <= Ntubes; j++)
			{
				xcol = ((i - 1) / Nrow) + 1;						//Identifies the column number of the observed tube
				xrow = i - (Nrow*(xcol - 1));						//Identifies the row number of the observed tube

				ycol = ((j - 1) / Nrow) + 1;						//Identifies the column number of the connected tube
				yrow = j - (Nrow * (ycol - 1));						//Identifies the row number of the connected tube

				if (abs(xcol - ycol) > 2)
				{
					if (tube_mtrx[i][j] != 0)
					{
						feasible = false;
						cout << "\nCONNECTION NOT ALLOWED. Long U-bend between tube[" << i << "] and tube [" << j << "]." << endl;
					}
				}

				if (abs(xrow - yrow) > 2)
				{
					if (tube_mtrx[i][j] != 0)
					{
						feasible = false;
						cout << "\nCONNECTION NOT ALLOWED. Long U-bend between tube[" << i << "] and tube [" << j << "]." << endl;
					}
				}
			}
		}

		//Internal Looping// 
		//Breadth Depth Search
		do
		{
			c_mass = 0;

			if (start == true)										//Starts at a tube in a new (parallel) circuit
			{
				Nloop++;											//Counts the number of parallel circuits

				for (int i = 1; i <= Ntubes; i++)
				{
					bool check = true;

					for (int j = 1; j <= Ntubes; j++)
					{
						if (mass[j] == i)							//Checks if the observed tube is already in the dummy array
						{
							check = false;
							break;
						}
					}

					if (check == true)
					{
						c_m = i;
						mass_initial(c_m, Ntubes);						//Puts (unique) tubes in the dummy array
						break;
					}
				}

				for (int i = 1; i <= Nloop; i++)
				{
					if (loop[i] == 0) {loop[i] = c_m;}
				}

				start = false;
				assign = true;
			}

			if (assign == true)									//Searches for all the connected tubes to the observed tube
			{
				for (int i = 1; i <= Ntubes; i++)
				{
					if (tube_mtrx[c_m][i] == 1 || tube_mtrx[c_m][i] == -1) {mass_initial(i, Ntubes);}

					if (i == Ntubes)
					{
						for (int j = 1; j <= Ntubes; j++)
						{
							if (mass_done[j] == 0)
							{
								mass_done[j] = c_m;					//Observed tube is completely explored
								break;
							}
						}
					}
				}

				assign = false;
				find = true;
			}

			if (find == true)										//Selects the next tube to be explored in the dummy array
			{
				bool check = true;

				for (int i = 1; i <= Ntubes; i++)
				{
					for (int j = 1; j <= Ntubes; j++)
					{
						if (mass[i] == mass_done[j])				//Will not select the tube that has been completely explored
						{
							break;
						}

						if (j == Ntubes)							//Will continue on the tube that has not been completely explored
						{
							c_m = mass[i];
							find = false;
							assign = true;
							check = false;
						}
					}

					if (check == false) {break;}

					if (i == Ntubes)								//If all tubes in the dummy array were completely explored, this indicates there is another (parallel) circuit
					{
						find = false;
						start = true;
					}
				}
			}

			for (int i = 1; i <= Ntubes; i++)
			{
				if (mass_done[i] != 0) {c_mass++;}
			}

		} while (c_mass < Ntubes);
		
		if (Npair > Ntubes - Nloop)
		{
			feasible = false; check_loop = false;
			cout << "\nINFEASIBLE. There is an internal looping." << endl;
		}

		if (check_loop == true)
		{
			assign = true, find = false;

			//Re-initialization of BFS
			for (int i = 1; i <= Ntubes; i++)
			{
				mass[i] = 0;
				mass_done[i] = 0;
			}

			for (int i = 1; i <= Ntubes; i++)
			{
				if (c_neg[i] == 0)
				{
					for (int j = 1; j <= Ntubes; j++)
					{
						if (tube_mtrx[i][j] == 1)
						{
							checker[i][j] = 1;
							checker[j][i] = 1;
						}
					}

					for (int j = 1; j <= Ntubes; j++)
					{

						if (mass[j] == 0)
						{
							mass[j] = i;
							mass_initial(i, Ntubes);
							break;
						}
					}
				}
			}
			
			c_m = mass[1];

			do
			{
				c_mass = 0;

				if (assign == true)
				{
					for (int i = 1; i <= Ntubes; i++)
					{
						if (tube_mtrx[c_m][i] == 1)
						{
							mass_initial(i, Ntubes);
						
							for (int j = 1; j <= Ntubes; j++)
							{
								if (tube_mtrx[c_m][j] == -1)
								{
									if (checker[c_m][j] == 1)
									{
										checker[c_m][i] = -1;
										checker[i][c_m] = -1;
									}

									else
									{
										checker[c_m][i] = 1;
										checker[i][c_m] = 1;
									}
								}
							}
						}

						if (i == Ntubes)
						{
							for (int j = 1; j <= Ntubes; j++)
							{
								if (mass_done[j] == 0)
								{
									mass_done[j] = c_m;
									break;
								}
							}
						}
					}

					assign = false;
					find = true;
				}

				if (find == true)
				{
					bool check = true;

					for (int i = 1; i <= Ntubes; i++)
					{
						for (int j = 1; j <= Ntubes; j++)
						{
							if (mass[i] == mass_done[j])
							{
								break;
							}

							if (j == Ntubes)
							{
								if (c_neg[mass[i]] > 1)
								{
									int c_merge = c_neg[mass[i]];

									for (int k = 1; k <= Ntubes; k++)
									{
										if (tube_mtrx[mass[i]][mass_done[k]] == -1)
										{
											c_merge--;
										}
									}

									if (c_merge == 0)
									{
										c_m = mass[i];
										assign = true;
										find = false;
										check = false;
									}

									else {break;}
								}

								else
								{
									c_m = mass[i];
									assign = true;
									find = false;
									check = false;
								}
							}
						}

						if (check == false) {break;}
					}
				}

				for (int i = 1; i <= Ntubes; i++)
				{
					if (mass_done[i] != 0)
					{
						c_mass++;
					}
				}

			} while (c_mass < Ntubes);
		}

		//Merging or (final) output tubes are not on the same side
		for (int i = 1; i <= Ntubes; i++)
		{
			temp_merge = 0;

			if (c_pos[i] == 0)
			{
				for (int j = 1; j <= Ntubes; j++)
				{
					if (tube_mtrx[i][j] == -1)
					{
						temp_out += checker[i][j];
					}
				}
			}

			if (c_neg[i] > 1)
			{
				for (int j = 1; j <= Ntubes; j++)
				{
					if (tube_mtrx[i][j] == -1)
					{
						temp_merge += checker[i][j];
					}
				}

				if ((c_neg[i] - abs(temp_merge)) != 0)
				{
					feasible = false;
					cout << "\nINFEASIBLE. Merging at tube " << i << " is at the opposite ends.";
				}
			}
		}

		if (Noutlet - abs(temp_out) != 0)
		{
			feasible = false;
			cout << "\nINFEASIBLE. Outlet tubes are at the opposite ends.";
		}

		cout << "\n";

	} while (!feasible);


	//Printout of the Tube-Tube Connectivity Matrix
	cout << "\nTube-Tube Connectivity Matrix\n";
	cout << "\t";

	for (int i = 1; i <= Ntubes; i++)
	{
		cout << "\t" << i;
	}

	for (int i = 1; i <= Ntubes; i++)
	{
		cout << "\n\t" << i;

		for (int j = 1; j <= Ntubes; j++) {cout << "\t" << tube_mtrx[i][j];}
	}

	cout << "\n\n";
	 
	for (int i = 1; i <= Ntubes; i++)
	{
		if (c_neg[i] == 0) {cout << "The input tube is at tube [" << i << "]." << endl;}
		if (c_pos[i] == 0) {cout << "The output tube is at tube [" << i << "]." << endl;}
		if (c_neg[i] > 1) {cout << "Merging at tube [" << i << "]." << endl;}
		if (c_pos[i] > 1) {cout << "Splitting at tube [" << i << "]." << endl;}
	}
	
	//===Input Conditions===//
	P_in_r1 = 1500.0;													//refrigerant inlet pressure
	P_out_r1 = 1499.0;													//refrigerant outlet pressure (arbitrarily assigned value for Set 2)
	h_in_r1 = 250.0;													//refrigerant inlet specific enthalpy
	G_r1 = 0.015;														//refrigerant mass flowrate

	ref.change_pure(1);
	ref.state_ph(P_in_r1, h_in_r1);
	T_in_r1 = ref.Rc.T;													//refrigerant inlet temperature

	T_air1 = 35.0;														//air inlet temperature
	G_air1 = 0.27;														//air mass flowrate

	ref.change_pure(2);
	ref.state_tp(T_air1, P_air1);
	double rho_a1 = ref.Rc.rho;											//air inlet density

	V_air1 = G_air1 / ((tube_L - (int)(tube_L / fin_S) * fin_T) * tube_adjV *double(Nrow)* rho_a1);		//air inlet velocity

	Sd = pow(((pow(tube_adjH, 0.5)) + (pow((tube_adjV / 2), 0.5))), 0.5);

	if (Sd < ((tube_adjV + tube_D_o) / 2))
	{
		Vmax_a = (tube_adjV / (2 * (Sd - tube_D_o))) * V_air1;
	}

	else
	{
		Vmax_a = (tube_adjV / (tube_adjV - tube_D_o)) * V_air1;
	}

	//===Local Variables Initialization===//
	//Mass initialization
	for (int i = 1; i <= Ntubes; i++)
	{
		G_r[i] = 0;
	}

	for (int i = 1; i <= Ntubes; i++)
	{
		//Air-side parameters
		T_in_air[mass_done[i]] = T_air1;
		T_out_air[mass_done[i]] = T_air1;
		P_in_air[mass_done[i]] = P_air1;
		G_air[mass_done[i]] = G_air1 / double(Nrow);

		//Refrigerant-side parameters
		h_in_r[mass_done[i]] = h_in_r1;
		h_out_r[mass_done[i]] = h_in_r1;
		T_out_r[mass_done[i]] = T_in_r1;
		T_in_r[mass_done[i]] = T_in_r1;
		T_out_r[mass_done[i]] = T_in_r1;


		if (c_neg[mass_done[i]] == 0)
		{
			G_r[mass_done[i]] = G_r1 / Ninlet;
			P_in_r[mass_done[i]] = P_in_r1;

			ref.change_pure(1);
			ref.state_ph(P_in_r1, h_in_r1);
			double x_r = ref.Rc.x;

			if (x_r < 0 || x_r > 1)
			{
				ref.change_pure(1);
				ref.state_tp(T_in_r1, P_in_r1);
				myu_r1 = ref.Rc.visc / 1000000;
				rho_r1 = ref.Rc.rho;
			}

			else
			{
				ref.change_pure(1);
				ref.sat_p(P_in_r1);
				double myu_rl = ref.Rl.visc / 1000000;					//in Pa-s
				double myu_rv = ref.Rv.visc / 1000000;
				double rho_rl = ref.Rl.rho;								//in kg/m^3
				double rho_rv = ref.Rv.rho;

				myu_r1 = (x_r * myu_rv) + ((1 - x_r) * myu_rl);
				double VF = 1 / (1 + (((1 - x_r) / x_r) * (rho_rv / rho_rl)));
				rho_r1 = (VF * rho_rv) + ((1 - VF) * rho_rl);
			}

			DP = 128 * tube_L * G_r[mass_done[i]] * myu_r1 / (PI * rho_r1 * pow(tube_D_i, 4.0) * 1000);		//in kPa

			P_out_r[mass_done[i]] = P_in_r[mass_done[i]] - DP;
		}

		else
		{
			for (int j = 1; j <= Ntubes; j++)
			{
				if (tube_mtrx[mass_done[i]][j] == -1)
				{
					if (c_neg[mass_done[i]] > 1)
					{
						G_r[mass_done[i]] = G_r[mass_done[i]] + G_r[j];

						double temp_P_in = P_out_r[j];

						for (int k = 1; k <= Ntubes; k++)								//Search for the lowest outlet pressure among the merging flows
						{
							if (tube_mtrx[mass_done[i]][k] == -1 && k != j)
							{
								if (P_out_r[k] < temp_P_in && P_out_r[k] > 0)
								{
									temp_P_in = P_out_r[k];
								}
							}
						}

						for (int k = 1; k <= Ntubes; k++)
						{
							if (tube_mtrx[mass_done[i]][k] == -1) {P_out_r[k] = temp_P_in;}
						}

						P_in_r[mass_done[i]] = temp_P_in;
						P_out_r[mass_done[i]] = P_in_r[mass_done[i]] - DP;
					}
					
					else
					{
						if (c_pos[j] > 1)
						{
							G_r[mass_done[i]] = G_r[j] / c_pos[j];
							P_in_r[mass_done[i]] = P_out_r[j];
							P_out_r[mass_done[i]] = P_in_r[mass_done[i]] - DP;
						}

						else
						{
							G_r[mass_done[i]] = G_r[j];
							P_in_r[mass_done[i]] = P_out_r[j];
							P_out_r[mass_done[i]] = P_in_r[mass_done[i]] - DP;
						}
					}		
				}
			}				
		}
	}

	if (Noutlet > 1)
	{
		for (int i = 1; i <= Ntubes; i++)
		{
			if (c_pos[i] == 0)
			{
				double temp_P_out = P_out_r[i];

				for (int j = 1; j <= Ntubes; j++)
				{
					if (c_pos[j] == 0 && j != i)
					{
						if (P_out_r[j] < temp_P_out) {temp_P_out = P_out_r[j];}
					}
				}

				for (int j = 1; j <= Ntubes; j++)
				{
					if (c_pos[j] == 0) {P_out_r[j] = temp_P_out;}
				}
			}
		}
	}

	for (int i = 1; i <= Ntubes; i++)
	{
		ref.change_pure(1);
		ref.state_ph(P_in_r1, h_in_r1);
		x_in_r[i] = ref.Rc.x;
		x_out_r[i] = ref.Rc.x;
	}

	//Printout of initialized variables
	for (int i = 1; i <= Ntubes; i++)
	{
		cout << "\n\nAt tube_" << i << "\trefrigerant's\tP_in = " << P_in_r[i] << "\t\tP_out = " << P_out_r[i] << "\t\th_in = " << h_in_r[i] << "\t\tx_in = " << x_in_r[i] << "\t\tT_out = " << T_out_r[i] << "\t\tG_r = " << G_r[i];
		cout << "\n\t\tair's\t\tP_in = " << P_in_air[i] << "\t\tT_out = " << T_out_air[i] << "\t\tT_in = " << T_in_air[i];
	}
	cout << endl << "\n\n";

	//===Newton-Rahpson Calculation===//
	CNewtonRaphsonMethod mnm;
	mnm.setup(2000, 1 + 1e-12, 1e-4);

	p = 0;

	for (int i = 1; i <= Ntubes; i++)
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

	mnm.setAcc(1);
	mnm.initial();

	for (mnm.main_loop_init(); mnm.main_loop_check(); mnm.main_loop_reinit())
	{
		for (mnm.sub_loop_init(); mnm.sub_loop_check(); mnm.sub_loop_reinit())
		{
			p = 0;

			for (int i = 1; i <= Ntubes; i++)
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
			for (int a = 1; a <= Ntubes; a++)
			{
				int i = mass_done[a];

				//Fin efficiency
				eta_fin[i] = 0.77;

				//Fouling resistance factor
				Rf = 0.0;

				//Fin-contact resistance factor
				Rc = 0.0;

				//Refrigerant- and air-heat transfer coefficients
				ref.change_pure(1);
				ref.sat_p(P_in_r[i]);
				double myu_rl = ref.Rl.visc / 1000000;					//in Pa-s
				double myu_rv = ref.Rv.visc / 1000000;
				double cp_rl = ref.Rl.cp * 1000;						//in J/kg-K
				double cp_rv = ref.Rv.cp * 1000;
				double k_rl = ref.Rl.thc;								//in W/m-k
				double k_rv = ref.Rv.thc;
				double rho_rl = ref.Rl.rho;								//in kg/m^3
				double rho_rv = ref.Rv.rho;
				double h_rl = ref.Rl.h;									//in kJ/kg
				double h_rv = ref.Rv.h;
				double st = ref.Rl.st;									//in N/m

				ref.change_pure(1);
				ref.state_ph(P_in_r[i], h_in_r[i]);
				x_in_r[i] = ref.Rc.x;

				double rho_r, k_r, myu_r;

				if (x_in_r[i] < 0 || x_in_r[i] > 1)
				{
					myu_r = ref.Rc.visc / 1000000;
					k_r = ref.Rc.thc;
					double cp_r = ref.Rc.cp * 1000;
					rho_r = ref.Rc.rho;

					Re_r[i] = G_r[i] * 4.0 / PI / tube_D_i / myu_r;
					Prandtl_r[i] = myu_r * cp_r / k_r;
				}

				else
				{
					myu_r = (x_in_r[i] * myu_rv) + ((1 - x_in_r[i]) * myu_rl);
					double VF = 1 / (1 + (((1 - x_in_r[i]) / x_in_r[i]) * (rho_rv / rho_rl)));
					rho_r = (VF * rho_rv) + ((1 - VF) * rho_rl);
				}


				ref.change_pure(2);
				ref.state_tp(T_in_air[i], P_in_air[i]);

				double rho_a = ref.Rc.rho;
				double myu_a = ref.Rc.visc / 1000000;
				double k_a = ref.Rc.thc;
				double cp_a = ref.Rc.cp * 1000;

				Re_a[i] = rho_a * V_air1 * tube_D_o / myu_a;
				Re_a_max[i] = rho_a * Vmax_a * tube_D_o / myu_a;
				Prandtl_a[i] = myu_a * cp_a / k_a;


				//Water-side heat transfer coefficient from Dittus-Boelter

				if (x_in_r[i] < 0 || x_in_r[i] > 1)
				{
					alfa_r[i] = 0.023 * k_r * pow(Re_r[i], 0.8) * pow(Prandtl_r[i], 0.3) / (1000 * tube_D_i);
				}

				else
				{
					double F;
					double Xtt = (pow(((1 - x_in_r[i]) / x_in_r[i]), 0.9)) * (pow((rho_rv / rho_rl), 0.5)) * (pow((myu_rl / myu_rv), 0.1));

					if ((1 / Xtt) < 0.1) { F = 1; }

					else
					{
						F = 2.35 * (pow((0.213 + (1 / Xtt)), 0.736));
					}

					double Re_l = (4 * (1 - x_in_r[i]) * G_r[i]) / (PI * tube_D_i * myu_rl);
					double Pr_l = myu_rl * cp_rl / k_rl;
					double hl = 0.023 * pow(Re_l, 0.8) * pow(Pr_l, 0.4) * k_rl / tube_D_i;
					double Re_tp = Re_l * pow(F, 1.25);
					double S = 1 / (1 + (0.00000253 * pow(Re_tp, 1.17)));
					double hfz = (0.00122 * pow(k_rl, 0.79) * pow(cp_rl, 0.45) * pow(rho_rl, 0.49) * pow((T_in_air[i] - T_in_r[i]), 0.24) * pow((0.01 * (P_in_r[i] - P_out_r[i])), 0.75)) / (pow(st, 0.5) * pow(myu_rl, 0.29) * pow(((h_rv - h_rl) * 1000), 0.24) * pow(rho_rl, 0.24));

					alfa_r[i] = (S * hfz) + (F * hl);
				}


				//Air-side heat transfer coefficient from Kim-Youn-Webb
				jcf = (rho_a*Vmax_a*cp_a) * pow(Prandtl_a[i], 1.5) * 0.163 * pow(Re_a_max[i], -0.369) * pow((tube_adjV / tube_adjH), 0.106) * pow((fin_S / tube_D_o), 0.0138) * pow((tube_adjV / tube_D_o), 0.13);

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

				//Equations for iteration							//Observes the downstream (destination) tube(s)
				G_res = 0.0;

				if (c_pos[i] == 0)									//Outlet tube
				{
					if (Noutlet > 1)								//Multiple outlet tubes
					{
						if (i == ref_out)							//Reference outlet tube
						{
							G_res = G_r1;

							for (int j = 1; j <= Ntubes; j++)
							{
								if (c_pos[j] == 0 && j != i)
								{
									G_res = G_res - G_r[j];
								}
							}

							mnm.setError(p, G_r[i], G_res);
							p++;

							f[i] = 64.0 * Re_r[i];

							//Hagen-Poiseulle equation
							mnm.setError(p, -128.0 * tube_L * G_r[i] * myu_r / (PI * rho_r * pow(tube_D_i, 4.0)), (P_out_r[i] - P_in_r[i]) * 1000.0);
							p++;
						}

						else										//Non-reference outlet tube(s)
						{
							f[i] = 64.0 * Re_r[i];

							//Hagen-Poiseulle equation
							mnm.setError(p, G_r[i], (P_out_r[i] - P_in_r[i]) * 1000.0 * PI * rho_r * pow(tube_D_i, 4.0) / (-128.0 * tube_L * myu_r));
							p++;

							mnm.setError(p, P_out_r[i], P_out_r[ref_out]);		//Should have the same outlet pressure with the reference outlet tube
							p++;
						}
					}

					else											//Single outlet tube
					{
						mnm.setError(p, G_r[i], G_r1);
						p++;

						f[i] = 64.0 * Re_r[i];

						//Hagen-Poiseulle equation
						mnm.setError(p, -128.0 * tube_L * G_r[i] * myu_r / (PI * rho_r * pow(tube_D_i, 4.0)), (P_out_r[i] - P_in_r[i]) * 1000.0);
						p++;
					}
				}

				else
				{
					for (int j = 1; j <= Ntubes; j++)					//Observes the downstream (destination) flow(s)
					{
						if (tube_mtrx[i][j] == 1)
						{
							G_res = G_res + G_r[j];
						}
					}

					mnm.setError(p, G_r[i], G_res);
					p++;

					f[i] = 64.0 * Re_r[i];

					//Hagen-Poiseulle equation
					mnm.setError(p, -128.0 * tube_L * G_r[i] * myu_r / (PI * rho_r * pow(tube_D_i, 4.0)), (P_out_r[i] - P_in_r[i]) * 1000.0);
					p++;
				}

				h_out_r[i] = h_in_r[i] + Q[i] / G_r[i];

				ref.change_pure(1);
				ref.state_ph(P_out_r[i], h_out_r[i]);
				double T_out_rx = ref.Rc.T;

				mnm.setError(p, T_out_r[i], T_out_rx);
				p++;

				ref.change_pure(2);
				ref.state_tp(T_in_air[i], P_in_air[i]);
				double cp_ax = ref.Rc.cp;

				DT = Q[i] / G_air[i] / cp_ax;

				mnm.setError(p, T_out_air[i], T_in_air[i] - DT);
				p++;

			}

			mnm.prt_sum();
		}

		//Input to the next tube/column								//Observes the downstream (destination) flow(s)
		for (int a = 1; a <= Ntubes; a++)
		{
			int i = mass_done[a];
			//int i = trial[a];
			//Air-side
			if (i <= Ntubes - Nrow)
			{
				T_in_air[i + Nrow] = T_out_air[i];
				P_in_air[i + Nrow] = P_out_air[i];
			}
			else {}

			//Refrigerant-side
			//Mass and energy (enthalpy) balance at the junction
			G_res = 0.0;
			H_res = 0.0;

			if (c_pos[i] == 0)
			{
				for (int j = 1; j <= Ntubes; j++)
				{
					if (tube_mtrx[i][j] == -1)
					{
						P_out_r[j] = P_in_r[i];
					}
				}
			}

			else
			{
				for (int j = 1; j <= Ntubes; j++)
				{
					G_res = G_r[i];
					H_res = G_res * h_out_r[i];

					if (tube_mtrx[i][j] == 1)							//Downstream flow
					{
						if (c_neg[j] > 1)								//Connected tube has a merge
						{
							for (int k = 1; k <= Ntubes; k++)
							{
								if (tube_mtrx[j][k] == -1)
								{
									if (k == i) {}

									else
									{
										G_res = G_res - G_r[i];
										H_res = H_res - (G_r[i] * h_out_r[i]);
									}
								}
							}
						}
					}

					//Assigning to connected tube
					for (int j = 1; j <= Ntubes; j++)
					{
						if (tube_mtrx[i][j] == 1)
						{
							h_in_r[j] = H_res / G_res;
							P_in_r[j] = P_out_r[i];

							ref.change_pure(1);
							ref.state_ph(P_in_r[j], h_in_r[j]);
							double T_rt = ref.Rc.T;

							T_in_r[j] = T_rt;
						}
					}
				}
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
	for (int i = 1; i <= Ntubes; i++)
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

		cout << "\n\nAt tube_" << i << "\trefrigerant's\tP_in = " << P_in_r[i] << "\t\tP_out = " << P_out_r[i] << "\t\th_out = " << h_out_r[i] << "\t\th_in = " << h_in_r[i] << "\t\tx_out = " << x_out_r[i] << "\t\tx_in = " << x_in_r[i] << "\t\tT_out = " << T_out_r[i] << "\t\tG_r = " << G_r[i];
		cout << "\n\t\tair's\t\tP_in = " << P_in_air[i] << "\t\tT_out = " << T_out_air[i] << "\t\tT_in = " << T_in_air[i];
	}
	cout << endl;

	Atot = (tube_D_o * PI * (tube_L - (int)(tube_L / fin_S) * fin_T) + (int)(tube_L / fin_S) * 2.0 * (tube_adjV * tube_adjH - PI / 4.0 * tube_D_o * tube_D_o))* double(Ntubes);
	cout << "\nQtot = " << Qtot << " kW\tDP = " << DP_i + DP_o << " kPa\tGtot = " << Gtot << " kg/s\tAtot = " << Atot << " m^2";

	end_time = clock();
	double run_time = (double)(end_time - start_time) / 1000;			//[sec]
	cout << "\n\nTotal run time: " << run_time << " s";

	//Writing results on .csv file
	cout << "\n\nWriting to .csv file...";

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

	///1
	int i = 0;
	for (int c = 0; c < Ncol; c++)
	{
		for (int r = 0; r < Nrow; r++)
		{
			++i;
			T_in_r_mtrx[r][c] = T_in_r[i];
		}
	}

	complex << "\nRefrigerant Inlet Temperature\n";
	for (int r = 0; r < Nrow; r++)
	{
		for (int c = 0; c < Ncol; c++)
		{
			complex << T_in_r_mtrx[r][c] << ",";
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
			T_out_r_mtrx[r][c] = T_out_r[i];
		}
	}

	complex << "\nRefrigerant Outlet Temperature\n";
	for (int r = 0; r < Nrow; r++)
	{
		for (int c = 0; c < Ncol; c++)
		{
			complex << T_out_r_mtrx[r][c] << ",";
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
			P_in_r_mtrx[r][c] = P_in_r[i];
		}
	}

	complex << "\nRefrigerant Inlet Pressure\n";
	for (int r = 0; r < Nrow; r++)
	{
		for (int c = 0; c < Ncol; c++)
		{
			complex << P_in_r_mtrx[r][c] << ",";
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
			P_out_r_mtrx[r][c] = P_out_r[i];
		}
	}

	complex << "\nRefrigerant Outlet Pressure\n";
	for (int r = 0; r < Nrow; r++)
	{
		for (int c = 0; c < Ncol; c++)
		{
			complex << P_out_r_mtrx[r][c] << ",";
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
			T_in_air_mtrx[r][c] = T_in_air[i];
		}
	}

	complex << "\nAir Inlet Temperature\n";
	for (int r = 0; r < Nrow; r++)
	{
		for (int c = 0; c < Ncol; c++)
		{
			complex << T_in_air_mtrx[r][c] << ",";
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
			T_out_air_mtrx[r][c] = T_out_air[i];
		}
	}

	complex << "\nAir Outlet Temperature\n";
	for (int r = 0; r < Nrow; r++)
	{
		for (int c = 0; c < Ncol; c++)
		{
			complex << T_out_air_mtrx[r][c] << ",";
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
			h_in_r_mtrx[r][c] = h_in_r[i];
		}
	}

	complex << "\nRefrigerant Inlet Specific Enthalpy\n";
	for (int r = 0; r < Nrow; r++)
	{
		for (int c = 0; c < Ncol; c++)
		{
			complex << h_in_r_mtrx[r][c] << ",";
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
			h_out_r_mtrx[r][c] = h_out_r[i];
		}
	}

	complex << "\nRefrigerant Outlet Specific Enthalpy\n";
	for (int r = 0; r < Nrow; r++)
	{
		for (int c = 0; c < Ncol; c++)
		{
			complex << h_out_r_mtrx[r][c] << ",";
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
			G_r_mtrx[r][c] = G_r[i];
		}
	}

	complex << "\nRefrigerant Flowrate\n";
	for (int r = 0; r < Nrow; r++)
	{
		for (int c = 0; c < Ncol; c++)
		{
			complex << G_r_mtrx[r][c] << ",";
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
			complex << tube_mtrx[i][j] << ",";
		}
	}

	cout << "\n\n";
	system("pause");
	return 0;
}

void mass_initial(int K, int Ntubes)								//Allocates a value into the path array
{
	bool check = false;

	for (int i = 1; i <= Ntubes; i++)
	{
		if (mass[i] == 0)											//Puts K in an "empty" slot of the path array
		{
			for (int j = 1; j <= i; j++)							//Checks if K has already been in the path array
			{
				if (mass[j] == K)
				{
					check = true;
					break;
				}

				else {}
			}

			if (check == false)
			{
				mass[i] = K;
			}

			else {}

			break;
		}

		else {}
	}
}*/
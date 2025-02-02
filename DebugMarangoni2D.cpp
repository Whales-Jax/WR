/*#include<math.h>
#include <iostream>
#include<random>
#include <string>
#include <fstream>
using namespace std;
#include "CNewtonRaphsonMethodPlus.h"
#include"CFluidParameter.h"
#include"PropertyLithiumBromide_ver05.h"
#include"PropertyWater_ver05.h"

PropertyLithiumBromide libr;
PropertyWater wat;

double Tmax(double a, double b, double c) {
	double max;
	max = a;
	if (a <= b) {
		max = b;
	}
	if (max <= c) {
		max = c;
	}
	return max;
}

int main()
{
	int N;              /////Angle's mesh
	int M;              /////Thickness's mesh
	int E;              /////Time-span
	int i;
	int j;
	int t;
	int p;
	int c = 0;
	int e = 0;
	double p_v;
	double T_sat;
	double X_in;
	double T_e;
	double T_in;
	double T_w;
	double X_e;
	double m;
	double T_ave;
	double X_ave;
	double W;
	double PI;
	double beta;
	double rho;
	double myu;
	double g;
	double k;
	double d;
	double cp;
	double a;
	double Pr;
	double Sc;
	double Re;
	double h;
	double L;
	double st;
	double dT;
	double dw;
	double ds_T;
	double ds_w;
	double Ma_T;
	double Ma_w;
	double time;
	double dt;
	double dX;
	double dY;
	double m_w;
	double mm_w;
	double hv;
	double Dx;
	double Dy;
	double uu[60][35];
	double vv[60][35];
	double tt[60][35];
	double xx[60][35];
	double PP[60][35];
	double u[60][35];
	double v[60][35];
	double T[60][35];
	double X[60][35];
	double P[60][35];
	double aP0;
	double aW_U;
	double aE_U;
	double aN_U;
	double aS_U;
	double aP_U;
	double aW_V;
	double aE_V;
	double aN_V;
	double aS_V;
	double aP_V;
	double aW_T;
	double aE_T;
	double aN_T;
	double aS_T;
	double aP_T;
	double aW_X;
	double aE_X;
	double aN_X;
	double aS_X;
	double aP_X;


	N = 40;              /////Angle's mesh
	M = 20;              /////Thickness's mesh
	E = 12000;              //////Time-span

	ofstream fout("2DMarangoniConvection20210302wa55%40c_ho10_X,rho,MaT改.csv");
	if (!fout) {
		cout << "ファイルをオープンできませんでした" << endl;
		return 1;
	}
	else
		cout << "ファイルをオープンしました" << endl;

	////////////////////Experiment condition////////////////////

	//p_v = 5.52;		                             //////[kPa]
	//T_sat = wat.p_t(p_v);                        //////[oC]
	T_in = 30.0;                        //////[oC]
	X_in = 0.60;		                         //////[-]
	T_sat = libr.sc_Tsat_XT(X_in, T_in) + 30.0;
	p_v = wat.t_p(T_sat);
	T_e = libr.sc_T_XTsat(X_in, T_sat);   //////[oC]
	//T_e = 40;
	//X_in = 0.55;
	//T_sat = libr.sc_Tsat_XT(X_in, T_e)+10.0;
	//p_v = wat.t_p(T_sat);
	//T_in = T_e - 20.0;                            //////[oC]
	T_w = T_in;								     //////[oC]
	X_e = libr.sc_X_TTsat(T_in, T_sat);          //////[-]
	m = 0.01808;		                         //////[kg/s]
	T_ave = T_in;
	X_ave = X_in;

	//cout << "80c" << libr.sc_rho_XT(0.55, 80) << endl;

	fout << "N" << "," << N << endl;
	fout << "M" << "," << M << endl;
	fout << "p_v" << "," << p_v << endl;
	fout << "T_sat" << "," << T_sat << endl;
	fout << "X_in" << "," << X_in << endl;
	fout << "T_in" << "," << T_in << endl;
	fout << "T_w" << "," << T_w << endl;
	fout << "T_e" << "," << T_e << endl;
	fout << "X_e" << "," << X_e << endl;
	fout << "m" << "," << m << endl;


	////////////////////Geometrical condition////////////////////

	W = 1.0;				//////Surface width [m]
	PI = (6 * asin(0.5));	//////[-]
	beta = PI / 10.0;		//////Surface inclination [rad]


	////////////////////Properties calculation////////////////////	

	rho = libr.sc_rho_XT(X_ave, T_ave);//////////////////////////////////////[kg/m3]
	myu = libr.sc_visc_XT(X_ave, T_ave);/////////////////////////////////////[Pas]	
	g = 9.81;/////////////////////////////////////////////////////////////[m/s2]
	k = libr.sc_thc_XT(X_ave, T_ave);////////////////////////////////////////[W/mK]
	d = libr.sc_d_XT(X_ave, T_ave);//////////////////////////////////////////[m2/s]
	cp = libr.sc_cp_XT(X_ave, T_ave) * 1000.0;/////////////////////////////////[J/kgK]
	a = k / rho / cp;
	Pr = myu * cp / k ;
	Sc = myu / rho / d ;
	Re = 4.0 * m / W / myu;
	h = 0.001;///////////////////////////////////////////////////////////////film thickness [m]
	L = 2.0 * h;//0.05;////////////////////////////////////////////////////////////Surface stream-wise length [m]
	st = libr.sc_st_XT(X_ave, T_ave);///////////////////////////////////////////[N/m]

	dX = 1.0 / double(N);
	dY = (h / L / double(M));
	Dx = dY / dX;
	Dy = dX / dY;

	//cout << "myu" << "," << myu << endl;
	fout << "Pr" << "," << Pr << endl;
	fout << "Sc" << "," << Sc << endl;
	fout << "Re" << "," << Re << endl;
	fout << "h" << "," << h << endl;
	fout << "L" << "," << L << endl;
	fout << "st" << "," << st << endl;

	/////////////////////////////Representative infinitesimal field disturbances
	dT = 0.05;
	dw = 0.0005;


	////////////////////// Condition with additive (Hozawa 1991)/////////////////


	if (30.0 <= T_in && T_in <= 40) {
		ds_w = -(-0.172250 * T_in + 9.88250) / 5.0 * 0.001 * 100;
	}
	else if (40.0 < T_in && T_in <= 60.0) {
		ds_w = -(-0.0387083 * T_in + 4.54083) / 5.0 * 0.001 * 100;
	}
	else if (60.0 < T_in) {
		ds_w = -(-0.0195167 * T_in + 3.38933) / 5.0 * 0.001 * 100;
	}
	else {
		ds_w = 0;
	}


	ds_T = -Tmax(-(0.0225 * T_in - 1.6667), 0.07575, 0) * 0.001;

	

	////////////////////// Condition without additive////////////////////////////

	ds_T = (libr.sc_st_XT(X_ave, T_ave + dT) - libr.sc_st_XT(X_ave, T_ave)) / (dT);
	ds_w = (libr.sc_st_XT(X_ave + dw, T_ave) - libr.sc_st_XT(X_ave, T_ave)) / (dw);

	/////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////// Marangoni numbers///////////////////////////////
	Ma_T = ds_T * L * rho * T_e / myu / myu;
				//Ma_w = ds_w * L * rho * (X_in - X_e) / myu / myu;
	Ma_w = ds_w * L * rho * (X_in) / myu / myu;
	
	//////////////////////////// Hozawa 1991////////////////////////////////////
	//Ma_T =6472;
		//Ma_w =3366;

	fout << "ds_T" << "," << ds_T << endl;
	fout << "ds_w" << "," << ds_w << endl;
	fout << "Ma_T" << "," << Ma_T << endl;
	fout << "Ma_w" << "," << Ma_w << endl;


	//////////////////////////////////time-step//////////////////////////////////
	///////////real time [s]
	time = 0;
	///////////dimensionless time step 
	dt = 0.0004; 
	///////////discretized time step coefficient
	aP0 = dX * dY / (dt * myu / (rho * L * L));
	///////////discretized dimensionless time step coefficient
	//aP0 = dX * dY / (dt);

	fout << "time[s]" << "," << dt << endl;


	//////////////////////////////////Initializing condition//////////////////////////////////

	//random_device seed_gen;
	//default_random_engine engine(seed_gen());

	//// 平均0.0、標準偏差1.0で分布させる
	//normal_distribution<> dist(0.0, 0.00005);
	//normal_distribution<> dist2(1.0, 0.01);
	//normal_distribution<> dist3(100, 10);

	/////////////////////////////////////////////////////////////////
	ifstream stream("random2.csv");
	string line;
	const string delim = ",";

	int row = 0;
	int col;
	while (getline(stream, line)) {
		col = 0;
		// delimを区切り文字として切り分け、doubleに変換してxx[][]に格納する
		for (string::size_type spos, epos = 0;
			(spos = line.find_first_not_of(delim, epos)) != string::npos;) {
			string token = line.substr(spos, (epos = line.find_first_of(delim, spos)) - spos);
			xx[col++][row] = stod(token);
		}
		++row;
	}
	/////////////////////////////////////////////////////////////////////
	/////////////////////intialization of previous time step fields///////////////////////////////////////////////
	for (j = 0; j <= M + 1; j++) {
		for (i = 0; i <= N + 1; i++) {
			uu[i][j] = 0.0;          
			vv[i][j] = 0.0;          
			//tt[i][j] = (T_in - T_w) / (T_e - T_w);  
			tt[i][j] = T_in / T_e;  
			xx[i][j] = 1.0;	
			PP[i][j] = 0.0;              
		}
	}

	for (j = 0; j <= M + 1; j++) {
		for (i = 0; i <= N + 1; i++) {
			fout << vv[i][j] << ",";
		}
		fout << 0 << endl;
	}

	for (j = 0; j <= M + 1; j++) {
		for (i = 0; i <= N + 1; i++) {
			fout << PP[i][j] << ",";
		}
		fout << 0 << endl;
	}

	for (j = 0; j <= M + 1; j++) {
		for (i = 0; i <= N + 1; i++) {
			fout << tt[i][j] << ",";
		}
		fout << 0 << endl;
	}

	for (j = 0; j <= M + 1; j++) {
		for (i = 0; i <= N + 1; i++) {
			fout << xx[i][j] << ",";
		}
		fout << 0 << endl;
	}


	////////////////////////////Iterated variables initialization////////////////////////////////


	for (j = 0; j <= M + 1; j++) {
		for (i = 0; i <= N + 1; i++) {
			u[i][j] = 0.0000001; 
			v[i][j] = 0.0000001; 
			//T[i][j] = (T_in - T_w) / (T_e - T_w);   
			T[i][j] = T_in / T_e;   
			X[i][j] = 0.9999999; 
			P[i][j] = 0.0;          
		}
	}

	hv = wat.sat_hv(T_sat) * 1000.0;                

	/////////////////////////////Time Loop/////////////////////////////////////////////////

	for (t = 1; t <= E; t++) {

		//time = time + dt / rho / L / L * myu;
		time = time + dt;
		if (t % 1 == 0) {
			fout << time << ",";
		}
		CNewtonRaphsonMethodPlus nrm;		// CNewtonRaphsonMethodクラスの宣言

		nrm.Initialize();

		p = 0;

		for (j = 0; j <= M + 1; j++) {
			for (i = 0; i <= N; i++) {
				nrm.SetValiable(p, u[i][j]);
				p++;
			}
		}



		for (j = 0; j <= M + 1; j++) {
			for (i = 0; i <= N + 1; i++) {
				nrm.SetValiable(p, P[i][j]);             //////check
				p++;
				nrm.SetValiable(p, T[i][j]);
				p++;
				nrm.SetValiable(p, X[i][j]);
				p++;
			}
		}

		for (j = 0; j <= M; j++) {
			for (i = 0; i <= N + 1; i++) {
				nrm.SetValiable(p, v[i][j]);
				p++;
			}
		}

		nrm.SetAcc(1.0);				//加速度勾配の入力　※なくても可
		nrm.SetAccAuto(true, 0.01, 100, 1.2);

		nrm.SetPrint(true, false);       //1ループごとのエラー値のディスプレイ設定
										//nrm.SetPrint( 全体エラー表示 , 変数ごとのエラー表示 );

		nrm.SetDelta(1.000001);           //微小移動量の設定(def1.0001)

		nrm.SetError(0.000001);       	//収束判定の設定(def0.0001)

		nrm.SetEndReport(true);      	//ニュートン法の終了時に最終エラーをディスプレイに表示する

		for (nrm.MainLoopInit(); nrm.MainLoopCheck(); nrm.MainLoopReinit()) {		// おまじない
			for (nrm.SubLoopInit(); nrm.SubLoopCheck(); nrm.SubLoopReinit()) {	    // おまじない

				p = 0;

				for (j = 0; j <= M + 1; j++) {
					for (i = 0; i <= N; i++) {
						nrm.GetValiable(p, u[i][j]);
						p++;
					}
				}



				for (j = 0; j <= M + 1; j++) {
					for (i = 0; i <= N + 1; i++) {
						nrm.GetValiable(p, P[i][j]);             //////check
						p++;
						nrm.GetValiable(p, T[i][j]);
						p++;
						nrm.GetValiable(p, X[i][j]);
						p++;
					}
				}

				for (j = 0; j <= M; j++) {
					for (i = 0; i <= N + 1; i++) {
						nrm.GetValiable(p, v[i][j]);
						p++;
					}
				}



				p = 0;
				///////////////////////Main////////////////////////////




				//rewrite equations

				for (j = 0; j <= M + 1; j++) {
					for (i = 0; i <= N; i++) {

						if (i == 0) {
							nrm.SetResult(p, u[i][j], 0.0);
							p++;
						}
						else if (i == N) {
							nrm.SetResult(p, u[i][j], 0.0);
							p++;
						}
						else if (j == 0) {
							nrm.SetResult(p, u[i][j], 0.0);
							p++;
						}
						else if (j == M + 1) {
							nrm.SetResult(p, u[i][j] * dX - ((u[i][j - 1] ) / 1.0) * dX, Ma_T * (T[i + 1][j] * dY - T[i][j] * dY) + Ma_w * (X[i + 1][j] * dY - X[i][j] * dY));
							//nrm.SetResult(p, u[i][j] * dX - ((u[i][j - 1]) / 1) * 2.0 * dX, Ma_T * (T[i + 1][j] * dY - T[i - 1][j] * dY) + Ma_w * (X[i + 1][j] * dY - X[i - 1][j] * dY));
							p++;
						}
						else {
							if (i == 1) {
								nrm.SetResult(p, (u[i][j] - u[i - 1][j])* dY + (v[i][j] - v[i][j - 1]) * dX, 0.0);
								p++;
							}
							else if (i == N-1) {
								nrm.SetResult(p, (u[i][j] - u[i -1][j]) * dY + (v[i][j] - v[i][j - 1]) * dX, 0.0);
								p++;
							}
							else if (j == 1) {
								nrm.SetResult(p, (u[i][j] - u[i - 1][j]) * dY + (v[i][j] - v[i][j - 1]) * dX, 0.0);
								p++;
							}
							else if (j == M ) {
								nrm.SetResult(p, (u[i][j] - u[i - 1][j])* dY + (v[i][j] - v[i][j-1]) * dX, 0.0);
								p++;
							}
							else {
								nrm.SetResult(p, (u[i][j] - u[i - 1][j]) * dY + (v[i][j] - v[i][j - 1]) * dX, 0.0);
								p++;
							}
						}

						
					}
				}



				for (j = 0; j <= M + 1; j++) {
					for (i = 0; i <= N + 1; i++) {
						

						if (i == 0) {
							nrm.SetResult(p, P[i][j], P[i + 1][j]);
							p++;
							nrm.SetResult(p, T[i][j], T[i + 1][j]);
							//nrm.SetResult(p, T[i][j], T_w / T_e);
							p++;
							nrm.SetResult(p, X[i][j], X[i + 1][j]);
							p++;
						}
						else if (i == N+1) {
							nrm.SetResult(p, P[i][j], P[i - 1][j]);
							p++;
							nrm.SetResult(p, T[i][j], T[i - 1][j]);
							//nrm.SetResult(p, T[i][j], T_w / T_e);
							p++;
							nrm.SetResult(p, X[i][j], X[i - 1][j]);
							p++;
						}
						else if (j == 0) {
							nrm.SetResult(p, P[i][j], P[i][j + 1]);
							p++;
							nrm.SetResult(p, T[i][j], T_w / T_e);
							p++;
							nrm.SetResult(p, X[i][j], X[i][j + 1]);
							p++;
						}
						else if (j == M + 1) {
							if (i == N) {
								nrm.SetResult(p, P[i][j], P[i][j-1]);
							p++;
								
								nrm.SetResult(p, T[i][j], libr.sc_T_XTsat((X[i][j] * (X_in)), wat.p_t(p_v)) / T_e);
								p++;
								nrm.SetResult(p, k* (T[i][j] - ((T[i][j - 1]) / 1.0)), -(hv * rho * d * X_in / (T_e)) * (X[i][j] - ((X[i][j - 1]) / 1.0)));////check
								p++;

							}
							else {
								nrm.SetResult(p, P[i][j], P[i][j - 1]);
							p++;
								//nrm.SetResult(p, P[i + 1][j] * dX, P[i][j] * dX + (Ma_w * (X[i + 1][j] - 2.0 * X[i][j] + X[i - 1][j]) + Ma_T * (T[i + 1][j] - 2.0 * T[i][j] + 1.0 * T[i - 1][j])));
								//p++;
								nrm.SetResult(p, T[i][j], libr.sc_T_XTsat((X[i][j] * (X_in)), wat.p_t(p_v)) / T_e);
								p++;
								//nrm.SetResult(p, k* (T[i][j] * dX - ((T[i][j - 1] + T[i][j - 2]) / 2.0) * dX), -(hv * rho * d * X_in / (T_e - T_w) * (X[i][j] * dX - ((X[i][j - 1] + X[i][j - 2]) / 2.0) * dX)));////check
								//nrm.SetResult(p, k* (T[i][j] * dX - ((T[i][j - 1]) / 1.0) * dX), -(hv * rho * d / (T_e * X[i][j]) )* (X[i][j] * dX - ((X[i][j - 1]) / 1.0) * dX));////check

								//////////////Babadi
								//nrm.SetResult(p, k* (T[i][j] - ((T[i][j - 1]) / 1.0)), -(hv * rho * d / (T_e * X[i][j])) * (X[i][j] - ((X[i][j - 1]) / 1.0)));////check

								///////////////Hozawa
								nrm.SetResult(p, k * (T[i][j] - ((T[i][j - 1]) / 1.0)), -(hv * rho * d * X_in / (T_e)) * (X[i][j] - ((X[i][j - 1]) / 1.0)));////check
								p++;
							}
							
							
						}
						else {
							
							if (i == 1) {
								aW_U = Tmax((u[i - 1][j] + u[i][j]) * dY / 2.0, 2.0*Dx + (u[i - 1][j] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
								//aW_U = 2.0 * Dx + (u[i - 1][j] + u[i][j]) * dY / (2.0 * 2.0);
								aE_U = Tmax(-(u[i][j] + u[i + 1][j]) * dY / 2.0, Dx - (u[i][j] + u[i + 1][j]) * dY / (2.0 * 2.0), 0.0);
								//aE_U = Dx - (u[i][j] + u[i + 1][j]) * dY / (2.0 * 2.0);
								aS_U = Tmax((v[i][j - 1] + v[i + 1][j - 1]) * dX / 2.0, Dy + (v[i][j - 1] + v[i + 1][j - 1]) * dX / (2.0 * 2.0), 0.0);
								//aS_U = Dy + (v[i][j - 1] + v[i + 1][j - 1]) * dX / (2.0 * 2.0);
								aN_U = Tmax(-(v[i][j] + v[i + 1][j]) * dX / 2.0, Dy - (v[i][j] + v[i + 1][j]) * dX / (2.0 * 2.0), 0.0);
								//aN_U = Tmax(-(v[i][j] + v[i + 1][j]) * dX / 2.0, Dy - (v[i][j] + v[i + 1][j]) * dX / (2.0 * 2.0), 0.0);
								aP_U = aW_U + aE_U + aS_U + aN_U + aP0;
								//aP_U = aW_U + aE_U + aS_U + aN_U + aP0 + (u[i - 1][j] + u[i][j]) * dY / 2.0 - (u[i][j] + u[i + 1][j]) * dY / 2.0 + (v[i][j - 1] + v[i + 1][j - 1]) * dX / 2.0 - (v[i][j] + v[i + 1][j]) * dX / 2.0;
								nrm.SetResult(p, aP_U * u[i][j], aW_U * u[i - 1][j] + aE_U * u[i + 1][j] + aS_U * u[i][j - 1] + aN_U * u[i][j + 1] + aP0 * uu[i][j] + (P[i][j] - P[i + 1][j]) * dY);
								p++;
								if(j==1){
									aW_T = Tmax(u[i - 1][j] * dY, 2.0*Dx / Pr + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_T = 2.0 * Dx / Pr + u[i - 1][j] * dY / 2.0;
									aE_T = Tmax(-u[i][j] * dY, Dx / Pr - u[i][j] * dY / 2.0, 0.0);
									//aE_T = Dx / Pr - u[i][j] * dY / 2.0;
									aS_T = Tmax(v[i][j - 1] * dX, 2.0 * Dy / Pr + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_T = 2.0 * Dy / Pr + v[i][j - 1] * dX / 2.0;
									aN_T = Tmax(-v[i][j] * dX, Dy / Pr - v[i][j] * dX / 2.0, 0.0);
									//aN_T = Dy / Pr - v[i][j] * dX / 2.0;
									aP_T = aW_T + aE_T + aS_T + aN_T + aP0 ;
									//aP_T = aW_T + aE_T + aS_T + aN_T + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									aW_X = Tmax(u[i - 1][j] * dY, 2.0 * Dx / Sc + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_X = 2.0 * Dx / Sc + u[i - 1][j] * dY / 2.0;
									aE_X = Tmax(-u[i][j] * dY, Dx / Sc - u[i][j] * dY / 2.0, 0.0);
									//aE_X = Dx / Sc - u[i][j] * dY / 2.0;
									aS_X = Tmax(v[i][j - 1] * dX, 2.0 * Dy / Sc + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_X = 2.0 * Dy / Sc + v[i][j - 1] * dX / 2.0;
									aN_X = Tmax(-v[i][j] * dX, Dy / Sc - v[i][j] * dX / 2.0, 0.0);
									//aN_X = Dy / Sc - v[i][j] * dX / 2.0;
									aP_X = aW_X + aE_X + aS_X + aN_X + aP0 ;
									//aP_X = aW_X + aE_X + aS_X + aN_X + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									nrm.SetResult(p, aP_T* T[i][j], aW_T* T[i - 1][j] + aE_T * T[i + 1][j] + aS_T * T[i][j - 1] + aN_T * T[i][j + 1] + aP0 * tt[i][j]);
									p++;
									nrm.SetResult(p, aP_X* X[i][j], aW_X* X[i - 1][j] + aE_X * X[i + 1][j] + aS_X * X[i][j - 1] + aN_X * X[i][j + 1] + aP0 * xx[i][j]);
									p++;
								}
								else if (j == M) {
									aW_T = Tmax(u[i - 1][j] * dY, 2.0*Dx / Pr + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_T = 2.0 * Dx / Pr + u[i - 1][j] * dY / 2.0;
									aE_T = Tmax(-u[i][j] * dY, Dx / Pr - u[i][j] * dY / 2.0, 0.0);
									//aE_T = Dx / Pr - u[i][j] * dY / 2.0;
									aS_T = Tmax(v[i][j - 1] * dX, Dy / Pr + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_T = Dy / Pr + v[i][j - 1] * dX / 2.0;
									aN_T = Tmax(-v[i][j] * dX, 2.0 * Dy / Pr - v[i][j] * dX / 2.0, 0.0);
									//aN_T = 2.0 * Dy / Pr - v[i][j] * dX / 2.0;
									aP_T = aW_T + aE_T + aS_T + aN_T + aP0 ;
									//aP_T = aW_T + aE_T + aS_T + aN_T + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									aW_X = Tmax(u[i - 1][j] * dY, 2.0 * Dx / Sc + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_X = 2.0 * Dx / Sc + u[i - 1][j] * dY / 2.0;
									aE_X = Tmax(-u[i][j] * dY, Dx / Sc - u[i][j] * dY / 2.0, 0.0);
									//aE_X = Dx / Sc - u[i][j] * dY / 2.0;
									aS_X = Tmax(v[i][j - 1] * dX, Dy / Sc + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_X = Dy / Sc + v[i][j - 1] * dX / 2.0;
									aN_X = Tmax(-v[i][j] * dX, 2.0 * Dy / Sc - v[i][j] * dX / 2.0, 0.0);
									//aN_X = 2.0 * Dy / Sc - v[i][j] * dX / 2.0;
									aP_X = aW_X + aE_X + aS_X + aN_X + aP0 ;
									//aP_X = aW_X + aE_X + aS_X + aN_X + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									nrm.SetResult(p, aP_T* T[i][j], aW_T* T[i - 1][j] + aE_T * T[i + 1][j] + aS_T * T[i][j - 1] + aN_T * T[i][j + 1] + aP0 * tt[i][j]);
									p++;
									nrm.SetResult(p, aP_X* X[i][j], aW_X* X[i - 1][j] + aE_X * X[i + 1][j] + aS_X * X[i][j - 1] + aN_X * X[i][j + 1] + aP0 * xx[i][j]);
									p++;
								}
								else {
									aW_T = Tmax(u[i - 1][j] * dY, 2.0*Dx / Pr + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_T = 2.0 * Dx / Pr + u[i - 1][j] * dY / 2.0;
									aE_T = Tmax(-u[i][j] * dY, Dx / Pr - u[i][j] * dY / 2.0, 0.0);
									//aE_T = Dx / Pr - u[i][j] * dY / 2.0;
									aS_T = Tmax(v[i][j - 1] * dX, Dy / Pr + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_T = Dy / Pr + v[i][j - 1] * dX / 2.0;
									aN_T = Tmax(-v[i][j] * dX, Dy / Pr - v[i][j] * dX / 2.0, 0.0);
									//aN_T = Dy / Pr - v[i][j] * dX / 2.0;
									aP_T = aW_T + aE_T + aS_T + aN_T + aP0 ;
									//aP_T = aW_T + aE_T + aS_T + aN_T + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									aW_X = Tmax(u[i - 1][j] * dY, 2.0 * Dx / Sc + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_X = 2.0 * Dx / Sc + u[i - 1][j] * dY / 2.0;
									aE_X = Tmax(-u[i][j] * dY, Dx / Sc - u[i][j] * dY / 2.0, 0.0);
									//aE_X = Dx / Sc - u[i][j] * dY / 2.0;
									aS_X = Tmax(v[i][j - 1] * dX, Dy / Sc + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_X = Dy / Sc + v[i][j - 1] * dX / 2.0;
									aN_X = Tmax(-v[i][j] * dX, Dy / Sc - v[i][j] * dX / 2.0, 0.0);
									//aN_X = Dy / Sc - v[i][j] * dX / 2.0;
									aP_X = aW_X + aE_X + aS_X + aN_X + aP0;
									//aP_X = aW_X + aE_X + aS_X + aN_X + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									nrm.SetResult(p, aP_T* T[i][j], aW_T* T[i - 1][j] + aE_T * T[i + 1][j] + aS_T * T[i][j - 1] + aN_T * T[i][j + 1] + aP0 * tt[i][j]);
									p++;
									nrm.SetResult(p, aP_X* X[i][j], aW_X* X[i - 1][j] + aE_X * X[i + 1][j] + aS_X * X[i][j - 1] + aN_X * X[i][j + 1] + aP0 * xx[i][j]);
									p++;
								}
							}
							else if (i == N) {

								aW_U = Tmax((u[i][j] + u[i-1][j]) * dY / 2.0, Dx + (u[i][j] + u[i-1][j]) * dY / (2.0 * 2.0), 0.0);
								aE_U = Tmax(-(u[i-1][j] + u[i][j]) * dY / 2.0, 2.0*Dx - (u[i-1][j] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
								aS_U = Tmax((v[i-1][j - 1] + v[i-1][j - 1]) * dX / 2.0, Dy + (v[i-1][j - 1] + v[i][j - 1]) * dX / (2.0 * 2.0), 0.0);
								aN_U = Tmax(-(v[i-1][j] + v[i][j]) * dX / 2.0, Dy - (v[i-1][j] + v[i][j]) * dX / (2.0 * 2.0), 0.0);
								aP_U = aW_U + aE_U + aS_U + aN_U + aP0;
								//aP_U = aW_U + aE_U + aS_U + aN_U + aP0 + (u[i][j] + u[i - 1][j]) * dY / 2.0 - (u[i - 1][j] + u[i][j]) * dY / 2.0 + (v[i - 1][j - 1] + v[i - 1][j - 1]) * dX / 2.0 - (v[i - 1][j] + v[i][j]) * dX / 2.0;
								nrm.SetResult(p, aP_U * u[i-1][j], aW_U * u[i - 2][j] + aE_U * u[i][j] + aS_U * u[i-1][j - 1] + aN_U * u[i-1][j + 1] + aP0 * uu[i-1][j] + (P[i][j] - P[i+1][j]) * dY);
								p++;
								if (j == 1) {
									aW_T = Tmax(u[i - 1][j] * dY,  Dx / Pr + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_T = Dx / Pr + u[i - 1][j] * dY / 2.0;
									aE_T = Tmax(-u[i][j] * dY, 2.0 * Dx / Pr - u[i][j] * dY / 2.0, 0.0);
									//aE_T = 2.0 * Dx / Pr - u[i][j] * dY / 2.0;
									aS_T = Tmax(v[i][j - 1] * dX, 2.0 * Dy / Pr + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_T = 2.0 * Dy / Pr + v[i][j - 1] * dX / 2.0;
									aN_T = Tmax(-v[i][j] * dX, Dy / Pr - v[i][j] * dX / 2.0, 0.0);
									//aN_T = Dy / Pr - v[i][j] * dX / 2.0;
									aP_T = aW_T + aE_T + aS_T + aN_T + aP0;
									//aP_T = aW_T + aE_T + aS_T + aN_T + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									aW_X = Tmax(u[i - 1][j] * dY,  Dx / Sc + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_X = Dx / Sc + u[i - 1][j] * dY / 2.0;
									aE_X = Tmax(-u[i][j] * dY, 2.0 * Dx / Sc - u[i][j] * dY / 2.0, 0.0);
									//aE_X = 2.0 * Dx / Sc - u[i][j] * dY / 2.0;
									aS_X = Tmax(v[i][j - 1] * dX, 2.0 * Dy / Sc + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_X = 2.0 * Dy / Sc + v[i][j - 1] * dX / 2.0;
									aN_X = Tmax(-v[i][j] * dX, Dy / Sc - v[i][j] * dX / 2.0, 0.0);
									//aN_X = Dy / Sc - v[i][j] * dX / 2.0;
									aP_X = aW_X + aE_X + aS_X + aN_X + aP0 ;
									//aP_X = aW_X + aE_X + aS_X + aN_X + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									nrm.SetResult(p, aP_T * T[i][j], aW_T * T[i - 1][j] + aE_T * T[i + 1][j] + aS_T * T[i][j - 1] + aN_T * T[i][j + 1] + aP0 * tt[i][j]);
									p++;
									nrm.SetResult(p, aP_X * X[i][j], aW_X * X[i - 1][j] + aE_X * X[i + 1][j] + aS_X * X[i][j - 1] + aN_X * X[i][j + 1] + aP0 * xx[i][j]);
									p++;
								}
								else if (j == M) {
									aW_T = Tmax(u[i - 1][j] * dY,  Dx / Pr + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_T = Dx / Pr + u[i - 1][j] * dY / 2.0;
									aE_T = Tmax(-u[i][j] * dY, 2.0 * Dx / Pr - u[i][j] * dY / 2.0, 0.0);
									//aE_T = 2.0 * Dx / Pr - u[i][j] * dY / 2.0;
									aS_T = Tmax(v[i][j - 1] * dX, Dy / Pr + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_T = Dy / Pr + v[i][j - 1] * dX / 2.0;
									aN_T = Tmax(-v[i][j] * dX, 2.0 * Dy / Pr - v[i][j] * dX / 2.0, 0.0);
									//aN_T = 2.0 * Dy / Pr - v[i][j] * dX / 2.0;
									aP_T = aW_T + aE_T + aS_T + aN_T + aP0 ;
									//aP_T = aW_T + aE_T + aS_T + aN_T + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									aW_X = Tmax(u[i - 1][j] * dY,  Dx / Sc + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_X = Dx / Sc + u[i - 1][j] * dY / 2.0;
									aE_X = Tmax(-u[i][j] * dY, 2.0 * Dx / Sc - u[i][j] * dY / 2.0, 0.0);
									//aE_X = 2.0 * Dx / Sc - u[i][j] * dY / 2.0;
									aS_X = Tmax(v[i][j - 1] * dX, Dy / Sc + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_X = Dy / Sc + v[i][j - 1] * dX / 2.0;
									aN_X = Tmax(-v[i][j] * dX, 2.0 * Dy / Sc - v[i][j] * dX / 2.0, 0.0);
									//aN_X = 2.0 * Dy / Sc - v[i][j] * dX / 2.0;
									aP_X = aW_X + aE_X + aS_X + aN_X + aP0 ;
									//aP_X = aW_X + aE_X + aS_X + aN_X + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									nrm.SetResult(p, aP_T * T[i][j], aW_T * T[i - 1][j] + aE_T * T[i + 1][j] + aS_T * T[i][j - 1] + aN_T * T[i][j + 1] + aP0 * tt[i][j]);
									p++;
									nrm.SetResult(p, aP_X * X[i][j], aW_X * X[i - 1][j] + aE_X * X[i + 1][j] + aS_X * X[i][j - 1] + aN_X * X[i][j + 1] + aP0 * xx[i][j]);
									p++;
								}
								else {
									aW_T = Tmax(u[i - 1][j] * dY,  Dx / Pr + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_T = Dx / Pr + u[i - 1][j] * dY / 2.0;
									aE_T = Tmax(-u[i][j] * dY, 2.0 * Dx / Pr - u[i][j] * dY / 2.0, 0.0);
									//aE_T = 2.0 * Dx / Pr - u[i][j] * dY / 2.0;
									aS_T = Tmax(v[i][j - 1] * dX, Dy / Pr + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_T = Dy / Pr + v[i][j - 1] * dX / 2.0;
									aN_T = Tmax(-v[i][j] * dX, Dy / Pr - v[i][j] * dX / 2.0, 0.0);
									//aN_T = Dy / Pr - v[i][j] * dX / 2.0;
									aP_T = aW_T + aE_T + aS_T + aN_T + aP0 ;
									//aP_T = aW_T + aE_T + aS_T + aN_T + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									aW_X = Tmax(u[i - 1][j] * dY,  Dx / Sc + u[i - 1][j] * dY / 2.0, 0.0);
									//aW_X = Dx / Sc + u[i - 1][j] * dY / 2.0;
									aE_X = Tmax(-u[i][j] * dY, 2.0 * Dx / Sc - u[i][j] * dY / 2.0, 0.0);
									//aE_X = 2.0 * Dx / Sc - u[i][j] * dY / 2.0;
									aS_X = Tmax(v[i][j - 1] * dX, Dy / Sc + v[i][j - 1] * dX / 2.0, 0.0);
									//aS_X = Dy / Sc + v[i][j - 1] * dX / 2.0;
									aN_X = Tmax(-v[i][j] * dX, Dy / Sc - v[i][j] * dX / 2.0, 0.0);
									//aN_X = Dy / Sc - v[i][j] * dX / 2.0;
									aP_X = aW_X + aE_X + aS_X + aN_X + aP0 ;
									//aP_X = aW_X + aE_X + aS_X + aN_X + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
									nrm.SetResult(p, aP_T * T[i][j], aW_T * T[i - 1][j] + aE_T * T[i + 1][j] + aS_T * T[i][j - 1] + aN_T * T[i][j + 1] + aP0 * tt[i][j]);
									p++;
									nrm.SetResult(p, aP_X * X[i][j], aW_X * X[i - 1][j] + aE_X * X[i + 1][j] + aS_X * X[i][j - 1] + aN_X * X[i][j + 1] + aP0 * xx[i][j]);
									p++;
								}
							}
							else if (j == 1) {
								aW_U = Tmax((u[i - 1][j] + u[i][j]) * dY / 2.0, Dx + (u[i - 1][j] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
								aE_U = Tmax(-(u[i][j] + u[i + 1][j]) * dY / 2.0, Dx - (u[i][j] + u[i + 1][j]) * dY / (2.0 * 2.0), 0.0);
								aS_U = Tmax((v[i][j - 1] + v[i + 1][j - 1]) * dX / 2.0, Dy*2.0 + (v[i][j - 1] + v[i + 1][j - 1]) * dX / (2.0 * 2.0), 0.0);
								aN_U = Tmax(-(v[i][j] + v[i + 1][j]) * dX / 2.0, Dy - (v[i][j] + v[i + 1][j]) * dX / (2.0 * 2.0), 0.0);
								aP_U = aW_U + aE_U + aS_U + aN_U + aP0 ;
								//aP_U = aW_U + aE_U + aS_U + aN_U + aP0 + (u[i - 1][j] + u[i][j]) * dY / 2.0 - (u[i][j] + u[i + 1][j]) * dY / 2.0 + (v[i][j - 1] + v[i + 1][j - 1]) * dX / 2.0 - (v[i][j] + v[i + 1][j]) * dX / 2.0;
								nrm.SetResult(p, aP_U* u[i][j], aW_U* u[i - 1][j] + aE_U * u[i + 1][j] + aS_U * u[i][j - 1] + aN_U * u[i][j + 1] + aP0 * uu[i][j] + (P[i][j] - P[i + 1][j]) * dY);
								p++;
								aW_T = Tmax(u[i - 1][j] * dY, Dx / Pr + u[i - 1][j] * dY / 2.0, 0.0);
								//aW_T = Dx / Pr + u[i - 1][j] * dY / 2.0;
								aE_T = Tmax(-u[i][j] * dY, Dx / Pr - u[i][j] * dY / 2.0, 0.0);
								//aE_T = Dx / Pr - u[i][j] * dY / 2.0;
								aS_T = Tmax(v[i][j - 1] * dX, Dy * 2.0 / Pr + v[i][j - 1] * dX / 2.0, 0.0);
								//aS_T = Dy * 2.0 / Pr + v[i][j - 1] * dX / 2.0;
								aN_T = Tmax(-v[i][j] * dX, Dy / Pr - v[i][j] * dX / 2.0, 0.0);
								//aN_T = Dy / Pr - v[i][j] * dX / 2.0;
								aP_T = aW_T + aE_T + aS_T + aN_T + aP0 ;
								//aP_T = aW_T + aE_T + aS_T + aN_T + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
								aW_X = Tmax(u[i - 1][j] * dY, Dx / Sc + u[i - 1][j] * dY / 2.0, 0.0);
								//aW_X = Dx / Sc + u[i - 1][j] * dY / 2.0;
								aE_X = Tmax(-u[i][j] * dY, Dx / Sc - u[i][j] * dY / 2.0, 0.0);
								//aE_X = Dx / Sc - u[i][j] * dY / 2.0;
								aS_X = Tmax(v[i][j - 1] * dX, Dy * 2.0 / Sc + v[i][j - 1] * dX / 2.0, 0.0);
								//aS_X = Dy * 2.0 / Sc + v[i][j - 1] * dX / 2.0;
								aN_X = Tmax(-v[i][j] * dX, Dy / Sc - v[i][j] * dX / 2.0, 0.0);
								//aN_X = Dy / Sc - v[i][j] * dX / 2.0;
								aP_X = aW_X + aE_X + aS_X + aN_X + aP0 ;
								//aP_X = aW_X + aE_X + aS_X + aN_X + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
								nrm.SetResult(p, aP_T* T[i][j], aW_T* T[i - 1][j] + aE_T * T[i + 1][j] + aS_T * T[i][j - 1] + aN_T * T[i][j + 1] + aP0 * tt[i][j]);
								p++;
								nrm.SetResult(p, aP_X* X[i][j], aW_X* X[i - 1][j] + aE_X * X[i + 1][j] + aS_X * X[i][j - 1] + aN_X * X[i][j + 1] + aP0 * xx[i][j]);
								p++;
							}
							else if (j == M) {
							aW_U = Tmax((u[i - 1][j] + u[i][j]) * dY / 2.0, Dx + (u[i - 1][j] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
							aE_U = Tmax(-(u[i][j] + u[i + 1][j]) * dY / 2.0, Dx - (u[i][j] + u[i + 1][j]) * dY / (2.0 * 2.0), 0.0);
							aS_U = Tmax((v[i][j - 1] + v[i + 1][j - 1]) * dX / 2.0, Dy + (v[i][j - 1] + v[i + 1][j - 1]) * dX / (2.0 * 2.0), 0.0);
							aN_U = Tmax(-(v[i][j] + v[i + 1][j]) * dX / 2.0, Dy * 2.0 - (v[i][j] + v[i + 1][j]) * dX / (2.0 * 2.0), 0.0);
							aP_U = aW_U + aE_U + aS_U + aN_U + aP0 ;
							//aP_U = aW_U + aE_U + aS_U + aN_U + aP0 + (u[i - 1][j] + u[i][j]) * dY / 2.0 - (u[i][j] + u[i + 1][j]) * dY / 2.0 + (v[i][j - 1] + v[i + 1][j - 1]) * dX / 2.0 - (v[i][j] + v[i + 1][j]) * dX / 2.0;
							nrm.SetResult(p, aP_U* u[i][j], aW_U* u[i - 1][j] + aE_U * u[i + 1][j] + aS_U * u[i][j - 1] + aN_U * u[i][j + 1] + aP0 * uu[i][j] + (P[i][j] - P[i + 1][j]) * dY);
							p++;
							aW_T = Tmax(u[i - 1][j] * dY, Dx / Pr + u[i - 1][j] * dY / 2.0, 0.0);
							//aW_T = Dx / Pr + u[i - 1][j] * dY / 2.0;
							aE_T = Tmax(-u[i][j] * dY, Dx / Pr - u[i][j] * dY / 2.0, 0.0);
							//aE_T = Dx / Pr - u[i][j] * dY / 2.0;
							aS_T = Tmax(v[i][j - 1] * dX, Dy / Pr + v[i][j - 1] * dX / 2.0, 0.0);
							//aS_T = Dy / Pr + v[i][j - 1] * dX / 2.0;
							aN_T = Tmax(-v[i][j] * dX, Dy * 2.0 / Pr - v[i][j] * dX / 2.0, 0.0);
							//aN_T = Dy * 2.0 / Pr - v[i][j] * dX / 2.0;
							aP_T = aW_T + aE_T + aS_T + aN_T + aP0 ;
							//aP_T = aW_T + aE_T + aS_T + aN_T + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
							aW_X = Tmax(u[i - 1][j] * dY, Dx / Sc + u[i - 1][j] * dY / 2.0, 0.0);
							//aW_X = Dx / Sc + u[i - 1][j] * dY / 2.0;
							aE_X = Tmax(-u[i][j] * dY, Dx / Sc - u[i][j] * dY / 2.0, 0.0);
							//aE_X = Dx / Sc - u[i][j] * dY / 2.0;
							aS_X = Tmax(v[i][j - 1] * dX, Dy / Sc + v[i][j - 1] * dX / 2.0, 0.0);
							//aS_X = Dy / Sc + v[i][j - 1] * dX / 2.0;
							aN_X = Tmax(-v[i][j] * dX, Dy * 2.0 / Sc - v[i][j] * dX / 2.0, 0.0);
							//aN_X = Dy * 2.0 / Sc - v[i][j] * dX / 2.0;
							aP_X = aW_X + aE_X + aS_X + aN_X + aP0 ;
							//aP_X = aW_X + aE_X + aS_X + aN_X + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
							nrm.SetResult(p, aP_T* T[i][j], aW_T* T[i - 1][j] + aE_T * T[i + 1][j] + aS_T * T[i][j - 1] + aN_T * T[i][j + 1] + aP0 * tt[i][j]);
							p++;
							nrm.SetResult(p, aP_X* X[i][j], aW_X* X[i - 1][j] + aE_X * X[i + 1][j] + aS_X * X[i][j - 1] + aN_X * X[i][j + 1] + aP0 * xx[i][j]);
							p++;
							}
							else {
								aW_U = Tmax((u[i - 1][j] + u[i][j]) * dY / 2.0, Dx + (u[i - 1][j] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
								aE_U = Tmax(-(u[i][j] + u[i + 1][j]) * dY / 2.0, Dx - (u[i][j] + u[i + 1][j]) * dY / (2.0 * 2.0), 0.0);
								aS_U = Tmax((v[i][j - 1] + v[i + 1][j - 1]) * dX / 2.0, Dy + (v[i][j - 1] + v[i + 1][j - 1]) * dX / (2.0 * 2.0), 0.0);
								aN_U = Tmax(-(v[i][j] + v[i + 1][j]) * dX / 2.0, Dy - (v[i][j] + v[i + 1][j]) * dX / (2.0 * 2.0), 0.0);
								aP_U = aW_U + aE_U + aS_U + aN_U + aP0 ;
								//aP_U = aW_U + aE_U + aS_U + aN_U + aP0 + (u[i - 1][j] + u[i][j]) * dY / 2.0 - (u[i][j] + u[i + 1][j]) * dY / 2.0 + (v[i][j - 1] + v[i + 1][j - 1]) * dX / 2.0 - (v[i][j] + v[i + 1][j]) * dX / 2.0;
								nrm.SetResult(p, aP_U * u[i][j], aW_U * u[i - 1][j] + aE_U * u[i + 1][j] + aS_U * u[i][j - 1] + aN_U * u[i][j + 1] + aP0 * uu[i][j] + (P[i][j] - P[i + 1][j]) * dY);
								p++;

								aW_T = Tmax(u[i - 1][j] * dY, Dx / Pr + u[i - 1][j] * dY / 2.0, 0.0);
								//aW_T = Dx / Pr + u[i - 1][j] * dY / 2.0;
								aE_T = Tmax(-u[i][j] * dY, Dx / Pr - u[i][j] * dY / 2.0, 0.0);
								//aE_T = Dx / Pr - u[i][j] * dY / 2.0;
								aS_T = Tmax(v[i][j - 1] * dX, Dy / Pr + v[i][j - 1] * dX / 2.0, 0.0);
								//aS_T = Dy / Pr + v[i][j - 1] * dX / 2.0;
								aN_T = Tmax(-v[i][j] * dX, Dy / Pr - v[i][j] * dX / 2.0, 0.0);
								//aN_T = Dy / Pr - v[i][j] * dX / 2.0;
								aP_T = aW_T + aE_T + aS_T + aN_T + aP0 ;
								//aP_T = aW_T + aE_T + aS_T + aN_T + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
								aW_X = Tmax(u[i - 1][j] * dY, Dx / Sc + u[i - 1][j] * dY / 2.0, 0.0);
								//aW_X = Dx / Sc + u[i - 1][j] * dY / 2.0;
								aE_X = Tmax(-u[i][j] * dY, Dx / Sc - u[i][j] * dY / 2.0, 0.0);
								//aE_X = Dx / Sc - u[i][j] * dY / 2.0;
								aS_X = Tmax(v[i][j - 1] * dX, Dy / Sc + v[i][j - 1] * dX / 2.0, 0.0);
								//aS_X = Dy / Sc + v[i][j - 1] * dX / 2.0;
								aN_X = Tmax(-v[i][j] * dX, Dy / Sc - v[i][j] * dX / 2.0, 0.0);
								//aN_X = Dy / Sc - v[i][j] * dX / 2.0;
								aP_X = aW_X + aE_X + aS_X + aN_X + aP0 ;
								//aP_X = aW_X + aE_X + aS_X + aN_X + aP0 + u[i - 1][j] * dY - u[i][j] * dY + v[i][j - 1] * dX - v[i][j] * dX;
								nrm.SetResult(p, aP_T * T[i][j], aW_T * T[i - 1][j] + aE_T * T[i + 1][j] + aS_T * T[i][j - 1] + aN_T * T[i][j + 1] + aP0 * tt[i][j]);
								p++;
								nrm.SetResult(p, aP_X * X[i][j], aW_X * X[i - 1][j] + aE_X * X[i + 1][j] + aS_X * X[i][j - 1] + aN_X * X[i][j + 1] + aP0 * xx[i][j]);
								p++;
							}
						}

					}
				}

				for (j = 0; j <= M; j++) {
					for (i = 0; i <= N + 1; i++) {
					
						if (i == 0) {
							nrm.SetResult(p, v[i][j], 0.0);
							p++;
						}
						else if (i == N+1) {
							nrm.SetResult(p, v[i][j], 0.0);
							p++;
						}
						else if (j == 0) {
							nrm.SetResult(p, v[i][j], 0.0);
							p++;
						}
						else if (j == M) {
							nrm.SetResult(p, v[i][j], 0.0);
							p++;
						}
						else {
							if (i == 1) {
								aW_V = Tmax((u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0, 2.0*Dx + (u[i - 1][j + 1] + u[i - 1][j]) * dY / (2.0 * 2.0), 0.0);
								aE_V = Tmax(-(u[i][j + 1] + u[i][j]) * dY / 2.0, Dx - (u[i][j + 1] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
								aS_V = Tmax((v[i][j] + v[i][j - 1]) * dX / 2.0, Dy + (v[i][j] + v[i][j - 1]) * dX / (2.0 * 2.0), 0.0);
								aN_V = Tmax(-(v[i][j] + v[i][j + 1]) * dX / 2.0, Dy - (v[i][j] + v[i][j + 1]) * dX / (2.0 * 2.0), 0.0);
								aP_V = aW_V + aE_V + aS_V + aN_V + aP0;
								//aP_V = aW_V + aE_V + aS_V + aN_V + aP0 + (u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0 - (u[i][j + 1] + u[i][j]) * dY / 2.0 + (v[i][j] + v[i][j - 1]) * dX / 2.0 - (v[i][j] + v[i][j + 1]) * dX / 2.0;
								nrm.SetResult(p, aP_V* v[i][j], aW_V* v[i - 1][j] + aE_V * v[i + 1][j] + aS_V * v[i][j - 1] + aN_V * v[i][j + 1] + aP0 * vv[i][j] + (P[i][j] - P[i][j + 1]) * dX);
								p++;
							}
							else if (i == N) {
								aW_V = Tmax((u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0, Dx + (u[i - 1][j + 1] + u[i - 1][j]) * dY / (2.0 * 2.0), 0.0);
								aE_V = Tmax(-(u[i][j + 1] + u[i][j]) * dY / 2.0, 2.0 * Dx - (u[i][j + 1] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
								aS_V = Tmax((v[i][j] + v[i][j - 1]) * dX / 2.0, Dy + (v[i][j] + v[i][j - 1]) * dX / (2.0 * 2.0), 0.0);
								aN_V = Tmax(-(v[i][j] + v[i][j + 1]) * dX / 2.0, Dy - (v[i][j] + v[i][j + 1]) * dX / (2.0 * 2.0), 0.0);
								aP_V = aW_V + aE_V + aS_V + aN_V + aP0; 
								//aP_V = aW_V + aE_V + aS_V + aN_V + aP0 + (u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0 - (u[i][j + 1] + u[i][j]) * dY / 2.0 + (v[i][j] + v[i][j - 1]) * dX / 2.0 - (v[i][j] + v[i][j + 1]) * dX / 2.0;
								nrm.SetResult(p, aP_V* v[i][j], aW_V* v[i - 1][j] + aE_V * v[i + 1][j] + aS_V * v[i][j - 1] + aN_V * v[i][j + 1] + aP0 * vv[i][j] + (P[i][j] - P[i][j + 1]) * dX);
								p++;
							}
							else if (j == 1) {
								aW_V = Tmax((u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0, Dx + (u[i - 1][j + 1] + u[i - 1][j]) * dY / (2.0 * 2.0), 0.0);
								aE_V = Tmax(-(u[i][j + 1] + u[i][j]) * dY / 2.0, Dx - (u[i][j + 1] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
								aS_V = Tmax((v[i][j] + v[i][j - 1]) * dX / 2.0, Dy + (v[i][j] + v[i][j - 1]) * dX / (2.0 * 2.0), 0.0);
								aN_V = Tmax(-(v[i][j] + v[i][j + 1]) * dX / 2.0, Dy - (v[i][j] + v[i][j + 1]) * dX / (2.0 * 2.0), 0.0);
								aP_V = aW_V + aE_V + aS_V + aN_V + aP0; 
								//aP_V = aW_V + aE_V + aS_V + aN_V + aP0 + (u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0 - (u[i][j + 1] + u[i][j]) * dY / 2.0 + (v[i][j] + v[i][j - 1]) * dX / 2.0 - (v[i][j] + v[i][j + 1]) * dX / 2.0;
								nrm.SetResult(p, aP_V* v[i][j], aW_V* v[i - 1][j] + aE_V * v[i + 1][j] + aS_V * v[i][j - 1] + aN_V * v[i][j + 1] + aP0 * vv[i][j] + (P[i][j] - P[i][j + 1]) * dX);
								p++;
							}
							else if (j == M) {
								aW_V = Tmax((u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0, Dx + (u[i - 1][j + 1] + u[i - 1][j]) * dY / (2.0 * 2.0), 0.0);
								aE_V = Tmax(-(u[i][j + 1] + u[i][j]) * dY / 2.0, Dx - (u[i][j + 1] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
								aS_V = Tmax((v[i][j] + v[i][j - 1]) * dX / 2.0, Dy + (v[i][j] + v[i][j - 1]) * dX / (2.0 * 2.0), 0.0);
								aN_V = Tmax(-(v[i][j] + v[i][j + 1]) * dX / 2.0, Dy - (v[i][j] + v[i][j + 1]) * dX / (2.0 * 2.0), 0.0);
								aP_V = aW_V + aE_V + aS_V + aN_V + aP0; 
								//aP_V = aW_V + aE_V + aS_V + aN_V + aP0 + (u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0 - (u[i][j + 1] + u[i][j]) * dY / 2.0 + (v[i][j] + v[i][j - 1]) * dX / 2.0 - (v[i][j] + v[i][j + 1]) * dX / 2.0;
								nrm.SetResult(p, aP_V* v[i][j], aW_V* v[i - 1][j] + aE_V * v[i + 1][j] + aS_V * v[i][j - 1] + aN_V * v[i][j + 1] + aP0 * vv[i][j] + (P[i][j] - P[i][j + 1]) * dX);
								p++;
							}
							else {
								aW_V = Tmax((u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0, Dx + (u[i - 1][j + 1] + u[i - 1][j]) * dY / (2.0 * 2.0), 0.0);
								aE_V = Tmax(-(u[i][j + 1] + u[i][j]) * dY / 2.0, Dx - (u[i][j + 1] + u[i][j]) * dY / (2.0 * 2.0), 0.0);
								aS_V = Tmax((v[i][j] + v[i][j - 1]) * dX / 2.0, Dy + (v[i][j] + v[i][j - 1]) * dX / (2.0 * 2.0), 0.0);
								aN_V = Tmax(-(v[i][j] + v[i][j + 1]) * dX / 2.0, Dy - (v[i][j] + v[i][j + 1]) * dX / (2.0 * 2.0), 0.0);
								aP_V = aW_V + aE_V + aS_V + aN_V + aP0; 
								//aP_V = aW_V + aE_V + aS_V + aN_V + aP0 + (u[i - 1][j + 1] + u[i - 1][j]) * dY / 2.0 - (u[i][j + 1] + u[i][j]) * dY / 2.0 + (v[i][j] + v[i][j - 1]) * dX / 2.0 - (v[i][j] + v[i][j + 1]) * dX / 2.0;
								nrm.SetResult(p, aP_V * v[i][j], aW_V * v[i - 1][j] + aE_V * v[i + 1][j] + aS_V * v[i][j - 1] + aN_V * v[i][j + 1] + aP0 * vv[i][j] + (P[i][j] - P[i][j + 1]) * dX);
								p++;
							}
						}

					}
				}
				//fout << 0 << endl;
				/////////////////////






			}
		}
		if (t % 1 == 0) {
			for (j = 0; j <= M + 1; j++) {
				for (i = 0; i <= N; i++) {
					fout << u[i][j] << ",";
				}
				fout << 0 << endl;
			}

			for (j = 0; j <= M; j++) {
				for (i = 0; i <= N + 1; i++) {
					fout << v[i][j] << ",";
				}
				fout << 0 << endl;
			}

			//for (j = 0; j <= M + 1; j++) {
				//for (i = 0; i <= N + 1; i++) {
					//fout << P[i][j] << ",";
				//}
				//fout << 0 << endl;
			//}

			for (j = 0; j <= M + 1; j++) {
				for (i = 0; i <= N + 1; i++) {
					fout << T[i][j] << ",";
				}
				fout << 0 << endl;
			}
			for (j = 0; j <= M + 1; j++) {
				for (i = 0; i <= N + 1; i++) {
					fout << X[i][j] << ",";
				}
				fout << 0 << endl;
			}

			mm_w = m_w;
			m_w = 0;

			////////////////////////////////Double check!!!
			for (j = 1; j <= M; j++) {
				for (i = 1; i <= N; i++) {
					//m_w += libr.sc_rho_XT(X[i][j] * (X_in - X_e) + X_e, T[i][j] * T_e) * (h / M) * (L / N) * (X_in - (X[i][j] * (X_in - X_e) + X_e));
					m_w += libr.sc_rho_XT(X[i][j] * (X_in) ,T[i][j]*T_e) * (h / M) * (L / N) * (X_in-(X[i][j]*X_in));  
				}
			}
			fout << m_w-mm_w << ",";
			fout << m_w << endl;                 //////[kg/m]=[g/mm]
		}
		for (i = 0; i <= N + 1; i++) {
			for (j = 0; j <= M + 1; j++) {
				uu[i][j] = u[i][j];
				vv[i][j] = v[i][j];
				tt[i][j] = T[i][j];
				xx[i][j] = X[i][j];
				PP[i][j] = P[i][j];
			}
		}
	}
}*/
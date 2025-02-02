/*
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "CNewtonRaphsonMethod.h"


int main()
{
	int p, i, j, N, M, X, n;
	N = 8000;/////
	M = 1;/////
	X = 200;

	//double	A, B, D, Dh, d, st, slp, Sgla, Sgls, sgla, sgls, ug, ugf, ul, ulf, ugla, ugls, Ug, Ul, Up, Ugla, Ugls, tet, L, G, rhog, rhol, rhop, myug, myul;
	//double  Mi, r, rg, rl, Rel, Reg, Re, g, CA,  X, Xi, P, X0, f, E, W, WD, We, Weg, Wel, beta, Alf, Alff, PI, m, eps;
	
	double  qout, Qout, qroom, Qroom, L, dx, dt, t, Ti, Troom, Tout;
	double  As, a, k, V, rho, hout, hroom, c;
	double  Tp[300], T[300];
	
	ofstream fout("1DwallConduction.csv");
	if (!fout) {
		cout << "ファイルをオープンできませんでした" << endl;
		return 1;
	}
	else
		cout << "ファイルをオープンしました" << endl;

	//PI = (6 * asin(0.5));//π

////////////////////Simulation condition////////////////////

	t = 0.0;				///////time [s]
	dt = 25.0;				///////time step [s]
	Ti = 35.0;				///////initial wall temperature [oC]
	Tout = 35.0;			///////external ambient temperature [oC]
	Troom = 27.0;			///////room temperature [oC]
	L = 0.2;				///////wall thickness [m]
	dx = L / double(X);		///////mesh size [m]

	///////////////wall properties//////////////////////////

	V = 12.890;				///////vlume [m3]
	k = 1.5;				///////thermal conductivity [J/m-K]
	rho = 2000.0;			///////density [kg/m3]
	c = 800.0;				///////specific heat [J/kg-K]
	hout = 5.0;				///////external heat transfer coefficeint [W/m2-K]
	hroom = 5.0;			///////internal heat transfer coefficeint [W/m2-K]
	a = k / rho / c;		///////thermal diffusivity [m2/s]
	As = 100.0;				///////heat transfer area [m2]

	fout << t << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";

	///////////////temperature profile initialization////////
	for (p = 0; p <= X; p++) {

		T[p] = Ti;
		Tp[p] = Ti;
		//T[p] = 35.0 - 7 / L * p * dx;
		//Tp[p] = 35.0 - 7 / L * p * dx;
		fout << T[p] << ",";

	}

	fout << 0 << endl;


	for (i = 1; i <= N; i++) {///////////////

		t = t+dt;
		fout << t << ",";
		n = i * dt / 1800;
		if (n % 2 == 0)
			Troom = 27.0;
		else Troom = 20.0;
		//Troom = 24.0 + 8.0 * sin(2.0 * 3.14 * dt * i / 1800);


		
		CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
		mnm.setup(2000, 1 + 1e-12, 1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

		p = 0;

		for (p = 0; p <= X; p++) {
			mnm.setValue(p, T[p]);
			
		}

		//			mnm.setValue( p , hmin );
		//		p++;

		mnm.setAcc(0.7);				//加速度勾配の入力　※なくても可

		mnm.initial();					//計算を開始する前に必ず初期化してください
		for (mnm.main_loop_init(); mnm.main_loop_check(); mnm.main_loop_reinit()) {		// おまじない
			for (mnm.sub_loop_init(); mnm.sub_loop_check(); mnm.sub_loop_reinit()) {	// おまじない

				p = 0;
				//値の設定
				for (p = 0; p <= X; p++) {

					T[p] = mnm.getValue(p);
					
				
				}

				//			hmin = mnm.getValue(p);
				//		p++;

				p = 0;

				///////////////outdoor boundary condition////////
				//mnm.setError(p, 0.0, -k*(T[p + 1] - T[p]) / (dx) - hout * (Tout - T[p]));
			
				mnm.setError(p, 0.0, T[p] -Tout);
				p++;

					for (p = 1; p <= X - 1; p++) {

						///////////////1-D energy transport equation////////
						mnm.setError(p, 0.0, (T[p+1]-2*T[p]+T[p-1])/(dx*dx)-(1/a)*(T[p]-Tp[p])/dt);
					
					}
				
				///////////////room-side boundary condition////////
				mnm.setError(p, 0.0, k * (T[p] - T[p - 1]) / (dx)- hroom * (Troom - T[p]));
				p++;

				//			mnm.setError( p ,  );
				//		p++;


						// mnm.prt();				//エラー表示
						mnm.prt_sum();			//エラーの合計を表示

			}
		}



		qroom= -k * (T[X] - T[X - 1]) / (dx);
		fout << qroom << ",";
		Qroom = -k *As* (T[X] - T[X - 1]) / (dx);
		fout << Qroom << ",";

		qout = -k * (T[1] - T[0]) / (dx);
		fout << qout << ",";
		Qout = -k * As * (T[1] - T[0]) / (dx);
		fout << Troom << ",";

		///////////////assign next-time-step conditions////////
		for (p = 0; p <= X; p++) {

			Tp[p] = T[p];
			fout << T[p] << ",";

		}

		fout << 0 << endl;
		


	}




}*/

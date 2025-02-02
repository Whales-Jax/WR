/*#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "CNewtonRaphsonMethod.h"


int main()
{
	int p, i, j, N, Xtot, X1, X2, X3, X4, M, n;
	N = 5000;/////
	M = 1;/////
	X1 = 60;
	X2 = 25;
	X3 = 67;
	X4 = 22;
	Xtot = X1 + X2 + X3 + X4;


	//double	A, B, D, Dh, d, st, slp, Sgla, Sgls, sgla, sgls, ug, ugf, ul, ulf, ugla, ugls, Ug, Ul, Up, Ugla, Ugls, tet, L, G, rhog, rhol, rhop, myug, myul;
	//double  Mi, r, rg, rl, Rel, Reg, Re, g, CA,  X, Xi, P, X0, f, E, W, WD, We, Weg, Wel, beta, Alf, Alff, PI, m, eps;

	double  QQout, QQroom, qout, Qout, qroom, Qroom, Ltot, L, L1, L2, L3, L4, dx1, dx2, dx3, dx4, dt, t, Ti, Troom, Tset, Tout;

	L1 = 0.06;				///////wall thickness [m]
	dx1 = L1 / double(X1);		///////mesh size [m]
	L2 = 0.025;				///////wall thickness [m]
	dx2 = L2 / double(X2);		///////mesh size [m]
	L3 = 0.067;				///////wall thickness [m]
	dx3 = L3 / double(X3);		///////mesh size [m]
	L4 = 0.022;				///////wall thickness [m]
	dx4 = L4 / double(X4);		///////mesh size [m]
	Ltot = L1 + L2 + L3 + L4;

	t = 0.0;				///////time [s]
	dt = 100.0;				///////time step [s]


	double  As, a1, k1, V, rho1, hout, hroom, c1, a2, k2, rho2, c2, a3, k3, rho3, c3, a4, k4, rho4, c4;
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

	Ti = 35.0;				///////initial wall temperature [oC]
	Tout = 35.0;			///////external ambient temperature [oC]
	Tset = 21.0;			///////room temperature [oC]

	///////////////wall properties//////////////////////////

	V = 12.890;				///////vlume [m3]
	k1 = 1.5;				///////thermal conductivity [J/m-K]
	rho1 = 2000.0;			///////density [kg/m3]
	c1 = 800.0;				///////specific heat [J/kg-K]
	k2 = 0.023;				///////thermal conductivity [J/m-K]
	rho2 = 35.0;			///////density [kg/m3]
	c2 = 1714.29;				///////specific heat [J/kg-K]
	k3 = 0.0242;				///////thermal conductivity [J/m-K]
	rho3 = 1.225;			///////density [kg/m3]
	c3 = 1006.43;				///////specific heat [J/kg-K]
	k4 = 0.22;				///////thermal conductivity [J/m-K]
	rho4 = 750.0;			///////density [kg/m3]
	c4 = 1106.67;				///////specific heat [J/kg-K]

	//k = 0.023;				///////thermal conductivity [J/m-K]
	//rho = 35.0;			///////density [kg/m3]
	//c = 1714.3;				///////specific heat [J/kg-K]
	hout = 5.0;				///////external heat transfer coefficeint [W/m2-K]
	hroom = 5.0;			///////internal heat transfer coefficeint [W/m2-K]
	a1 = k1 / rho1 / c1;		///////thermal diffusivity [m2/s]
	a2 = k2 / rho2 / c2;		///////thermal diffusivity [m2/s]
	a3 = k3 / rho3 / c3;		///////thermal diffusivity [m2/s]
	a4 = k4 / rho4 / c4;		///////thermal diffusivity [m2/s]
	As = 96.0;				///////heat transfer area [m2]

	fout << t << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";

	///////////////temperature profile initialization////////
	for (p = 0; p <= Xtot; p++) {

		T[p] = Tout;
		Tp[p] = Tout;
		//T[p] = 32.0 - 2.0 / L * p * dx;
		//Tp[p] = 32.0 - 2.0 / L * p * dx;
		//fout << T[p] << ",";

	}
	
	L = 0;
	fout << L << ",";
	for (p = 1; p <= X1; p++) {

		L = L+dx1;
		fout << L << ",";

	}
	for (p = 1; p <= X2; p++) {

		L = L + dx2;
		fout << L << ",";

	}
	for (p = 1; p <= X3; p++) {

		L = L + dx3;
		fout << L << ",";

	}
	for (p = 1; p <= X4; p++) {

		L = L + dx4;
		fout << L << ",";

	}

	fout << 0 << endl;


	for (i = 1; i <= N; i++) {///////////////

		t = t + dt;
		fout << t << ",";
		n = i * dt / 1800;
		//if (n % 2 == 0)
		Troom = Tset;
		//else Troom = 20.0;
		//Troom = Tset + 1.0 * sin(2.0 * 3.14 * dt * double(i) / 6400.0);



		CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
		mnm.setup(2000, 1 + 1e-12, 1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

		p = 0;

		for (p = 0; p <= Xtot; p++) {
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
				for (p = 0; p <= Xtot; p++) {

					T[p] = mnm.getValue(p);


				}

				//			hmin = mnm.getValue(p);
				//		p++;

				p = 0;

				///////////////outdoor boundary condition////////
				mnm.setError(p, 0.0, hout * (Tout - T[p]) + k1 * (T[p + 1] - T[p]) / dx1 - rho1 * c1 * dx1 / 2.0 * (T[p] - Tp[p]) / dt);

				//mnm.setError(p, 0.0, T[p] - Tout);
				p++;

				for (p = 1; p <= X1-1; p++) {

					///////////////1-D energy transport equation////////
					//mnm.setError(p, 0.0, (T[p + 1] - 2 * T[p] + T[p - 1]) / (dx * dx) - (1 / a) * (T[p] - Tp[p]) / dt);
					//mnm.setError(p, 0.0, k1 / dx1 * T[p - 1] - (2 * k1 / dx1 + rho1 * c1 * dx1 / dt) * T[p] + k1 / dx1 * T[p + 1] + rho1 * c1 * dx1 / dt * Tp[p]);
					mnm.setError(p, 0.0, - k1 / dx1 * (T[p]-T[p - 1]) - (- k1 / dx1  * (T[p + 1] - T[p])) - rho1 * c1 * dx1 / dt * (T[p] - Tp[p]));

				}
				///////////////1st layer boundary condition////////
				mnm.setError(p, 0.0, -k1 / dx1 * (T[p] - T[p - 1]) - (-k2 / dx2 * (T[p + 1] - T[p])) - (rho1 * c1 * dx1 / 2 + rho2 * c2 * dx2 / 2) / dt * (T[p] - Tp[p]));

				//mnm.setError(p, 0.0, T[p] - Tout);
				p++;

				for (p = X1+1; p <= Xtot-X4-X3-1; p++) {

					///////////////1-D energy transport equation////////
					//mnm.setError(p, 0.0, (T[p + 1] - 2 * T[p] + T[p - 1]) / (dx * dx) - (1 / a) * (T[p] - Tp[p]) / dt);
					mnm.setError(p, 0.0, k2 / dx2 * T[p - 1] - (2 * k2 / dx2 + rho2 * c2 * dx2 / dt) * T[p] + k2 / dx2 * T[p + 1] + rho2 * c2 * dx2 / dt * Tp[p]);

				}
				///////////////2nd layer boundary condition////////
				mnm.setError(p, 0.0, -k2 / dx2 * (T[p] - T[p - 1]) - (-k3 / dx3 * (T[p + 1] - T[p])) - (rho3 * c3 * dx3 / 2 + rho2 * c2 * dx2 / 2) / dt * (T[p] - Tp[p]));

				//mnm.setError(p, 0.0, T[p] - Tout);
				p++;

				for (p = X1+X2+1; p <= Xtot-X4-1; p++) {

					///////////////1-D energy transport equation////////
					//mnm.setError(p, 0.0, (T[p + 1] - 2 * T[p] + T[p - 1]) / (dx * dx) - (1 / a) * (T[p] - Tp[p]) / dt);
					mnm.setError(p, 0.0, k3 / dx3 * T[p - 1] - (2 * k3 / dx3 + rho3 * c3 * dx3 / dt) * T[p] + k3 / dx3 * T[p + 1] + rho3 * c3 * dx3 / dt * Tp[p]);

				}
				///////////////3rd layer boundary condition////////
				mnm.setError(p, 0.0, -k3 / dx3 * (T[p] - T[p - 1]) - (-k4 / dx4 * (T[p + 1] - T[p])) - (rho3 * c3 * dx3 / 2 + rho4 * c4 * dx4 / 2) / dt * (T[p] - Tp[p]));

				//mnm.setError(p, 0.0, T[p] - Tout);
				p++;

				for (p = X1+X2+X3+1; p <= Xtot-1; p++) {

					///////////////1-D energy transport equation////////
					//mnm.setError(p, 0.0, (T[p + 1] - 2 * T[p] + T[p - 1]) / (dx * dx) - (1 / a) * (T[p] - Tp[p]) / dt);
					mnm.setError(p, 0.0, k4 / dx4 * T[p - 1] - (2 * k4 / dx4 + rho4 * c4 * dx4 / dt) * T[p] + k4 / dx4 * T[p + 1] + rho4 * c4 * dx4 / dt * Tp[p]);

				}


				///////////////room-side boundary condition////////
				//mnm.setError(p, 0.0, k * (T[p] - T[p - 1]) / (dx)-hroom * (Troom - T[p]));
				mnm.setError(p, 0.0, hroom * (T[p] - Troom) + k4 * (T[p] - T[p - 1]) / (dx4)+rho4 * c4 * dx4 / 2.0 * (T[p] - Tp[p]) / dt);

				p++;

				//			mnm.setError( p ,  );
				//		p++;


						// mnm.prt();				//エラー表示
				mnm.prt_sum();			//エラーの合計を表示

			}
		}



		qroom = -k4 * (T[Xtot] - T[Xtot - 1]) / (dx4);
		fout << qroom << ",";
		Qroom = -k4 * As * (T[Xtot] - T[Xtot - 1]) / (dx4);
		fout << Qroom << ",";

		qout = -k1 * (T[1] - T[0]) / (dx1);
		fout << qout << ",";
		Qout = -k1 * As * (T[1] - T[0]) / (dx1);
		fout << Qout << ",";
		QQout = hout * As * (Tout - T[0]);
		fout << QQout << ",";
		QQroom = hroom * As * (T[Xtot] - Troom);
		fout << QQroom << ",";
		fout << Troom << ",";
		///////////////assign next-time-step conditions////////
		for (p = 0; p <= Xtot; p++) {

			Tp[p] = T[p];
			fout << T[p] << ",";

		}

		fout << 0 << endl;



	}




}

*/
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "CNewtonRaphsonMethod.h"


int main()
{
	int p, i, j, N, M, X, n;
	N = 8640;///// calculation times. N*dt= total time
	M = 1;/////
	X = 200;

	//double	A, B, D, Dh, d, st, slp, Sgla, Sgls, sgla, sgls, ug, ugf, ul, ulf, ugla, ugls, Ug, Ul, Up, Ugla, Ugls, tet, L, G, rhog, rhol, rhop, myug, myul;
	//double  Mi, r, rg, rl, Rel, Reg, Re, g, CA,  X, Xi, P, X0, f, E, W, WD, We, Weg, Wel, beta, Alf, Alff, PI, m, eps;

	double  QQout, QQroom, qout, Qout, qroom, Qroom, L, dx, dt, t, Ti, Troom, Tset, Tout;
	double  As, a, k, V, rho, hout, hroom, c;
	double  Tp[400], T[400];
	double  dT, Cair, Croom, Qac,Qin, kroom, Vair, cair, rhoair;

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
	dt = 10.0;				///////time step [s]
	Ti = 35.0;				///////initial wall temperature [oC]
	Tout = 35.0;			///////external ambient temperature [oC]
	Tset = 27.0;			///////room temperature [oC]
	Troom = Tset;           ///////room temperature [oC] initialization
	L = 0.20;				///////wall thickness [m]
	dx = L / double(X);		///////mesh size [m]

	///////////////wall properties////////////////////////// 

	V = 10;				///////volume [m3]
	//k = 1.7;				///////thermal conductivity [J/m-K]
	//rho = 2200.0;			///////density [kg/m3]
	//c = 840.0;				///////specific heat [J/kg-K]
	k = 1.2;				///////thermal conductivity [J/m-K]
	rho = 2310;			    ///////density [kg/m3]
	c = 879;				///////specific heat [J/kg-K]
	hout = 5.0;				///////external heat transfer coefficeint [W/m2-K]
	hroom = 5.0;			///////internal heat transfer coefficeint [W/m2-K]
	a = k / rho / c;		///////thermal diffusivity [m2/s]
	As = 45;				///////heat transfer area [m2]

	///////////////air properties//////////////////////////

	dT = 0;                 ///////Temperature change [oC] initialization
	Vair = 27;              ///////volume of room air[m3]
	rhoair = 1.2;           ///////density [kg/m3]
	cair = 1005;            ///////specific heat [J/kg-K]
	Cair = cair * Vair * rhoair;        ///////Heat capacity of room's air [J/K]
	kroom = 1000;                ////////coefficient of Croom divide Cair      
	Croom = kroom * Cair;      ///////Heat capacity of whole indoor space [J/K]

	fout << t << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";
	fout << 0 << ",";

	///////////////temperature profile initialization////////
	for (p = 0; p <= X; p++) {

		//T[p] = Tout;
		//Tp[p] = Tout;
		//T[p] = 32.0 - 2.0 / L * p * dx;
		//Tp[p] = 32.0 - 2.0 / L * p * dx;
		//T[p] = 32.178 - 11.765 * p * dx;           //this is for linear T distribution
		//Tp[p] = 32.178 - 11.765 * p * dx;          //this is for linear T distribution

		T[p] = 33.457 - 10.816 * p * dx;								//this is for another linear T distribution
		Tp[p] = 33.457 - 10.816 * p * dx;								//this is for another linear T distribution
		//T[p] = 27;									//this is for uniform T distribution
		//Tp[p] = 27;									//this is for uniform T distribution
		fout << T[p] << ",";

	}

	fout << 0 << endl;


	for (i = 1; i <= N; i++) {///////////////

		
		t = t + dt;
		fout << t << ",";
		n = i * dt / 1800;
		//if (n % 2 == 0)
			//Troom = Tset;

		//else Troom = 20.0;
		//Troom = Tset + 1.0 * sin(2.0 * 3.14 * dt * double(i) / 4500.0);



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
				mnm.setError(p, 0.0, hout * (Tout - T[p]) + k * (T[p + 1] - T[p]) / dx - rho * c * dx / 2.0 * (T[p] - Tp[p]) / dt);

				//mnm.setError(p, 0.0, T[p] - Tout);
				p++;

				for (p = 1; p <= X - 1; p++) {

					///////////////1-D energy transport equation////////
					//mnm.setError(p, 0.0, (T[p + 1] - 2 * T[p] + T[p - 1]) / (dx * dx) - (1 / a) * (T[p] - Tp[p]) / dt);
					mnm.setError(p, 0.0, k / dx * T[p - 1] - (2 * k / dx + rho * c * dx / dt) * T[p] + k / dx * T[p + 1] + rho * c * dx / dt * Tp[p]);

				}

				///////////////room-side boundary condition////////
				//mnm.setError(p, 0.0, k * (T[p] - T[p - 1]) / (dx)-hroom * (Troom - T[p]));
				mnm.setError(p, 0.0, hroom * (T[p] - Troom) + k * (T[p] - T[p - 1]) / (dx)+rho * c * dx / 2.0 * (T[p] - Tp[p]) / dt);
				
				p++;

				//			mnm.setError( p ,  );
				//		p++;


						// mnm.prt();				//エラー表示
				mnm.prt_sum();			//エラーの合計を表示

			}
		}



		qroom = -k * (T[X] - T[X - 1]) / (dx);        //conduction heat per area
		fout << qroom << ",";
		Qroom = -k * As * (T[X] - T[X - 1]) / (dx);        //conduction heat totally 
		fout << Qroom << ",";

		qout = -k * (T[1] - T[0]) / (dx);          //conduction
		fout << qout << ",";
		Qout = -k * As * (T[1] - T[0]) / (dx);        //convection heat totally
		fout << Qout << ",";
		QQout = hout * As * (Tout - T[0]);
		fout << QQout << ",";
		QQroom = hroom * As * (T[X] - Troom);
		fout << QQroom << ",";

		Qac = 0;										//air-conditioner heat flux (power)
		Qin = Qroom + QQroom + Qac;						//heat flux to room in total
		dT = Qin * dt / Croom;							//temperature change
		Troom = Troom + dT;
		fout << Qin << ",";
		fout << Troom << ",";
		///////////////assign next-time-step conditions////////
		for (p = 0; p <= X; p++) {

			Tp[p] = T[p];
			fout << T[p] << ",";

		}

		fout << 0 << endl;



	}




}

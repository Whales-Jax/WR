#include <iostream>
#include <fstream>
#include <cmath>
#include "multi_newton.h"
#include "Property.h"
using namespace std;

	
//****************************************************************************
//void関数　各装置での計算式を設定
//****************************************************************************


//*********************
//①蒸発器（EVAPORATER)
//*********************


//****
//引数
//****

//入口・・・流量、圧力、比エンタルピ
//出口・・・流量、圧力、比エンタルピ（ポインタ）
//充填量M


 void Eva(double Eva_Rin_G, double Eva_Rin_P, double Eva_Rin_h,
		 double *Eva_Rout_G, double *Eva_Rout_P, double *Eva_Rout_h, double *M){


//*******
//変数設定
//*******

	//熱交換量
	double Eva_Ref_Q;

	//仮の比エンタルピ
	double Eva_Ref_h_if;  
	
	//出口における冷媒充填量	
	double Eva_Rout_M = 0;  


	//配列設定　※分割数20の分布定数系を用いるため、配列の要素の数を20とする.
	//流量
	double Eva_Ref_G[21];
	
	//比エンタルピ
	double Eva_Ref_h[21];
	
	//温度
	double Eva_Ref_T[21];
	
	//密度
	double Eva_Ref_p[21];

	//蒸発器入口空気温度（12℃=285K)
	double Eva_Air_T = 12.0;


//******
//設計値
//******

	//時間
//	double Eva_ele_t = 0.0;

	//円周率π
	double Eva_ele_pi = atan(1.0) * 4.0;
	
	//伝熱管長さ
	double Eva_ele_L = 1.0;
	
	//伝熱管を20分割した時の長さ
	double Eva_ele_dx = Eva_ele_L/20;

	//伝熱管直径
	double Eva_ele_D = 0.01;

	//蒸発器熱通過率
	double Eva_ele_K = 13.2;

	//伝熱管断面積
	double Eva_ele_a = Eva_ele_pi * Eva_ele_D * Eva_ele_D / 4;	


//******	
//物性値
//******

	//Fluid関数名：Eva_r
	Fluid Eva_r;

	//R134Aの物性値呼び出し
	Property ref("refprop00.DLL");  
	ref.setup("R134A.FLD","IIR");  


//******
//計算式
//******
	
	//蒸発器入口
	//0番目の比エンタルピにvoid関数の引数Eva_Rin_hを代入.
	Eva_Ref_h[0] = Eva_Rin_h;

	//0番目の流量にvoid関数の引数Eva_Rin_Gを代入.
	Eva_Ref_G[0] = Eva_Rin_G;
	

	//多変数ニュートン法
	for(int i=0 ; i<20 ; i++){
				multi_newton mnm;
				mnm.setup(10);

				//初期値の設定  ※仮の比エンタルピに代入する事となる.
				mnm.set_assumed_value( 0 , 300 );

				mnm.initial();

				for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){
					for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){

					//加速度勾配
					//mnm.acc = 0.8;

					//加速度勾配自動設定（収束しやすいが遅い）
					//mnm.spcacc();
					//変数に初期値を代入

					//仮の比エンタルピに初期値を代入.
					Eva_Ref_h_if = mnm.set_value[0];

					//流量一定の式
					//(i番目の流量)＝(i+1番目の流量)
					Eva_Ref_G[i+1] = Eva_Ref_G[i];

					//Eva_r関数の呼び出し
					//入力：入口圧力、i+1番目の仮の比エンタルピ
					ref.state_ph(Eva_Rin_P, Eva_Ref_h_if, Eva_r);
					
					//Eva_r関数の出力結果を、i+1番目の温度と密度とする.
					Eva_Ref_T[i+1] = Eva_r.T;
					Eva_Ref_p[i+1] = Eva_r.rho;

					//i+1番目の場所における、空気と冷媒の熱交換量
					Eva_Ref_Q = Eva_ele_K * Eva_ele_pi * Eva_ele_D * 
								Eva_ele_dx*(Eva_Ref_T[i+1] - Eva_Air_T); //<0

					//i+1番目の比エンタルピ
					Eva_Ref_h[i+1]=(Eva_Ref_G[i]*Eva_Ref_h[i]-Eva_Ref_Q)/Eva_Ref_G[i+1];

					//エラー値
					//仮の比エンタルピとi+1番目のエンタルピを比較
					mnm.set_error2( 0 , Eva_Ref_h[i+1], Eva_Ref_h_if );
					
					//エラー表示
					//mnm.prt();

					//エラーの合計を表示
					//mnm.prt_sum();
					}
				}

				//蒸発器出口における冷媒充填量
				Eva_Rout_M = Eva_Rout_M + Eva_Ref_p[i+1] * Eva_ele_pi * Eva_ele_D * Eva_ele_D * Eva_ele_dx/4;
	}

	Eva_Rin_G = Eva_Ref_G[0];
	Eva_Rin_h = Eva_Ref_h[0];

	//ポインタ
	//*Eva_Rout_GにEva_Ref_G[20]を収納.
	*Eva_Rout_G = Eva_Ref_G[20];

	//*Eva_Rout_PにEva_Rin_Pを収納.
	*Eva_Rout_P = Eva_Rin_P;

	//*Eva_Rout_hにEva_Ref_h[20]を収納.
	*Eva_Rout_h = Eva_Ref_h[20];

	//*MにEva_Rout_Mを収納.
	*M = Eva_Rout_M;

	//蒸発器出口のGとｈから　同じ出口におけるT,P,pを呼び出す.
	cout << "Eva has caliculated" << endl;

	ref.state_ph(Eva_Rin_P , Eva_Rin_h , Eva_r);
	cout << "入口T = " << Eva_r.T<<endl;

//出口P,hからTを求める.
	ref.state_ph(*Eva_Rout_P , *Eva_Rout_h , Eva_r);
	cout << "出口T = " << Eva_r.T << " " << "出口P = " << Eva_r.P << " " << "出口rho = " << Eva_r.rho << endl;
	cout << "出口h = " << Eva_r.h << " " << "出口s = " << Eva_r.s << endl;
	cout << "出口G = " << *Eva_Rout_G << " " << "出口M = " << Eva_Rout_M << endl;
}


//*********************
//②圧縮機（COMPRESSOR)
//*********************

void Com(double Com_Rin_G, double Com_Rin_P, double Com_Rin_h,
		 double *Com_Rout_h , double *Com_Rout_P, double *Com_Rout_G, double *Com_Rout_M ){

//*******
//変数設定
//*******

	//入口 比エントロピー
	double Com_Rin_s ;
	
	//出口 比エントロピー
	double Com_Rout_s ;
	
	//出口 比エンタルピ（断熱過程）
	double Com_Rout_had ;

	//流量一定
	*Com_Rout_G = Com_Rin_G ;
	
	//設計値
	double Com_ele_W = 0.45;
	double Com_ele_V = 0.00002;
	double Com_ele_eff = 1.0;
	Fluid Com_r;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");
	
	multi_newton mnm;
	mnm.setup(10);

	//出口圧力を仮に設定
	mnm.set_assumed_value( 0 , 1000 );

	mnm.initial();
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){

		//加速度勾配
		//mnm.acc = 0.8;

		//加速度勾配自動設定（収束しやすいが遅い）
		//mnm.spcacc();

		//*出口圧力に初期値代入
		*Com_Rout_P = mnm.set_value[0];
		
		//Com_r関数の呼び出し
		//Pとhからsを求める.(at入口）
		ref.state_ph(Com_Rin_P , Com_Rin_h , Com_r);
		Com_Rin_s = Com_r.s;

		//*出口エンタルピ
		*Com_Rout_h = (Com_Rin_G * Com_Rin_h + Com_ele_W)/(*Com_Rout_G);
		
		//断熱効率を考慮した出口エンタルピ
		Com_Rout_had = Com_ele_eff * (*Com_Rout_h - Com_Rin_h) +Com_Rin_h;

		//Pとhadからsを求める.(at出口)
		ref.state_ph(*Com_Rout_P , Com_Rout_had , Com_r);
		Com_Rout_s = Com_r.s;

		//エラー値
		//入口と出口のエントロピー比較.
		mnm.set_error2(0 , Com_Rin_s , Com_Rout_s );
		}
	} 

	//*Rout_Pと*Rout_hからT,P,rho,h,sを求める.
	ref.state_ph(*Com_Rout_P , *Com_Rout_h , Com_r);
	cout << "Com has caliculated（全て出口の値）" << endl;
	cout << "T = " << Com_r.T << " " << "P = " << Com_r.P << " " << "rho = " << Com_r.rho << endl;
	cout << "h = " << Com_r.h << " " << "s = " << Com_r.s << endl;

	//*出口冷媒充填量
	*Com_Rout_M = Com_ele_V * Com_r.rho;
	cout << "G = " << *Com_Rout_G << " " << "M = " << Com_ele_V*Com_r.rho << endl;
}



//*********************
//③凝縮器（CONDENSER)
//*********************

void Con(double Con_Rin_h, double Con_Rin_P, double Con_Rin_G,
		 double *Con_Rout_h, double *Con_Rout_P, double *Con_Rout_G, double *M1){

	double Con_Ref_Q;
	double Con_Ref_h_if;
	double Con_Rout_M = 0;

	double Con_Ref_G[21];
	double Con_Ref_h[21];
	double Con_Ref_T[21];
	double Con_Ref_p[21];

	//凝縮器入口空気温度(32.0℃=305.0K)
	double Con_Air_T = 32.0;

	//設計値
//	double Con_ele_t = 0.0;
	double Con_ele_pi = atan(1.0) * 4.0;
	double Con_ele_L = 1.0;
	double Con_ele_dx = Con_ele_L/20;
	double Con_ele_D = 0.01;
	double Con_ele_K = 13.3;
	double Con_ele_a = Con_ele_pi * Con_ele_D * Con_ele_D / 4;
	
	//物性値
	Fluid Con_r;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");

	//h,Gの設定at凝縮器入口
	Con_Ref_h[0] = Con_Rin_h;
	Con_Ref_G[0] = Con_Rin_G;

	Con_Rout_M = 0;

			for(int i=0 ; i<20 ; i++){
				multi_newton mnm;
				mnm.setup(10);
				//初期値の設定
				mnm.set_assumed_value( 0 , 300 );

				mnm.initial();

				for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){
					for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){

					//加速度勾配
					//mnm.acc = 0.8;

					//加速度勾配自動設定（収束しやすいが遅い）
					//mnm.spcacc();

					//変数に初期値を代入
					//仮のエンタルピに初期値代入.
					Con_Ref_h_if = mnm.set_value[0];
					
					//流量一定
					Con_Ref_G[i+1] = Con_Ref_G[i];
					
					//入口Pと仮エンタルピからi+1番目のTとpを求める.
					ref.state_ph(Con_Rin_P, Con_Ref_h_if, Con_r);
					Con_Ref_T[i+1] = Con_r.T;
					Con_Ref_p[i+1] = Con_r.rho;

					//交換熱量
					Con_Ref_Q = Con_ele_K * Con_ele_pi * Con_ele_D * 
								Con_ele_dx*(Con_Ref_T[i+1] - Con_Air_T);	//<0

					//冷媒エンタルピat(i+1)番目
					Con_Ref_h[i+1]=(Con_Ref_G[i]*Con_Ref_h[i]-Con_Ref_Q)/Con_Ref_G[i+1];

					//エラー値
					//i+1番目エンタルピと仮エンタルピ比較
					mnm.set_error2( 0 , Con_Ref_h[i+1], Con_Ref_h_if );
					
					//エラー表示
					//mnm.prt();

					//エラーの合計を表示
					//mnm.prt_sum();
					}
				}
				
				//出口冷媒充填量
				Con_Rout_M = Con_Rout_M + Con_Ref_p[i+1] * Con_ele_pi * Con_ele_D * 
								Con_ele_D * Con_ele_dx/4;
			}

	Con_Rin_G = Con_Ref_G[0];
	Con_Rin_h = Con_Ref_h[0];

	//出口における状態量に20番目の値をそれぞれ代入.
	*Con_Rout_G = Con_Ref_G[20];
	*Con_Rout_P = Con_Rin_P;
	*Con_Rout_h = Con_Ref_h[20];
	*M1 = Con_Rout_M;

	//出口PとhからT,rho,s,G,Mを求める.


	//計算結果出力
	cout << "Con has caliculated" << endl;
	ref.state_ph(Con_Rin_P , Con_Rin_h , Con_r);
	cout << "入口T = " << Con_r.T << endl;
	
	ref.state_ph(*Con_Rout_P , *Con_Rout_h , Con_r);
	cout << "出口T = " << Con_r.T << " " << "出口P = " << Con_r.P << " " << "出口rho = " << Con_r.rho << endl;
	cout << "出口h = " << Con_r.h << " " << "出口s = " << Con_r.s << endl;
	cout << "出口G = " << *Con_Rout_G << " " << "出口M = " << Con_Rout_M << endl;	
}


//*******************************
//④アキュムレーター（ACCUMULATOR)
//*******************************

void Acc(double Acc_Vin_P, double Acc_Vin_G, double Acc_Lin_L, double *Acc_Vin_h,
		 double *Acc_Vout_G, double *Acc_Vout_P, double *Acc_Vout_h, double *Acc_Lout_M, double *Acc_Vout_M){

	//入口
//	double Acc_Vin_G_int;
//	double Acc_Vin_P_int;
	double Acc_Vin_h_int;

	double Acc_Lout_rhol = 0;
	double Acc_Vout_rhov = 0;

	//設計値
	double Acc_ele_a = 0.01 ;
	double Acc_ele_L = 1.0;
	double Acc_ele_g = 9.80665 ;

	Fluid Acc_rv;
	Fluid Acc_rl;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");


	ref.sat_p(Acc_Vin_P , Acc_rl , Acc_rv);

	//入口Pから入口hを算出.
	*Acc_Vin_h = Acc_rv.h;

	//入口蒸気比エンタルピ
	Acc_Vin_h_int = *Acc_Vin_h;

	//出口h=入口h
	*Acc_Vout_h = Acc_Vin_h_int;

	//流量一定
	*Acc_Vout_G = Acc_Vin_G;

	ref.sat_p(Acc_Vin_P , Acc_rl , Acc_rv);
	Acc_Vout_rhov = Acc_rv.rho;
	Acc_Lout_rhol = Acc_rl.rho;

//	*Acc_Lout_P = Acc_Vin_P_int + (Acc_Vout_rhol) / 1000 * Acc_ele_g * Acc_Lin_L;
	*Acc_Lout_M = Acc_ele_a * Acc_Lin_L * Acc_Lout_rhol;
	*Acc_Vout_M = Acc_ele_a * (Acc_ele_L -  Acc_Lin_L) * Acc_Vout_rhov;
	*Acc_Vout_P = Acc_Vin_P;

	cout << "Acc has caliculated（全て出口の値）" << endl;
	cout << "G = " << *Acc_Vout_G << " " << "P = " << *Acc_Vout_P << endl;
	cout << "h = " << *Acc_Vin_h << " " << "M = " << *Acc_Lout_M + *Acc_Vout_M << endl;
}


//**************************
//⑤膨張弁（EXPANSION VALVE)
//**************************

void Exp(double Exp_Rin_h, double Exp_Rin_P, double Exp_Rin_G,
		 double *Exp_Rout_h, double *Exp_Rout_P, double *Exp_Rout_G){
	
	double Exp_Rin_p ;

	//設計値
	double Exp_ele_Cv = 0.6;
	double Exp_ele_a = 0.0000269;     //2.69 * 1e-5;

	Fluid Exp_r;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");

	//計算
	ref.state_ph(Exp_Rin_P , Exp_Rin_h , Exp_r);
	Exp_Rin_p = Exp_r.rho;

	//流量一定
	*Exp_Rout_G = Exp_Rin_G;
	
	*Exp_Rout_h = Exp_Rin_h;
	*Exp_Rout_P = Exp_Rin_P-(Exp_Rin_G*Exp_Rin_G)/(2 * Exp_Rin_p * Exp_ele_Cv * Exp_ele_Cv * Exp_ele_a * Exp_ele_a); 

	ref.state_ph(*Exp_Rout_P , *Exp_Rout_h , Exp_r);
	cout << "Exp has caliculated（出口の値）" << endl;
	cout << "T = " << Exp_r.T << " " << "P = " << Exp_r.P << " " << "rho = " << Exp_r.rho << endl;
	cout << "h = " << Exp_r.h << " " << "s = " << Exp_r.s << endl;
	cout << "G = " << *Exp_Rout_G << endl;
}


//****************************************************************************
//main関数　void関数で作った各装置の計算式を実行.（ここで初期値を入力値として使用）
//****************************************************************************

int main(){

	Fluid Eva_r;
	Fluid Con_r;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");

	double Eva_Rin_G, Eva_Rin_P, Eva_Rin_h ;
	double Com_Rin_G, Com_Rin_P, Com_Rin_h ;
	double Con_Rin_G, Con_Rin_P, Con_Rin_h ;
	double Acc_Vin_G, Acc_Vin_P, Acc_Lin_L ;
	double Exp_Rin_G, Exp_Rin_P, Exp_Rin_h ;
	

	double ALL_M;
	double Eva_Rout_G = 0;
	double Eva_Rout_P = 0;
	double Eva_Rout_h = 0;
	double Eva_Rout_M = 0;
	double Com_Rout_G = 0;
	double Com_Rout_P = 0;
	double Com_Rout_h = 0;
	double Com_Rout_M = 0;
	double Con_Rout_G = 0;
	double Con_Rout_P = 0;
	double Con_Rout_h = 0;
	double Con_Rout_M = 0;
	double Acc_Vout_G = 0;
	double Acc_Vout_P = 0;
	double Acc_Vout_h = 0;
	double Acc_Vin_h = 0;
	double Acc_Lout_M = 0;
	double Acc_Vout_M = 0;
	double Exp_Rout_G = 0;
	double Exp_Rout_P = 0;
	double Exp_Rout_h = 0;
				
			cout << "SYSTEM 静特性" << endl;

			//静特性
			multi_newton mnm;
			mnm.setup(20);
			//初期値の設定
			mnm.set_assumed_value( 0  , 0.02 );
			mnm.set_assumed_value( 1  , 350.0 );
			mnm.set_assumed_value( 2  , 250.0 );
			mnm.set_assumed_value( 3  , 0.02 );
			mnm.set_assumed_value( 4  , 1000.0 );
			mnm.set_assumed_value( 5  , 425.0 );
			mnm.set_assumed_value( 6  , 0.02 );
			mnm.set_assumed_value( 7  , 1000.0 );
			mnm.set_assumed_value( 8  , 250.0 );
			mnm.set_assumed_value( 9  , 0.02 );
			mnm.set_assumed_value( 10 , 350.0 );
			mnm.set_assumed_value( 11 , 400.0 );
			mnm.set_assumed_value( 12 , 0.02 );
			mnm.set_assumed_value( 13 , 350.0 );
			mnm.set_assumed_value( 14 , 0.1 );
			mnm.initial();

			for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){
				for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){

				//加速度勾配
				//mnm2.acc = 0.8;

				//加速度勾配自動設定（収束しやすいが遅い）
				//mnm2.spcacc();

				//変数に初期値を代入
				Eva_Rin_G = mnm.set_value[0];
				Eva_Rin_P = mnm.set_value[1];
				Eva_Rin_h = mnm.set_value[2];
				Con_Rin_G = mnm.set_value[3];
				Con_Rin_P = mnm.set_value[4];
				Con_Rin_h = mnm.set_value[5];
				Exp_Rin_G = mnm.set_value[6];
				Exp_Rin_P = mnm.set_value[7];
				Exp_Rin_h = mnm.set_value[8];
				Com_Rin_G = mnm.set_value[9];
				Com_Rin_P = mnm.set_value[10];
				Com_Rin_h = mnm.set_value[11];
				Acc_Vin_G = mnm.set_value[12];
				Acc_Vin_P = mnm.set_value[13];
				Acc_Lin_L = mnm.set_value[14];


				//計算　＆Xで変数Xのアドレスを示し、&Xを読み込む事で変数Xを呼び出せる.

				Eva(Eva_Rin_G, Eva_Rin_P, Eva_Rin_h, &Eva_Rout_G, &Eva_Rout_P, 
						&Eva_Rout_h, &Eva_Rout_M);
				Acc(Acc_Vin_P, Acc_Vin_G, Acc_Lin_L, &Acc_Vin_h, &Acc_Vout_G, 
						&Acc_Vout_P, &Acc_Vout_h, &Acc_Lout_M, &Acc_Vout_M);
				Com(Com_Rin_G, Com_Rin_P, Com_Rin_h, &Com_Rout_h, &Com_Rout_P, 
						&Com_Rout_G, &Com_Rout_M);
				Con(Con_Rin_h, Con_Rin_P, Con_Rin_G, &Con_Rout_h, &Con_Rout_P, 
						&Con_Rout_G, &Con_Rout_M);
				Exp(Exp_Rin_h, Exp_Rin_P, Exp_Rin_G, &Exp_Rout_h, &Exp_Rout_P, 
						&Exp_Rout_G);
				ALL_M = Eva_Rout_M + Con_Rout_M + Com_Rout_M + Acc_Lout_M + Acc_Vout_M;

				cout << "ALL_M =0." << ALL_M << endl;
	cout << "," << endl;

				//エラー値
				mnm.set_error2( 0  , ALL_M , 1.5);
				mnm.set_error2(  1 , Con_Rout_P , Exp_Rin_P);
				mnm.set_error2(  2 , Con_Rout_h , Exp_Rin_h);
				mnm.set_error2(  3 , Exp_Rout_G , Eva_Rin_G);
				mnm.set_error2(  4 , Exp_Rout_P , Eva_Rin_P);
				mnm.set_error2(  5 , Exp_Rout_h , Eva_Rin_h);
				mnm.set_error2(  6 , Eva_Rout_G , Acc_Vin_G);
				mnm.set_error2(  7 , Eva_Rout_P , Acc_Vin_P);
				mnm.set_error2(  8 , Eva_Rout_h , Acc_Vin_h);
				mnm.set_error2(  9 , Acc_Vout_G , Com_Rin_G);
				mnm.set_error2( 10 , Acc_Vout_P , Com_Rin_P);
				mnm.set_error2( 11 , Acc_Vout_h , Com_Rin_h);
				mnm.set_error2( 12 , Com_Rout_G , Con_Rin_G);
				mnm.set_error2( 13 , Com_Rout_P , Con_Rin_P);
				mnm.set_error2( 14 , Com_Rout_h , Con_Rin_h);
				
				//エラー表示
				//mnm.prt();
				//エラーの合計を表示
				mnm.prt_sum();
				}
			}

//**************
//計算結果の出力
//**************

			cout << "\t" << "G" << "\t" << "\t" << "P" << "\t" << "\t" << "h" << "\t" << endl;

			cout << "EVA" << "\t" << Eva_Rin_G << "\t" << Eva_Rin_P << "\t" <<Eva_Rin_h << "\t" << endl;

			cout << "COM" << "\t" << Com_Rin_G << "\t" << Com_Rin_P << "\t" << Com_Rin_h << "\t" << endl;

			cout << "CON" << "\t" << Con_Rin_G << "\t" << Con_Rin_P << "\t" << Con_Rin_h << "\t" << endl;

			cout << "Acc" << "\t" << Acc_Vin_G << "\t" << Acc_Vin_P << "\t" <<Acc_Vin_h << "\t" << endl;

			cout << "EXP" << "\t" << Exp_Rin_G << "\t" << Exp_Rin_P << "\t" <<Exp_Rin_h << "\t" << "\t" << endl;
			cout << "," <<endl;
			cout << "EVA" << "\t" << Eva_Rout_G << "\t" << Eva_Rout_P << "\t" <<Eva_Rout_h << "\t" << endl;

			cout << "COM" << "\t" << Com_Rout_G << "\t" << Com_Rout_P << "\t" <<Com_Rout_h << "\t" << endl;

			cout << "CON" << "\t" << Con_Rout_G << "\t" << Con_Rout_P << "\t" <<Con_Rout_h << "\t" << endl;

			cout << "Acc" << "\t" << Acc_Vout_G << "\t" << Acc_Vout_P << "\t" <<Acc_Vout_h << "\t" << endl;

			cout << "EXP" << "\t" << Exp_Rout_G << "\t" << Exp_Rout_P << "\t" <<Exp_Rout_h << "\t" << endl;
			cout << "ALL_M" << "\t" << ALL_M << endl;
			cout << "Acc_ele.L = " << Acc_Lin_L << endl;

			ref.state_ph(Eva_Rin_P , Eva_Rin_h , Eva_r);
			ref.state_ph(Con_Rin_P , Con_Rin_h , Con_r);

			cout << "T_Eva_In = " << Eva_r.T << "\t" << "T_Con_In = " << Con_r.T <<endl;

			ref.state_ph(Eva_Rout_P , Eva_Rout_h , Eva_r);
			ref.state_ph(Con_Rout_P , Con_Rout_h , Con_r);

			cout << "T_Eva_Out = " << Eva_r.T << "\t" << "T_Con_Out = " << Con_r.T <<endl;
		
		

}
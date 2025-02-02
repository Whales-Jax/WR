//熱交換器
#include"DXProperty_ver06.h"
#include<math.h>
using namespace std;
#include "CNewtonRaphsonMethod.h"
#include"CFluidParameter.h"
void main()
{
	int n;
	double P;	// 圧力[kPa]
	double h;	// エンタルピ[kJ/kg]
	double T;	// 温度[℃]
	double x;   //乾き度
	double rho1;
	double hl2;
	double rho2;
	double A;double UA;
	double mc;
	double mh;
    double ans;
	double Ti,To,Q,cpi,cpo,ho,dt;
	DXProperty ref;					// DXPropertyクラスの宣言

	CFluidParameter H_t1;

	mc=0.3;
	mh=0.2;
    A=2.0*3.14*0.01*10.0;
    UA=A*1.0;
	// セットアップ
	ref.reibai = "R410A.PPF";		// 冷媒の設定
	ref.joutai = "IIR";				// 基準状態の設定 ( 普通はIIR, 水はNBP )
	ref.LoadDLL("refprop.DLL");		// 読み込むDLLファイルの設定 ( 2種類以上の冷媒物性が必要な場合，それぞれの冷媒に対してDLLファイルが必要 )
	ref.setup();					// セットアップ

	//初期値
	Ti =10.0;
	To = 10.0;
	Q = 0;
	
	CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
	mnm.setup(10);					//セットアップ　※引数はニュートン法の変数以上にする

	// 初期値の代入
	mnm.setValue( 1 , To );
	mnm.setValue( 2 , Q );

	//mnm.setAcc(0.8);				//加速度勾配の入力　※なくても可

	mnm.initial();					//計算を開始する前に必ず初期化してください
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// おまじない
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// おまじない

			//値の設定
			To = mnm.getValue(0);
			Q = mnm.getValue(1);
            
	// 物性値の計算
	//ref.state_ph( 1000, 350 );	// 任意状態関数( 圧力，エンタルピ ) -> 物性値はRcクラスに格納される
	ref.state_tp( 40, 4000 );	// 単相状態関数( 温度, 圧力 )		-> 物性値はRcクラスに格納される
	//ref.sat_p( P );			// 飽和状態関数( 圧力 ) -> 飽和液の物性値はRlクラスに, 飽和蒸気の物性値はRvクラスに格納される
	//ref.sat_t( T );			// 飽和状態関数( 温度 ) -> 飽和液の物性値はRlクラスに, 飽和蒸気の物性値はRvクラスに格納される

	// 物性値の呼び出し
	hl2 = ref.Rc.h;			//Rcクラスの温度を呼び出すとき
	rho1 = ref.Rc.rho;
	cout<< hl2<<endl;
	cout<< rho1<<endl;



	// 物性値の計算
	//ref.state_ph( 1000, 350 );	// 任意状態関数( 圧力，エンタルピ ) -> 物性値はRcクラスに格納される
	ref.state_tp( Ti, 4000 );	// 単相状態関数( 温度, 圧力 )		-> 物性値はRcクラスに格納される
	//ref.sat_p( P );			// 飽和状態関数( 圧力 ) -> 飽和液の物性値はRlクラスに, 飽和蒸気の物性値はRvクラスに格納される
	//ref.sat_t( T );			// 飽和状態関数( 温度 ) -> 飽和液の物性値はRlクラスに, 飽和蒸気の物性値はRvクラスに格納される
		        ho = ref.Rv.h;
	// 物性値の呼び出し
	hl2 = ref.Rc.h;			//Rcクラスの温度を呼び出すとき
	rho2 = ref.Rc.rho;
	cout<< hl2<<endl;
	cout<< rho2<<endl;

			ref.state_ph( 4000, hl2 );
			Ti= ref.Rc.T;
            ref.state_ph( 4000, ho );
			To= ref.Rc.T;
			//エラー値(ここに式を入力する)
			mnm.setError( 0 , Q/mc+hl2 , ho );
		    dt=(40-10)/(40-To);
			mnm.setError( 1 , UA*(To-10)/(log10(dt)) , Q );
			mnm.prt();				//エラー表示
			mnm.prt_sum();			//エラーの合計を表示
		}
	}

	cout << "ho = " << ho << endl;
}

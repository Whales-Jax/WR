/*#include"DXProperty_ver06.h"
#include<math.h>
using namespace std;
#include "CNewtonRaphsonMethod.h"
#include"CFluidParameter.h"
void main()
{//膨張弁→01,蒸発器→23,アキュムレータ→45,圧縮機→67,凝縮器→89	
	int i,z;
	double L,A,V,a,Cv,Q[1],W,G[10],P[10],h[10],rho[10],M[5],T[10],l,dL,s,s2,Tg[100],Tg2[100],hg[100],hg2[100],dp,dh;
	DXProperty ref;					// DXPropertyクラスの宣言
	// セットアップ
	ref.reibai = "R134A.FLD";		// 冷媒の設定
	ref.joutai = "IIR";				// 基準状態の設定 ( 普通はIIR, 水はNBP )
	ref.LoadDLL("refprop.DLL");		// 読み込むDLLファイルの設定 ( 2種類以上の冷媒物性が必要な場合，それぞれの冷媒に対してDLLファイルが必要 )
	ref.setup();
	L=0.11;
	Cv=0.6;
	a=2.69e-5;
	A=3.1415*0.01*1.0/20.0;
	V=3.1415*0.005*0.005*1.0;
	W=0.45;
	P[0]=1000.0;
	T[0]=20.0;//膨張弁温度初期値
    ref.state_tp( T[0], P[0]);
    h[0] = ref.Rc.h;
	Q[0]=-1.0;
	Q[1]=1.0;
	z=0;
	ref.state_ph( P[0], h[0]);
	rho[0] = ref.Rc.rho;
	G[0]=0.03;
		P[5]=P[0]-0.5*rho[0]*pow(G[0]/(rho[0]*Cv*a),2.0)*1.0e-3;
	dL=1.0;
	dp=1.0;
	dh=1.0;
		ofstream fout("圧縮サイクル.csv");
	if(!fout){
	cout << "ファイルをオープンできませんでした"<< endl;
	}
	else{
		cout<< "ファイルをオープンしました" << endl;
	}

	while(L>1.0||L<0||dp>1.0e-6||dL>1.0e-6||dh>1.0e-6){
//膨張弁
		fout<<endl;
		fout<< "///////////"<<z+1<<"周目///////////" << endl;
		ref.state_ph( P[0], h[0]);
		rho[0] = ref.Rc.rho;
		G[0]=rho[0]*Cv*a*pow((2.0*(P[0]-P[5])*1.0e3)/rho[0],0.5);
		P[1]=P[5];
		G[1]=G[0];
		h[1]=h[0];
		M[0]=0.0;
		ref.state_ph( P[1], h[1]);
		rho[1] = ref.Rc.rho;
		T[1]=ref.Rc.T;

//接続-膨張弁→凝縮器
		G[2]=G[1];
		h[2]=h[1];
		P[2]=P[1];
		rho[2] = rho[1];
		T[2]=T[1];
//凝縮器
	//初期値	
        P[3]=P[2];
        G[3]=G[2];
			Tg[0]=T[2];
			ref.state_tp( Tg[0], P[2]);
			hg[0] = ref.Rc.h;

		for(i=0;i<20;i++){
			CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
			mnm.setup(10);					//セットアップ　※引数はニュートン法の変数以上にする
// 初期値の代入
			mnm.setValue( 0 , hg[i+1] );
			mnm.setValue( 1 , Q[0] );
			mnm.setAcc(0.8);
			mnm.initial();					//計算を開始する前に必ず初期化してください
			for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// おまじない
				for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// おまじない
					hg[i+1] = mnm.getValue(0);
					Q[0] = mnm.getValue(1);
				ref.state_ph( P[3],hg[i+1]);
			Tg[i]=ref.Rc.T;				mnm.setError( 0 , Q[0]/G[2]+hg[i], hg[i+1] );
					mnm.setError( 1 , 13.3*A*(32.0-ref.Rc.T) , Q[0] );
					mnm.prt();				//エラー表示
					mnm.prt_sum();			//エラーの合計を表示
				}
			}
			hg[i+2]=hg[i+1];
			fout<<Q[0]<<",";
		}
		fout<<"Q"<<endl;
		for(i=0;i<20;i++){
			fout<<Tg[i]<<",";
		}
		fout<<"T"<<endl;

		for(i=0;i<20;i++){
			fout<<hg[i]<<",";
		}
		fout<<"h"<<endl;

		h[3]=hg[19];
		ref.state_ph( P[3],h[3]);
		T[3]=ref.Rc.T;
		rho[3] = ref.Rc.rho;
		M[1]=rho[3]*V;
        T[3]=ref.Rc.T;
		G[3]=G[2];
//接続-凝縮器→アキュムレータ
		G[4]=G[3];
		h[4]=h[3];
		P[4]=P[3];
		rho[4] = rho[3];
        T[4]=T[3];
//あきゅむれーた
		G[5]=G[4];
		ref.sat_t( T[4]);
		rho[5]=ref.Rl.rho;//中身は飽和液体です
		P[5]=ref.Rv.P;
		dp=abs(P[5]-P[4]);
		M[2]=0.01*L*rho[5];
		T[5]=T[4];
		ref.state_tp( T[5],P[5]);
		h[5]=ref.Rv.h;
//接続-アキュムレータ→圧縮機
		G[6]=G[5];
		h[6]=h[5];
		P[6]=P[5];
        T[6]=T[5];
//圧縮機
		if(z==0){
			P[7]=P[6]*2.2;
		}
		G[7]=G[6];
		ref.state_ph( P[6],h[6]);
		rho[6] = ref.Rv.rho;
		s=ref.Rc.s;
		h[7]=(G[6]*h[6]+W)/G[7];
		CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
		mnm.setup(10);					//セットアップ　※引数はニュートン法の変数以上にする
	// 初期値の代入
		mnm.setValue( 0 , s );
		mnm.setValue( 1 , P[7] );
		mnm.setAcc(0.5);
		mnm.initial();					//計算を開始する前に必ず初期化してください
		for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// おまじない
			for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// おまじない
				ref.state_ph( P[7],h[7]);
				s = mnm.getValue(0);
				P[7] = mnm.getValue(1);
				mnm.setError( 0 , ref.Rc.s , s );
				mnm.setError( 1 , ref.Rc.P , P[7] );
				mnm.prt();				//エラー表示
				mnm.prt_sum();			//エラーの合計を表示
			}
		}

			ref.state_ph( P[7], h[7]);
			T[7]=ref.Rc.T;
			rho[7]=ref.Rc.rho;
			M[3]=20.0e-6*rho[7];
//接続-圧縮器→蒸発器
		G[8]=G[7];
		h[8]=h[7];
		P[8]=P[7];
		G[9]=G[8];
		rho[8] = rho[7];
		T[8]=T[7];
//蒸発器
		P[9]=P[8];
			Tg2[0]=T[8];
			ref.state_tp( Tg2[0], P[8]);
			hg2[0] = ref.Rc.h;
		for(i=0;i<20;i++){

			CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
			mnm.setup(10);					//セットアップ　※引数はニュートン法の変数以上にする
// 初期値の代入
			mnm.setValue( 0 , hg2[i+1] );
			mnm.setValue( 1 , Q[1] );
			mnm.setAcc(0.8);
			mnm.initial();					//計算を開始する前に必ず初期化してください

			for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// おまじない
				for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// おまじない
					hg2[i+1] = mnm.getValue(0);
					Q[1] = mnm.getValue(1);
			ref.state_ph( P[8],hg2[i+1]);
			Tg2[i]=ref.Rc.T;
					mnm.setError( 0 , Q[1]/G[2]+hg2[i], hg2[i+1] );
					mnm.setError( 1 , 13.2*A*(12.0-ref.Rc.T) , Q[1] );
					mnm.prt();				//エラー表示
					mnm.prt_sum();			//エラーの合計を表示
				}
			}
			hg2[i+2]=hg2[i+1];
			fout<<Q[1]<<",";
		}
		fout<<"Q"<<endl;
		for(i=0;i<20;i++){
			fout<<Tg2[i]<<",";
		}
		fout<<"T"<<endl;

		for(i=0;i<20;i++){
			fout<<hg2[i]<<",";
		}
		fout<<"h"<<endl;

		h[9]=hg2[19];
		ref.state_ph( P[9],h[9]);
		T[9]=ref.Rc.T;
		rho[9] = ref.Rc.rho;
		M[4]=rho[9]*V;
		M[2]=1.5-M[1]-M[3]-M[4];
		l=M[2]/0.01/rho[5];
		T[9]=ref.Rc.T;
		G[0]=G[9];
		h[0]=h[9];
		P[0]=P[9];
		T[0]=T[9];
		rho[0]=rho[9];
		dL=abs(L-l);
		dh=abs(h[5]-h[4]);
		L=l;
for(i=0;i<10;i++){
	cout<<"rho["<<i<<"]"<<rho[i];
	cout<<"P["<<i<<"]"<<P[i];
	cout<<"G["<<i<<"]"<<G[i];
	cout<<"T["<<i<<"]"<<T[i];
	cout<<"h["<<i<<"]"<<h[i]<<endl;
	fout<<"rho["<<i<<"]"<<rho[i];
	fout<<"P["<<i<<"]"<<P[i];
	fout<<"G["<<i<<"]"<<G[i];
	fout<<"T["<<i<<"]"<<T[i];
	fout<<"h["<<i<<"]"<<h[i]<<endl;
	}

	cout<<"L"<<l<<endl;
	cout<<"膨張弁充てん量"<<M[0]<<endl;
	cout<<"凝縮器充てん量"<<M[1]<<endl;
	cout<<"アキュムレータ充てん量"<<M[2]<<endl;
	cout<<"圧縮機充てん量"<<M[3]<<endl;
	cout<<"蒸発器充てん量"<<M[4]<<endl;
	cout<<"冷媒流量"<<G[0]<<endl;
	cout<<"凝縮圧力"<<P[3]<<endl;
	cout<<"蒸発圧力"<<P[9]<<endl;
	cout<<"蒸発器入口エンタルピ"<<h[9]<<endl;
	cout<<"凝縮器入口温度"<<T[2]<<endl;
	cout<<"凝縮器出口温度"<<T[3]<<endl;
	cout<<"蒸発器入口温度"<<T[8]<<endl;
	cout<<"蒸発器出口温度"<<T[9]<<endl;
	cout<<"dp"<<dp<<endl;
	cout<<"dL"<<dL<<endl;
	cout<<"dh"<<dh<<endl;
	cout<<"z"<<z<<endl;
	fout<<"L"<<l<<endl;
	fout<<"膨張弁充てん量"<<M[0]<<endl;
	fout<<"凝縮器充てん量"<<M[1]<<endl;
	fout<<"アキュムレータ充てん量"<<M[2]<<endl;
	fout<<"圧縮機充てん量"<<M[3]<<endl;
	fout<<"蒸発器充てん量"<<M[4]<<endl;
	fout<<"冷媒流量"<<G[0]<<endl;
	fout<<"凝縮圧力"<<P[3]<<endl;
	fout<<"蒸発圧力"<<P[9]<<endl;
	fout<<"蒸発器入口エンタルピ"<<h[9]<<endl;
	fout<<"凝縮器入口温度"<<T[2]<<endl;
	fout<<"凝縮器出口温度"<<T[3]<<endl;
	fout<<"蒸発器入口温度"<<T[8]<<endl;
	fout<<"蒸発器出口温度"<<T[9]<<endl;
	fout<<"dp"<<dp<<endl;
	fout<<"dL"<<dL<<endl;
	fout<<"dh"<<dh<<endl;
	z++;
	for(i=0;i<=430;i++){
			ref.state_ph( 1749.5338139416667,430-i);
		T[9]=ref.Rc.T;
		cout<<T[9]<<endl;
	}

	}	
}*/
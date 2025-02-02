#include <iostream>
#include <fstream>
#include<math.h>
using namespace std;
#include "CNewtonRaphsonMethod.h"
#define M 10
#define N 10

int main(){

				ofstream fout("ヌッセルト6.csv");
	if(!fout){
	cout << "ファイルをオープンできませんでした"<< endl;
    return 1;
	}
	else
		cout<< "ファイルをオープンしました" << endl;
	



	//濃度60%	温度℃
	double u[M][N];
	double T[M][N];
	double Tv=96.44432489;		//液膜蒸気側温度
	double Ti=83.35 ;		//伝熱管温度
	double Ta=99.90 ;		//溶液入口温度
	double k=0.479938777;	//適当，ほんとは濃度，温度によって変わる
	double p=1711.790009 ;		//密度
	double m=0.3 ;
	double q=0.00650693170434067;	//粘性係数
	double cp=2069.194255;			//比熱
	double sig=0.090283664901;		//表面張力
	double c;				//μ/ρ
	c=q/p;
	double a;				//k/ρc
	a=k/(p*cp);
	int i,j;
	double g=9.80665;
	double r=1.8e-2;	//管半径
	double rr[N];		//film thickness[m]
	double dr[N];		//film thickness per a cell[m]
	double uave[N];		
	double sum=0;
	double h[N];
	double b[N];
	double bave;		//傾き平均
	double have;		//熱伝達率
	double rrave;		//液膜厚さ平均（全角度）
	double uuave;		//流速平均（全セル）
	double width[N];

	//contact angle 29.7[deg]
	double hoo=0.5285;	//minimum thickness[-]
	double ho[N];			//
	double WR[N];			//wetting ratio[-]
	double WRave;			//average of WR[-]
	double sine[N];
	double aa;
	aa=pow(p,3)*pow(g,2)/(15*pow(q,2)*sig);
	double aaa;

	for(j=0 ; j<=N ; j++){			//calculation of ho
		sine[j]=sin(3.14159*(2*j+1)/(2*N));
		aaa=sine[j]*sine[j]*aa;
		ho[j]=hoo/pow(aaa,0.2);
	}
	
//初期値
	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<=N ; j++){
			u[i][j] = 0.03;
		}
	}
	for(j=0 ; j<=N ; j++){
		dr[j] = 0.001/N;
		uave[j]=0.03;
	}
	
	CNewtonRaphsonMethod mnm;		//**************************速度分布******************************************************
	mnm.setup((M+2)*(N+2));					//セットアップ　引数はニュートン法の変数以上にしてください　※2:

	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			mnm.setValue( j+N*i , u[i][j] );
		}
	}


	//mnm.setAcc(0.8);				//加速度勾配の入力　なくても可　※3:

	mnm.initial();					//計算を開始する前に必ず初期化してください
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		//おまじない
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	//おまじない

			for(i=0 ; i<=M ; i++){
				for(j=0 ; j<N ; j++){
					u[i][j] = mnm.getValue(j+N*i);
				}
			}
		
			
			//エラー値
			for(j=0 ; j<N ; j++){
				for(i=0 ; i<=M ; i++){
					sum=sum+u[i][j];
				}
				uave[j]=sum/(M+1);
				sum=0;
			}
			for(j=0 ; j<N ; j++){		//液膜厚さ
				dr[j] = m/(M*p* uave[j]);
			}
			
			
			
			
			for(j=0 ; j<N ; j++){		//蒸気側
				mnm.setError( j , u[0][j] , u[1][j] );
			}

			for(i=1 ; i<M ; i++){
				for(j=0 ; j<N ; j++){
					mnm.setError( j+N*i , g*sin(3.14159*(2*j+1)/(2*N))+c*(u[i+1][j]-2*u[i][j]+u[i-1][j])/(dr[j]*dr[j]) , 0 );
				}
			}

			for(j=0 ; j<N ; j++){		//管側
				mnm.setError( j+N*N , u[M][j] , 0.0 );
			}			

			mnm.prt();				//エラー表示
			mnm.prt_sum();			//エラーの合計を表示
		}
	}//***********************************************************************************************************************
	





	fout << "速度分布" << endl <<",";
	for(j=0 ; j<N ; j++){
		fout<< j << "," ;
	}


	fout << endl;
	for(i=0 ; i<=M ; i++){
		fout << i <<",";
		for(j=0 ; j<N ; j++){
			fout << u[i][j] << ",";
		}
		fout << endl;
	}
	


	
	
	//初期値
	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			T[i][j] = 40.0;	
		}
	}

	CNewtonRaphsonMethod mnm2;		//*******************************温度分布**************************************************
	mnm2.setup((M+1)*(N+1));					//セットアップ　引数はニュートン法の変数以上にしてください　※2:

	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			mnm2.setValue( j+N*i , T[i][j] );
		}
	}
	
	
	//mnm.setAcc(0.8);				//加速度勾配の入力　なくても可　※3:

	mnm2.initial();					//計算を開始する前に必ず初期化してください
	for(mnm2.main_loop_init();mnm2.main_loop_check();mnm2.main_loop_reinit()){		//おまじない
		for(mnm2.sub_loop_init();mnm2.sub_loop_check();mnm2.sub_loop_reinit()){	//おまじない
		
			for(i=0 ; i<=M ; i++){
				for(j=0 ; j<N ; j++){
					T[i][j] = mnm2.getValue(j+N*i);
				}
			}
			int t=0;
			//エラー値
			for(j=0 ; j<N ; j++){		//蒸気側
				mnm2.setError( t , T[0][j] , Tv );
				t++;
			}
			
			for(i=1 ; i<M ; i++){		
				mnm2.setError( t , T[i][0] , Ta );
				t++;
			}
		
			for(i=1 ; i<M ; i++){
				for(j=1 ; j<N-1 ; j++){
					mnm2.setError( t ,(uave[j]*(T[i][j+1]-T[i][j-1]))/(r*(3.14159/M)), a*((T[i+1][j]-2*T[i][j]+T[i-1][j])/(dr[j]*dr[j])+(T[i][j+1]-2*T[i][j]+T[i][j-1])/(r*r*(3.14159/M)*(3.14159/M))) );
					t++;
				}
			}
		
			for(j=0 ; j<N ; j++){		//管側
				mnm2.setError( t , T[M][j] , Ti );
				t++;
			}
			for(i=1 ; i<M ; i++){		//
				mnm2.setError( t , T[i][N-1] , Ti+(Tv-Ti)*(M-i)/M);
				t++;
			}

			mnm2.prt();				//エラー表示
			mnm2.prt_sum();			//エラーの合計を表示
		}
	}



	

	fout << endl << "温度分布" << endl<< ",";
	for(j=0 ; j<N ; j++){
		fout<< j << "," ;
	}
	fout << endl;
	for(i=0 ; i<=M ; i++){
		fout << i <<",";
		for(j=0 ; j<N ; j++){
			fout << T[i][j] << ",";
		}
		fout << endl;
	}

	sum=0;
	fout << endl << "流速" <<endl;
	for(j=0 ; j<N ; j++){
		sum=sum+uave[j];
	}
	uuave=sum/N;

	fout << endl << "液膜厚さ" << endl;
	for(j=0 ; j<N ; j++){
		rr[j]=dr[j]*M;
		fout << rr[j] << "," << endl;
	}

	


	fout << endl << "傾き" << endl;
	for(j=0 ; j<N ; j++){
		b[j]=(T[M-1][j]-T[M][j])/dr[j];
		h[j]=k/(Tv-Ti)*b[j];
		fout<<T[M-1][j]<<"-"<<T[M][j]<<"/"<<dr[j]<<"=";
		fout<< b[j] << endl;
	}




	sum=0;
	for(j=1 ; j<N ; j++){
		sum=sum+b[j];
	}
	cout<< sum;
	bave=sum/(N-1);

	sum=0;
	for(j=1 ; j<N ; j++){
		sum=sum+h[j];
	}
	cout<< sum;
	have=sum/(N-1);

	sum=0;
	for(j=1 ; j<N ; j++){
		sum=sum+rr[j];
	}
	cout<< sum;
	rrave=sum/(N-1);





	//*******************************　厚さhoの時の速度計算　**************************************************
	CNewtonRaphsonMethod mnm3;		//**************************速度分布******************************************************
	mnm3.setup((M+2)*(N+2));					//セットアップ　引数はニュートン法の変数以上にしてください　※2:

	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			mnm3.setValue( j+N*i , u[i][j] );
		}
	}
	//mnm.setAcc(0.8);				//加速度勾配の入力　なくても可　※3:
	mnm3.initial();					//計算を開始する前に必ず初期化してください
	for(mnm3.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		//おまじない
		for(mnm3.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	//おまじない

			for(i=0 ; i<=M ; i++){
				for(j=0 ; j<N ; j++){
					u[i][j] = mnm3.getValue(j+N*i);
				}
			}
		
			
			//エラー値
			for(j=0 ; j<N ; j++){
				for(i=0 ; i<=M ; i++){
					sum=sum+u[i][j];
				}
				uave[j]=sum/(M+1);
				sum=0;
			}
			for(j=0 ; j<N ; j++){		//液膜厚さ
				dr[j] = ho[j]/M;
			}
			
			
			
			
			for(j=0 ; j<N ; j++){		//蒸気側
				mnm3.setError( j , u[0][j] , u[1][j] );
			}

			for(i=1 ; i<M ; i++){
				for(j=0 ; j<N ; j++){
					mnm3.setError( j+N*i , g*sin(3.14159*(2*j+1)/(2*N))+c*(u[i+1][j]-2*u[i][j]+u[i-1][j])/(dr[j]*dr[j]) , 0 );
				}
			}

			for(j=0 ; j<N ; j++){		//管側
				mnm3.setError( j+N*N , u[M][j] , 0 );
			}			

			mnm3.prt();				//エラー表示
			mnm3.prt_sum();			//エラーの合計を表示
		}
	}


	//*******************************　WR計算　**************************************************
	for(j=0 ; j<N ; j++){
		cout<<"ho"<<ho[j]<<endl;
		cout<<"r"<<rr[j]<<endl;

		if(rr[j]>=ho[j]){
			WR[j]=1.0;
		}else{
			width[j]=m/(uave[j]*p*ho[j]);
			cout<<"width"<<width[j]<<endl;
			WR[j]=0.4*width[j];
		}
		fout<<"WR,"<<WR[j]<<endl;
	}

	sum=0;
	for(j=1 ; j<N ; j++){
		sum=sum+WR[j];
	}
	WRave=sum/(N-1);

	fout << endl << "傾き平均,"<<bave;
	fout << endl << "液膜厚さ平均,"<<rrave;
	fout << endl << "流速平均," <<uuave;
	fout << endl << "熱伝達率,"<<have;
	fout << endl << "WR平均,"<<WRave;
	fout.close();

	return 0;
}
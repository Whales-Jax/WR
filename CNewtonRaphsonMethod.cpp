/**************************************************************************************************
*
*	多変数ニュートン法 ： CNewtonRaphsonMethod 
*
*	1G03 武田拓也
*	1G04 大野慶祐
*
*	tool_nc_multi_newton_methodがモジュール解析および起動停止を含めた動特性に耐えられないので
*	カスタマイズしたmulti_newtonの使い勝手をよくするためにさらにカスタマイズしました．
*	ガウスの消去法部分はわからないので全コピペです．
*	[変更点]
*	・分割コンパイル対応．
*	・変数英語化，配列2次元化済み．
*	・関数として独立していたgaussian_elimination,irekaeをクラスの中に格納．
*	・勾配を出す際の2点間の差を最小でもDELTA-1.0になるように修正．
*	・変化率，最大変化量の導入．
*	・変数の数を内部でカウントすることで変数を数えなくてもいいように改良．
*   ・ガウスの消去法に部分ピボッティングを追加
*   ・boostを使っています．インストールしてください．
*	
*	[予定]
*	gaussian_elimination中の__GEMの非マクロ化
*
*	[更新履歴]
*	2008/1/15:	一通りのバグ取り終了．
*   2008/11/20: 部分ピボット選択のガウスの消去法作成
*
*   [boostのインストール方法]
*   １．boostを入手する．共有フォルダにあるので場所はシスカンに聞きましょう
*   ２．c:\lib\boostに回答する．フォルダ構造は以下のとおりになるはずです．
*   c┬[lib]──────[boost]┬[boost]
*  　├[windows]　　　　　　　 ├[doc]
*  　├[prigram files]　　　　 ├[libs]
*    ：　　　　　　　　　　　　：
*   ３．「ツール」→「オプション」→「プロジェクトおよびソリューション」→「VC++ディレクトリ」
*        → 「ディレクトリを表示するプロジェクト」→「インクルード ファイル」に c:\lib\boost を加える。 
*   
*	[基本的な使用法]（裏技多数）
*


#include <iostream>
#include "CNewtonRaphsonMethod.h"

void main(void){

	//初期値
	double a = 5;
	double b = 5;
	double c = 5;
	
	//a * b = 10
	//b * c = 20
	//c * a = 30
	//を解くプログラム
	
	CNewtonRaphsonMethod mnm;		//classを宣言します　※1
	mnm.setup(10);					//セットアップ　引数はニュートン法の変数以上にしてください　※2:

	mnm.setValue( 0 , a );
	mnm.setValue( 1 , b );
	mnm.setValue( 2 , c );

	//mnm.setAcc(0.8);				//加速度勾配の入力　なくても可　※3:

	mnm.initial();					//計算を開始する前に必ず初期化してください
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		//おまじない
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	//おまじない

			a = mnm.getValue(0);
			b = mnm.getValue(1);
			c = mnm.getValue(2);

			//エラー値
			mnm.setError( 0 , a * b , 10 );
			mnm.setError( 1 , b * c , 20 );
			mnm.setError( 2 , c * a , 30 );

			mnm.prt();				//エラー表示
			mnm.prt_sum();			//エラーの合計を表示
		}
	}
	std::cout << " a = " << a << std::endl;
	std::cout << " b = " << b << std::endl;
	std::cout << " c = " << c << std::endl;

}

*	※1:	このmnmという名前が異なれば異なった多変数ニュートン法と認識されます．多変数ニュートン法の中に多変数ニュートン法を
*			構築したいときはここの名前だけ変えたものを大枠の多変数ニュートン法の中にコピペすればＯＫです．
*	※2:	setupする際，必要に応じて微小移動量(いわゆるDELTA)と収束誤差(いわゆるERROR)と最小ピボット値（ガウスの消去法におい
*			て特別な意味をもつ値．通常はこの値を気にする必要はありません．分かっている人だけその値を絶対値で設定してください）
*			を指定できます．指定したい場合はパラメータの数のあとに，微小移動量，収束誤差，最小ピボット値の順で記述します．つま
*			り
*			mnm.setup(2,1+1e-3,1e-3);     //この場合微小移動量:1+1e-3，収束誤差:1e-3//////////////////////////////////////////////
*			mnm.setup(2,DELTA,ERROR,1e-9);//この場合微小移動量:DELTA，収束誤差:ERROR,最小ピボット値:1e-9
*			などという様に書きます．未指定の場合は微小移動量=1.0+1e-6，収束誤差=1e-6，最小ピボット値=-1.0となります．上記の使
*			用法ではパラメータの数=2，微小移動量=1.0+1e-6，収束誤差=1e-6，最小ピボット値=1e-9となります．
*			ピボットが最小ピボット値より小さくなると計算途中でも終了します．（最小ピボット値未指定の場合は途中終了することはあ
*			りません）
*	※3:	加速度勾配を設定することでニュートン法による変化量を抑えることができます．例えば変化率を0.3にしたい場合は
*			mnm.acc = 0.3;
*			としてください．この場合は元の変化量を0.3倍します．デフォルトは1.0です．
*	※4:	mnm.setValueに初期値を代入し，mnm.getValueで値の呼び出しをします．例えば蒸発器内圧力P_vap_eva_inと凝縮器内圧
*			力P_vap_con_inの２つをパラメータとして多変数ニュートン法を行う場合，その各々の初期値を1.0,5.0とすると，
*			mnm.setValue( 0 , 1.0 );//P_vap_eva_in
*			mnm.setValue( 1 , 5.0 );//P_vap_con_in
*			mnm.initial();
*				…
*			P_vap_eva_in = mnm.getValue( 0 );//値の代入
*			P_vap_con_in = mnm.getValue( 1 );
*			となります．
*			初期値に続いて最大変化量を入れることでニュートン法による変化を抑制することができます．これは上述の変化率より優先さ
*			れます．例えば
*			mnm.setValue( 0 , 1.0 , 0.05 );//G_vap_eva_in
*			とするとこのパラメータは最大でも0.05しか変動しません．
*	※5:	mnm.initial()により配列のサイズを初期化します．配列はsetValue関数によりカウントするのでここより前に仮定値を入力
*			してください．
*	※6:	多変数ニュートン法のエラー値はsetError関数により設定します．当然パラメータの数だけsetError関数を記述します．エラ
*			ー値の設定の基本は
*			mnm.setError( エラー番号 , 変数　, 変数　);
*			となりますが，この二つの変数の値が一致するように収束計算が行われます．
*	※7:	エラー表示は個別の表示と全体の表示の2種類が使えます．要素数が少ないときなどはprt()を，多いときはprt_sum()とか使い分
*			けてください．
*
*
*
**************************************************************************************************/
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include "CNewtonRaphsonMethod.h"

double Newton(double X1,double X2,double Y1,double Y2)
{
	return X1-Y1*(X2-X1)/(Y2-Y1);
}

CNewtonRaphsonMethod::CNewtonRaphsonMethod()
{
	FlagSetup = false;

}
CNewtonRaphsonMethod::~CNewtonRaphsonMethod(void)
{

	if( FlagSetup == false ){
		setup(2);
	}

	int i;
	//配列用メモリ解放
	delete [] ErrorMethod;
	delete [] Error;
	delete [] answer;
	delete [] set_value;
	delete [] residual;
	delete [] MaxDisplacement;
	delete [] _MaxDisplacement;
	delete [] assumed_value;
	delete [] _assumed_value;
	delete [] step_value;
	for( i = 0 ; i < ArrayNum+1 ; i++ ){
		delete [] ErrorAbs[i];
		delete [] ErrorRel[i];
	}
	for( i = 0 ; i < ArrayNum ; i++ ){
		delete [] jacob[i];
	}
	delete [] ErrorAbs;
	delete [] ErrorRel;
	delete [] jacob;
}
//仮定値の代入
void CNewtonRaphsonMethod::setup(int n,double _delta,double _err,double _pibot)
{

	FlagSetup = true;

	//仮変数の数，微小移動分，許容エラー値
	_num  = 0;
	delta = _delta;
	err	  = _err;
	pibot = _pibot;
	num   = 0;
	ArrayNum = n;
	FlagInitial  = false;
	acc    = 1.0;
	count_div = 10;
	count_max = 10000;
	EscLoop = 30;
	Solver = 0;
	err_sum = 100.0;
	err_max = 100000;
	ErrorLoopOver = 1;
	MinimumEPS = _delta - 1.0;
	Bprt = false;
	Bprt_sum = false;
	//配列用メモリ確保

	int i ;
	ErrorMethod			= new int[n];
	Error				= new double[n];
	answer				= new double[n];
	residual			= new double[n];
	set_value			= new double[n];
	MaxDisplacement		= new double[n];
	_MaxDisplacement	= new double[n];
	MaxValue			= new double[n];
	MinValue			= new double[n];
	assumed_value		= new double[n];
	_assumed_value		= new double[n];
	step_value			= new double[n];
	jacob				= new double*[n];
	ErrorAbs			= new double*[n+1];
	ErrorRel			= new double*[n+1];
	for( i = 0 ; i < n+1 ; i++ ){
		ErrorAbs[i] = new double[n];
		ErrorRel[i] = new double[n];
	}
	for( i = 0 ; i < n ; i++ ){
		jacob[i] = new double[n];
	}

	for( i = 0 ; i < n ; i++ ){
		ErrorMethod[i] = 0;
		Error[i] = err;
	}

}
void CNewtonRaphsonMethod::setDeltaError( double _delta , double _err ){
	delta = _delta;
	err	  = _err;
	MinimumEPS = _delta - 1.0;
	for( int i = 0 ; i < ArrayNum ; i++ ){
		Error[i] = err;
	}
}

void CNewtonRaphsonMethod::setErrorMethod( int _ErrorMethod ){
	int i ;
	for( i = 0 ; i < ArrayNum ; i++ ){
		ErrorMethod[i] = _ErrorMethod;
	}
}
void CNewtonRaphsonMethod::setAcc( double _acc ){
	acc = _acc;
}
void CNewtonRaphsonMethod::setSolver( int _Solver ){
	Solver = _Solver;
}
void CNewtonRaphsonMethod::setEscLoop( int _EscLoop ){
	EscLoop = _EscLoop;
}
void CNewtonRaphsonMethod::reset(){
	FlagInitial = false;
	_num = 0;
}

void CNewtonRaphsonMethod::set_assumed_value(int i,double value,double _max)
{

	if( FlagInitial == true ){
		reset();
	}

	if( i+1 > _num ){
		_num = i+1;
	}
	_assumed_value[i] = value;
	_MaxDisplacement[i] = _max;

}
void CNewtonRaphsonMethod::setValue(int i,double value,double _max , double _MinValue , double _MaxValue )
{

	if( FlagInitial == true ){
		reset();
	}

	if( i+1 > _num ){
		_num = i+1;
	}
	_assumed_value[i] = value;
	_MaxDisplacement[i] = _max;

	MinValue[i] = _MinValue;
	MaxValue[i] = _MaxValue;


}
double CNewtonRaphsonMethod::getValue(int i){
	return set_value[i];
}

void CNewtonRaphsonMethod::initial(int n)
{
	int i;
	if( n == 0 ){
		num = _num;
	}else{
		num = n;
	}

//	matA = boost::numeric::ublas::matrix<double>( num , num );
//	vecX = boost::numeric::ublas::vector<double>( num );
//	vecB = boost::numeric::ublas::vector<double>( num );

	if( FlagInitial == false ){
		FlagInitial = true;
		for( i = 0 ; i < num ; i++ ){
			assumed_value[i] = _assumed_value[i];
			MaxDisplacement[i]	 = _MaxDisplacement[i];
		}
	}

	ErrorRel[num][0] = 100000000.0;
	ErrorAbs[num][0] = 100000000.0;
	err_sum = num * 100;
	err_ave = err_sum / num;
	ErrorLoopOver = 1;
}
void CNewtonRaphsonMethod::main_loop_init()
{
	//値の初期化
	if( FlagInitial == false ){
		std::cout<<"error : initialize"<<std::endl;
		getchar();
		exit(1);
	}
	check = false;
	check2 = false;
	count = 1;

	int i,j;
	for( i = 0 ; i < ArrayNum+1 ; i++ ){
		for( j = 0 ; j < ArrayNum ; j++ ){
			ErrorAbs[i][j] = 100.0;
			ErrorRel[i][j] = 100.0;
		}
	}

}
bool CNewtonRaphsonMethod::main_loop_check()//これが最後の処理
{


	if( count > EscLoop ){
		std::cout<< name << "loop max skip "<< "error_max = " << err_max << std::endl;
		ErrorLoopOver = 0;
		return (false);
	}

	return (!check);

/*
	//合計エラーの計算
	err_max = 0.0;
	err_sum = 0.0;
	for(int i=0;i<num;i++){
		if( ErrorMethod[i] == 0 ){
			err_sum  += fabs(ErrorRel[num][i]);
			if( fabs(ErrorRel[num][i]) > err_max ){
				err_max = fabs(ErrorRel[num][i]);
			}
		}else if( ErrorMethod[i] == 1 ){
			err_sum  += fabs(ErrorAbs[num][i]);
			if( fabs(ErrorAbs[num][i]) > err_max ){
				err_max = fabs(ErrorAbs[num][i]);
			}
		}
	}
	err_ave = err_sum / num;


	if( Bprt == true ){
		prt2();
	}
	if( Bprt_sum == true ){
		prt_sum2();
	}


	//収束判断
	check = true;
	if( count > EscLoop ){
		std::cout<< name << "loop max skip "<< "error_max = " << err_max << std::endl;
		ErrorLoopOver = 0;
		return (!check);
	}

	for( int i=0 ; i <= num-1 ; i++ )
	{
		if( ErrorMethod[i] == 0 ){
			if( fabs(ErrorRel[num][i] ) >= Error[i] ){
				check = false;
				break;
			}
		}else if( ErrorMethod[i] == 1 ){
			if( fabs(ErrorAbs[num][i] ) >= Error[i] ){
				check = false;
				break;
			}
		}

	}
	return (!check);
	*/

}
void CNewtonRaphsonMethod::main_loop_reinit()
{
	int i;
	//合計エラーの計算
	err_max = 0.0;
	err_sum = 0.0;
	for(i=0;i<num;i++){
		if( ErrorMethod[i] == 0 ){
			err_sum  += fabs(ErrorRel[num][i]);
			if( fabs(ErrorRel[num][i]) > err_max ){
				err_max = fabs(ErrorRel[num][i]);
			}
		}else if( ErrorMethod[i] == 1 ){
			err_sum  += fabs(ErrorAbs[num][i]);
			if( fabs(ErrorAbs[num][i]) > err_max ){
				err_max = fabs(ErrorAbs[num][i]);
			}
		}
	}
	err_ave = err_sum / num;

	if( Bprt == true ){
		prt2();
	}
	if( Bprt_sum == true ){
		prt_sum2();
	}

	//収束判断
	check = true;
	for(i=0 ; i <= num-1 ; i++ )
	{
		if( ErrorMethod[i] == 0 ){
			if( fabs(ErrorRel[num][i] ) >= Error[i] ){
				check = false;
				break;
			}
		}else if( ErrorMethod[i] == 1 ){
			if( fabs(ErrorAbs[num][i] ) >= Error[i] ){
				check = false;
				break;
			}
		}
	}
	if( check == true ){
		return;//収束していた場合は修正値を計算せずにmain_loop_check()へ
	}



	int gyou,retu;
	double tmp,sum;
	//偏微分の算出
	for( gyou = 0 ; gyou <= num-1 ; gyou++ )
	{
		residual[gyou] = ErrorAbs[num][gyou];

		for( retu = 0 ; retu <= num-1 ; retu++ )
		{
			jacob[gyou][retu] = ( ErrorAbs[num][gyou] - ErrorAbs[retu][gyou] )
				/ (assumed_value[retu] - step_value[retu] );
		}
	}

	//連立一次方程式
	if( Solver == 0 ){
		ErrorSolver = gaussian_elimination();
	}else if( Solver == 1 ){
		ErrorSolver = gaussian_elimination_2();
	}else if( Solver == 2 ){
//		ErrorSolver = lu_factorize_ublas();
	}else{
		ErrorSolver = gaussian_elimination();
	}


	if( ErrorSolver == 0 )
	{
		std::cout<<"error : solver"<<std::endl;
		getchar();
		exit(1);
	}
	//値の変更
	for( i=0; i<=num-1; i++ )
	{
		if(MaxDisplacement[i]==0.0){
			assumed_value[i] = set_value[i] - acc*answer[i];
		}
		else{
			double temp2;
			temp2 = acc * answer[i];
			tmp = fabs(temp2);
			if( tmp > MaxDisplacement[i] ){
				if( answer[i]>0 ){
					assumed_value[i] = set_value[i] - MaxDisplacement[i];
				}else{
					assumed_value[i] = set_value[i] + MaxDisplacement[i];
				}
			}else{
				assumed_value[i] = set_value[i] -  temp2;
			}
		}

		if( MinValue[i] != 0.0 && assumed_value[i] < MinValue[i] ){
			assumed_value[i] = MinValue[i];
		}
		if( MaxValue[i] != 0.0 && assumed_value[i] > MaxValue[i] ){
			assumed_value[i] = MaxValue[i];
		}



	}
	//カウンタをあげる
	if(count==count_max){
		sum = 0;
		for(i=0;i<num;i++){
			sum += fabs(ErrorAbs[num][i]);
		}
		if(sum<sqrt(err)){
			std::cout<<"loop over"<<std::endl;
			check=true;
		}
		else std::cout<<"error : loop over"<<std::endl;
		getchar();
		exit(1);
	}



	count++;
	
}
void CNewtonRaphsonMethod::spcacc(double maxacc){

	acc = maxacc * pow( 1.5 , count - 10.0 );
	if( acc > maxacc ){
		acc = maxacc;
	}
}
void CNewtonRaphsonMethod::sub_loop_init()
{
	//値の初期化
	i_sys   = 0;
	i_error = 0;
//	double MinimumEPS = 1.0;
	//値の代入
	for(int i=0; i<=num-1; i++ )
	{
		if(fabs(assumed_value[i])>MinimumEPS)	step_value[i] = assumed_value[i]*delta;
		else									step_value[i] = assumed_value[i]+(assumed_value[i]>0?MinimumEPS:-MinimumEPS)*(delta-1.0);
//		step_value[i] = assumed_value[i]*delta;

		set_value[i] = (i==i_sys)?step_value[i]:assumed_value[i];
	}
}
bool CNewtonRaphsonMethod::sub_loop_check()
{
	return ( i_sys <= num );
}
void CNewtonRaphsonMethod::sub_loop_reinit()
{
	i_sys++;
	//値の初期化
	i_error = 0;
	//値の代入
	if( i_sys != num + 1 )
	{
		for(int i=0; i<=num-1; i++ )
		{
			set_value[i] = (i==i_sys)?step_value[i]:assumed_value[i];
		}
	}

}

//エラー値の代入
void CNewtonRaphsonMethod::set_error(double atai,double bairitu)
{
	ErrorAbs[i_sys][i_error] = atai*bairitu;
	ErrorRel[i_sys][i_error] = ErrorAbs[i_sys][i_error];
	i_error++;
}
void CNewtonRaphsonMethod::set_error2(int bangou,double atai1,double atai2,double bairitu)
{
	ErrorAbs[i_sys][bangou] = ( atai1 - atai2 ) * bairitu;
/*
	if( fabs(atai1) < 1e-3 || fabs(atai2) < 1e-3 ){
		if( atai1 > 0.0 ){
			atai1 += 1e-3;
			atai2 += 1e-3;
		}else{
			atai1 -= 1e-3;
			atai2 -= 1e-3;
		}
	}*/

	ErrorRel[i_sys][bangou] = ( 1.0 - ( atai1 / atai2 ) ) * bairitu;
}

void CNewtonRaphsonMethod::setError(int i,double Value1,double Value2, int _ErrorMethod , double error )
{
	if( error > 0.0 ){
		Error[i] = error;
	}else{
		Error[i] = err;
	}

	ErrorMethod[i] = _ErrorMethod;

	ErrorAbs[i_sys][i] = ( Value1 - Value2 );

	//この処理が怪しいのでコメント化　by大野
	if( fabs(Value1) < 1e-3 || fabs(Value2) < 1e-3 ){
		if( Value2 > 0.0 ){
			Value1 += 1e-3;
			Value2 += 1e-3;
		}else{
			Value1 -= 1e-3;
			Value2 -= 1e-3;
		}
	}
	
	ErrorRel[i_sys][i] = ( 1.0 - ( Value1 / Value2 ) );


//	if( ErrorRel[i_sys][i] < -100000.0 || 100000 < ErrorRel[i_sys][i] ){
//		std::cout << "asdf";
//	}


}


//エラー値の表示
void CNewtonRaphsonMethod::prt(double value)
{
	int i;
	if( i_sys == num )
	{
		std::cout<<"\tnum of loop = "<<count<<std::endl;
		for( i=0; i<=num-1; i++ )
		{
			if( fabs(ErrorAbs[i_sys][i]) > value ){
				std::cout << std::setprecision(6);
				std::cout << "ErrorRel[" << i+1 << "] = " << fabs(ErrorRel[i_sys][i]) << "\t";
				std::cout << "ErrorAbs[" << i+1 << "] = " << fabs(ErrorAbs[i_sys][i]) << "\t";
				std::cout << std::endl;
			}
		}
	}
}
void CNewtonRaphsonMethod::prt_sum()
{
	if(i_sys==num)
	{
		std::cout.precision(6);
		std::cout << "\tnum of loop = "<< count; 
		std::cout << "\tacc = "<< acc;
		std::cout << "\terr_sum = "<< err_sum;
		std::cout << "\terr_ave = " << err_ave;
		std::cout << "\terr_max = " << err_max;
		std::cout << std::endl;
	}
}
void CNewtonRaphsonMethod::prt2(double value)
{
	int i;
	std::cout<<"\tnum of loop = "<<count<<std::endl;
	for( i=0; i<=num-1; i++ )
	{
		if( fabs(ErrorAbs[i_sys][i]) > value ){
			std::cout << std::setprecision(6);
			std::cout << "ErrorRel[" << i+1 << "] = " << fabs(ErrorRel[i_sys][i]) << "\t";
			std::cout << "ErrorAbs[" << i+1 << "] = " << fabs(ErrorAbs[i_sys][i]) << "\t";
			std::cout << std::endl;
		}
	}

}
void CNewtonRaphsonMethod::prt_sum2()
{

	std::cout.precision(6);
	std::cout << "\tnum of loop = "<< count; 
	std::cout << "\tacc = "<< acc;
	std::cout << "\terr_sum = "<< err_sum;
	std::cout << "\terr_ave = " << err_ave;
	std::cout << "\terr_max = " << err_max;
	std::cout << std::endl;

}

void CNewtonRaphsonMethod::matrixFileMake(){
	
	matrixFile.open("matrix.csv" , std::ios::out );

}
void CNewtonRaphsonMethod::matrixOut(){
	


	matrixFile << "jacob and answer" << std::endl;
	for( int i = 0 ; i < num ; i++ ){
		for( int j = 0 ; j < num ; j++ ){
			matrixFile << jacob[i][j] << ",";
		}
		matrixFile << answer[i] << ",";
		matrixFile << std::endl;
	}
	matrixFile << std::endl;

	matrixFile << "REL error" << std::endl;
	for( int i = 0 ; i < num ; i++ ){
		for( int j = 0 ; j < num ; j++ ){
			matrixFile << ErrorRel[i][j] << ",";
		}
		matrixFile << std::endl;
	}
	matrixFile << std::endl;

	matrixFile << "ABS error" << std::endl;
	for( int i = 0 ; i < num ; i++ ){
		for( int j = 0 ; j < num ; j++ ){
			matrixFile << ErrorAbs[i][j] << ",";
		}
		matrixFile << std::endl;
	}
	matrixFile << std::endl;


}

//ガウスの消去法　部分ピボット選択ver
int CNewtonRaphsonMethod::gaussian_elimination_2(){
	int i,j,k;
	double max;
	double keisu;
	int si,sj,maxi;//ピボット選択時に使用
	for( k = 0 ; k < num ; k++ ){
		//ピボット選択
		max = 0.0;
		for( si = k ; si < num ; si++ ){
			if( max <= fabs( jacob[si][k] ) ){
				maxi = si;
				max  = fabs( jacob[si][k] );
			}
			if( max > 1.0 ){
				break;
			}
		}
		if( max <= pibot ){
			return(1);
		}
		//入れ替え
		if( maxi != k ){
			for( sj = 0 ; sj < num ; sj++ ){
				irekae( jacob[k][sj] , jacob[maxi][sj] );
			}
			irekae( residual[k] , residual[maxi] );
		}
		for( i = k+1 ; i < num ; i++ ){
			if( fabs(jacob[i][k]) > 1e-12 ){
				keisu = jacob[i][k] / jacob[k][k];
				for( j = k+1 ; j < num ; j++ ){
					jacob[i][j] -= jacob[k][j] * keisu; 
				}
				residual[i] -= residual[k] * keisu;
				jacob[i][k] = 0.0;
			}
		}
	}
	//後退代入
	for( k = num-1 ; k >= 1 ; k-- ){
		for( i = k-1 ; i >= 0 ; i-- ){
			residual[i] -= residual[k] * jacob[i][k] / jacob[k][k];
			jacob[i][k] = 0.0;
		}
	}
	//答え計算
	for( k = 0 ; k < num ; k++ ){
		answer[k] = residual[k] / jacob[k][k];
	}
	return 1;
}
/*
int CNewtonRaphsonMethod::lu_factorize_ublas(){
	int gyou;
	int retu;
	for( gyou = 0 ; gyou < num ; gyou++ ){
		for( retu = 0 ; retu < num ; retu++ ){
			matA(gyou,retu) = jacob[gyou][retu];
		}
	}
	for( gyou = 0 ; gyou < num ; gyou++ ){
		vecB(gyou) = residual[gyou];
	}
	boost::numeric::ublas::matrix<double> lhs(matA);
	boost::numeric::ublas::vector<double> rhs(vecB);
	boost::numeric::ublas::permutation_matrix<> pm(matA.size1());
	boost::numeric::ublas::lu_factorize(lhs, pm);
	boost::numeric::ublas::lu_substitute(lhs, pm, rhs);
	vecX.assign_temporary(rhs);
	for( gyou = 0 ; gyou < num ; gyou++ ){
		answer[gyou] = vecX(gyou);
	}

	return 1;
}
*/




/*
 *  Numerical Calculation : ガウスの消去法（完全ピボッティング） ver5.0 G98堀内雄次
 *
 *  keisuu * x = y （keisuu,y:既知,x:未知）を解く．
 *
 *  ─────────┬───┬──
 *  引数              │型    │備考
 *  ─────────┼───┼──
 *  係数行列（keisuu）│type *│*1
 *  定数行列(y)       │type *│*2
 *  解(x)             │type *│*2
 *  係数行列の次元(n) │int   │
 *  最小ピボット値    │type  │
 *  ─────────┴───┴──
 *  *1:係数行列の次元と配列の次元を等しくしてください．
 *     つまり2x2行列ならa[2][2]などをいれてa[5][5]などは入れないで下さい
 */

void CNewtonRaphsonMethod::irekae(double &x,double &y ){double z=x;x=y;y=z;}
void CNewtonRaphsonMethod::irekae(int    &x,int    &y ){int    z=x;x=y;y=z;}

#define __GEM(i,j)  ( (j!=num) ? (jacob[i][j]) : (residual[i]) )
int CNewtonRaphsonMethod::gaussian_elimination()
{
	int i,j,k;
	double max;
	int si,sj,maxi,maxj;//ピボット選択時に使用
	int *x_narabi = new int [num];//解の並び順
	for( i=0; i<=num-1; i++ )
	{
		x_narabi[i]=i;
	}
	double *x_tmp = new double [num];//解を一時的にストック

	for( k=0; k<=num-1; k++ )
	{
		//完全ピボット選択
		max = 0.0;
		for( si=k; si<=num-1; si++ )
		{
			for( sj=k; sj<=num-1; sj++ )
			{
				if( max <= fabs( __GEM(si,sj) ) )
				{
					max  = fabs( __GEM(si,sj) );
					maxi = si;
					maxj = sj;
				}
			}
		}
		if( max <= pibot )
		{
			return (1);
		}
		if( maxi != k )//行入れ換え
		{
			for( sj=k; sj<=num; sj++ )
			{
				irekae( __GEM(k,sj), __GEM(maxi,sj) );
			}
		}
		if( maxj != k )//列入れ換え
		{
			for( si=0; si<=num-1; si++ )
			{
				irekae( __GEM(si,k), __GEM(si,maxj) );
			}
			irekae( x_narabi[k], x_narabi[maxj] );
		}

	
		//対象となる行を操作
		for( j=k+1; j<=num; j++ )
		{
			__GEM(k,j) /= __GEM(k,k);
		}
		//対象となる行より下方にある行を操作
		for( i=k+1; i<=num-1; i++ )
		{
			if( __GEM(i,k) != 0.0 ){
				for( j = k+1 ; j < num+1 ; j++ )
				{
					__GEM(i,j) -= __GEM(i,k) *__GEM(k,j);
				}
			}
		}
	}
	//解の算出
	for( k=num-1; k>=0; k-- )
	{
		x_tmp[k] = residual[k];
		for( j=k+1; j<=num-1; j++ )
		{
			x_tmp[k] -= x_tmp[j] * __GEM(k,j);
		}
	}
	for( i=0; i<=num-1; i++ )
	{
		answer[ x_narabi[i] ] = x_tmp[i];
	}
	//メモリ解放
	delete [] x_narabi;
	delete [] x_tmp;
	return 1;
}

#undef __GEM
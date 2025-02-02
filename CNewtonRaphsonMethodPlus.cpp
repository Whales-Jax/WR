/**************************************************************************************************
*
*	多変数ニュートン法 ： CNewtonRaphsonMethodPlus
*
*	1G04 大野慶祐
*
*	CNewtonRaphsonMethodが並列計算やヤコビ行列作成スキップ機能に対応しにくいため，作り直した．
*   
*   [主な変更点]
*   ・openmpによる並列計算に対応
*   ・ヤコビ行列の作成をスキップする機能に対応
*   ・内部の行列はvector型で宣言
*   
*	[基本的な使用法]
*

#include <iostream>
#include "CNewtonRaphsonMethodPlus.h"

void main(){

	//クラスを宣言
	CNewtonRaphsonMethodPlus nrm;

	double a;
	double b;
	double c;

//------------------------以下ノーマルな使用方法---------------------------//

	a = 1.0;
	b = 1.0;
	c = 1.0;

	//初期化
	nrm.Initialize();

	//変数のセット
	nrm.SetValiable( 0 , a );
	nrm.SetValiable( 1 , b );
	nrm.SetValiable( 2 , c );

	//各種セットアップ

	//加速係数の設定
	nrm.SetAcc( 0.5 );

	//1ループごとのエラー値のディスプレイ設定
	//nrm.SetPrint( 全体エラー表示 , 変数ごとのエラー表示 );
	nrm.SetPrint( true , true );

	//微小移動量の設定
	nrm.SetDelta( 1.0001 );

	//収束判定の設定
	nrm.SetError( 0.0001 );

	for(nrm.MainLoopInit();nrm.MainLoopCheck();nrm.MainLoopReinit()){
		for(nrm.SubLoopInit();nrm.SubLoopCheck();nrm.SubLoopReinit()){
			
			//修正値の受け取り
			nrm.GetValiable( 0 , a );
			nrm.GetValiable( 1 , b );
			nrm.GetValiable( 2 , c );

			//エラー値の代入
			nrm.SetResult( 0 , a * b , 10.0 );
			nrm.SetResult( 1 , b * c , 20.0 );
			nrm.SetResult( 2 , c * a , 30.0 );
		
		}
	}

	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << "c = " << c << std::endl;

//------------------------以上ノーマルな使用方法---------------------------//


//------------------------以下アドバンスドな使用方法---------------------------//

	a = 1.0;
	b = 1.0;
	c = 1.0;

	//初期化
	nrm.Initialize();

	//変数のセット
	nrm.SetValiable( 0 , a );
	nrm.SetValiable( 1 , b );
	nrm.SetValiable( 2 , c );

	//各種セットアップ

	//加速係数の設定
	nrm.SetAcc( 0.5 );

	//1ループごとのエラー値のディスプレイ設定
	nrm.SetPrint( false , false );

	//微小移動量の設定
	nrm.SetDelta( 1.0001 );

	//収束判定の設定
	nrm.SetError( 0.000001 );

	//ヤコブ行列の作成をスキップするステップ数を設定．
	//この場合では，ニュートン法5ループは同じヤコビ行列を使用する
	nrm.SetJacobStep( 5 );

	//ニュートン法の終了時に最終エラーをディスプレイに表示する
	nrm.SetEndReport( true ) ;


	for(nrm.MainLoopInit();nrm.MainLoopCheck();nrm.MainLoopReinit()){
		for(nrm.SubLoopInit();nrm.SubLoopCheck();nrm.SubLoopReinit()){
			
			//修正値の受け取り
			nrm.GetValiable( 0 , a );
			nrm.GetValiable( 1 , b );
			nrm.GetValiable( 2 , c );

			//エラー値の代入
			nrm.SetResult( 0 , a * b , 10.0 );
			nrm.SetResult( 1 , b * c , 20.0 );
			nrm.SetResult( 2 , c * a , 30.0 );
		
		}
	}

	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << "c = " << c << std::endl;

//------------------------以以上アドバンスドな使用方法---------------------------//
	
}

*	以下の文章はCNewtonRaphsonMethodのコピペですが，CNewtonRaphsonMethodPlusでもほぼ同じ機能をそろえています．
*
*	※1:	このnrmという名前が異なれば異なった多変数ニュートン法と認識されます．多変数ニュートン法の中に多変数ニュートン法を
*			構築したいときはここの名前だけ変えたものを大枠の多変数ニュートン法の中にコピペすればＯＫです．
*	※2:	setupする際，必要に応じて微小移動量(いわゆるDELTA)と収束誤差(いわゆるERROR)と最小ピボット値（ガウスの消去法におい
*			て特別な意味をもつ値．通常はこの値を気にする必要はありません．分かっている人だけその値を絶対値で設定してください）
*			を指定できます．指定したい場合はパラメータの数のあとに，微小移動量，収束誤差，最小ピボット値の順で記述します．つま
*			り
*			mnm.setup(2,1+1e-3,1e-3);     //この場合微小移動量:1+1e-3，収束誤差:1e-3
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
#include "CNewtonRaphsonMethodPlus.h"


#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)

/*
memo

for( 1 init ; 2 check ; 3 reinit ){
	4
}

12/432/432/432/432/432/43......


1初期化
2抜け出す判定→仮定値の計算
4モジュール計算
3収束判定→ガウスの消去法→修正値の計算




造るものリスト
・仮定値の履歴の保存
・説明文
・関数の注釈
・keisuを配列化する
・式と変数の数が合わなかった時のエラー
・サブループ回数を指定すると変数を出してくれる関数

*/





#include "CNewtonRaphsonMethodPlus.h"



int CNewtonRaphsonMethodPlus::MatrixOut(){


	int i;
	int j;


	MatrixFile.open(Name + "matrix.csv", std::ios::out);

	if (!MatrixFile){

		std::cout << "file open error" << std::endl;
		std::cout << Name << "のマトリックスファイルが開けません" << std::endl;
		system("pause");
		exit(0);

	}


	//ヤコビ行列のアウトプット

	MatrixFile << "JacobMatrix" << std::endl;

	MatrixFile << "#" << ",";
	for (j = 0; j < NumVariable; j++){
		MatrixFile << "#eq" << j << ",";
	}
	MatrixFile << "Results"<<std::endl;

	for ( i = 0; i < NumVariable; i++){
		MatrixFile << "#Variable_" << i << "_" << VariableName[i] << ",";
		for (j = 0; j < NumVariable + 1; j++){
			MatrixFile << JacobMatrix[i][j] << ",";
		}
		MatrixFile << std::endl;
	}
	MatrixFile << std::endl;





	//上三角行列のアウトプット

	MatrixFile << "UpperTriangularMatrix" << std::endl;

	MatrixFile << "#" << ",";
	for (j = 0; j < NumVariable; j++){
		MatrixFile << "#eq" << j << ",";
	}
	MatrixFile << "Results" << std::endl;

	for (i = 0; i < NumVariable; i++){
		MatrixFile << "#Variable_" << i << "_" <<VariableName[i] << ",";
		for ( j = 0; j < NumVariable + 1; j++){
			MatrixFile << UpperTriangularMatrix[i][j] << ",";
		}
		MatrixFile << std::endl;
	}
	MatrixFile << std::endl;



	//対角行列のアウトプット

	MatrixFile << "DiagonalMatrix" << std::endl;

	MatrixFile << "#" << ",";
	for (j = 0; j < NumVariable; j++){
		MatrixFile << "#eq" << j << ",";
	}
	MatrixFile << "Results" << std::endl;

	for (i = 0; i < NumVariable; i++){
		MatrixFile << "#Variable_" << i << "_" <<VariableName[i] << ",";
		for ( j = 0; j < NumVariable + 1; j++){
			MatrixFile << DiagonalMatrix[i][j] << ",";
		}
		MatrixFile << std::endl;
	}
	MatrixFile << std::endl;






	MatrixFile.close();



	return 1;
}



int CNewtonRaphsonMethodPlus::SetName( std::string _Name ){
	Name = _Name;
	return 1;
}

/*
 * 初期化を行う
 *
 *
 */
int CNewtonRaphsonMethodPlus::Initialize(){
	NumVariable = 0;
	FlagJacobSkip = false;
	return 1;
}

int CNewtonRaphsonMethodPlus::SetAccAuto(bool _FlagAccAuto, double _AccAutoMin , double _AccAutoMax , double _AccAutoStep ){
	FlagAccAuto = _FlagAccAuto;
	AccAutoMin = _AccAutoMin;
	AccAutoMax = _AccAutoMax;
	AccAutoStep = _AccAutoStep;


	return 1;

}

int CNewtonRaphsonMethodPlus::SetAcc( double _Acc ){
	int i;
	for( i = 0 ; i < NumVariable ; i++ ){
		SetAccDirect( i , _Acc );
	}
	return 1;
}

int CNewtonRaphsonMethodPlus::SetAccDirect( int _Num , double _Acc ){
	Acc[ _Num ] = _Acc;

	return 1;
}

int CNewtonRaphsonMethodPlus::SetEndReport( bool _Flag ){
	FlagEndReport = _Flag;
	return 1;
}

int CNewtonRaphsonMethodPlus::ValueContainerResize( int _NumVariable ){
	Variable.resize( _NumVariable+1 );//変数それ自体
	VariableName.resize( _NumVariable + 1 );
	VariableMemory.resize( _NumVariable+1 );//変数の履歴


	Acc.resize( _NumVariable+1 );//変数の加速係数
	MaxValue.resize( _NumVariable+1 );//変数の最大値
	MinValue.resize( _NumVariable+1 );//変数の最小値
	MaxDisplacement.resize( _NumVariable+1 );//変数の最大移動量

	MaxValue[_NumVariable] = 0.0;
	MinValue[_NumVariable] = 0.0;
	MaxDisplacement[_NumVariable] = 0.0;

	StepControl.resize( _NumVariable+1 );
	MinMaxControl.resize( _NumVariable+1 );
	StepControl[_NumVariable] = false;
	MinMaxControl[_NumVariable] = false;

	Step.resize( _NumVariable+1 );//変数の最大移動量


	for( int i = NumVariable ; i < _NumVariable+1 ; i++ ){
		Variable[ i ].resize(1);
		VariableMemory[ i ].resize(1);
		Acc[ i ] = 1.0;//Accは1で初期化する．
	}
	return 1;
}

int CNewtonRaphsonMethodPlus::SetValiable( int _NumVariable , double _Variable ){

	if( NumVariable < _NumVariable+1 ){
		ValueContainerResize( _NumVariable );
		NumVariable = _NumVariable+1;
	}
	
	Variable[ _NumVariable ][ 0 ] = _Variable;
	VariableMemory[ _NumVariable ][ 0 ] = _Variable;
//	StepControl[_NumVariable] = false;
//	MinMaxControl[_NumVariable] = false;	
	return 1;
}


int CNewtonRaphsonMethodPlus::SetValiableName(int _NumVariable, std::string _VariableName){

	VariableName[_NumVariable]= _VariableName;

	return 1;
}


int CNewtonRaphsonMethodPlus::SetValiable( int _NumVariable , double _Variable , double _MaxStep ){
	SetValiable(  _NumVariable ,  _Variable );
	SetMaxStep(  _NumVariable ,  _MaxStep );
	return 1;
}

int CNewtonRaphsonMethodPlus::SetValiable( int _NumVariable , double _Variable , double _MaxStep , double _MinValue , double _MaxValue ){
	SetValiable(  _NumVariable ,  _Variable );
	SetMaxStep(  _NumVariable ,  _MaxStep );
	SetMinMaxValue(  _NumVariable ,  _MinValue , _MaxValue );
	return 1;
}

int CNewtonRaphsonMethodPlus::SetMaxStep( int _NumVariable , double _MaxStep ){
	MaxDisplacement[ _NumVariable ] = _MaxStep;
//	StepControl[_NumVariable] = true;

	return 1;
}
int CNewtonRaphsonMethodPlus::SetMinMaxValue( int _NumVariable , double _MinValue , double _MaxValue ){
	MinValue[ _NumVariable ] = _MinValue;
	MaxValue[ _NumVariable ] = _MaxValue;
//	MinMaxControl[_NumVariable] = true;
	return 1;
}


int CNewtonRaphsonMethodPlus::GetValiable( int _NumVariable , double &_Variable ){
	_Variable = Variable[ _NumVariable ][ SubLoopCounter ];
	return 1;
}
int CNewtonRaphsonMethodPlus::GetValiableDirect( int _NumVariable , int _SubLoopCounter, double &_Variable ){
	_Variable = Variable[ _NumVariable ][ _SubLoopCounter ];
	return 1;
}

int CNewtonRaphsonMethodPlus::SetResult( int _NumVariable , double _VariableA , double _VariableB ){
	ResultMatrix[ _NumVariable ][ SubLoopCounter ] = _VariableA - _VariableB;
	return 1;
}
int CNewtonRaphsonMethodPlus::SetResultDirect( int _NumVariable , int _SubLoopCounter , double _VariableA , double _VariableB ){
	ResultMatrix[ _NumVariable ][ _SubLoopCounter ] = _VariableA - _VariableB;
	return 1;
}


int CNewtonRaphsonMethodPlus::SetPrint( bool _Print , bool _PrintAll ){
	FlagPrint = _Print;
	FlagPrintAll = _PrintAll;
	return 1;
}

int CNewtonRaphsonMethodPlus::SetMatrixOut(bool _MatrixOut){
	FlagMatrixOut = _MatrixOut;
	return 1;
}



int CNewtonRaphsonMethodPlus::MainLoopInit(){
	int i;


	time(&TimeStart);

	FlagConbergence = false;
	MainLoopCounter = 0;
	JacobSkipCounter = 0;
	JacobReset = false;



	//配列のサイズを決める

	for( i = 0 ; i < NumVariable ; i++ ){
		Variable[i].resize( NumVariable+1 );
	}


	ResultMatrix.resize( NumVariable );
	JacobMatrix.resize( NumVariable );
	UpperTriangularMatrix.resize( NumVariable );
	DiagonalMatrix.resize( NumVariable );
	UpperC.resize( NumVariable );

	RelError.resize( NumVariable );
	AbsError.resize( NumVariable );
	DeltaValue.resize( NumVariable );

	for( i = 0 ; i < NumVariable ; i++ ){
		ResultMatrix[i].resize( NumVariable+1 );
		JacobMatrix[i].resize( NumVariable+1 );
		UpperTriangularMatrix[i].resize( NumVariable+1 );
		DiagonalMatrix[i].resize( NumVariable+1 );
	}


	//最初の収束判定で飛ばされないようにエラー値に100を入れておく
	ResultMatrix[0][NumVariable] = 100.0;

	for( i = 0 ; i < NumVariable ; i++ ){
		Variable[i][ NumVariable ] = Variable[i][0];
		DiagonalMatrix[i][ NumVariable ] = 0.0;
	}

	return 1;
}

int CNewtonRaphsonMethodPlus::MainLoopCheck(){

	//Errorの値を見る
	MaxError = 0.0;
	AveError = 0.0;
	for( int i = 0 ; i < NumVariable ; i++ ){
		if( MaxError < fabs( ResultMatrix[i][NumVariable] ) ){
			MaxError = fabs( ResultMatrix[i][NumVariable] );
		}
		AveError += fabs( ResultMatrix[i][NumVariable] );
	}
	AveError /= NumVariable;


	//errorによるAccの自動調整
	ErrorMax0 = ErrorMax1;
	ErrorMax1 = MaxError;

	ErrorDif0 = ErrorDif0;
	ErrorDif1 = ErrorMax1 - ErrorMax0;

	if (FlagAccAuto){
		if (ErrorDif1 < 0.0){
		}else{
			double _Acc = Acc[0];
			if (_Acc / AccAutoStep < AccAutoMin){
				SetAcc( AccAutoMin);
			}
			else{
				SetAcc(_Acc / AccAutoStep);
			}
		}
	}


	//収束判定を行う
	FlagConbergence = true;
	if( Error < MaxError ){
		FlagConbergence = false;
	}else{
		FlagConbergence = true;
	}


	//エラーディスプレイ
	if( FlagPrint == true ){
		Print();
	}
	if( FlagPrintAll == true ){
		PrintAll();
	}


	//収束していた場合は0を返す．
	if( FlagConbergence == true ){
		if( FlagEndReport == true ){
			EndReport();
		}
		return 0;
	}


	//MaxLoopを超えていたら抜け出す．
	if( MaxLoop < MainLoopCounter ){
		if( FlagEndReport == true ){
			EndReport();
		}
		//std::cout << "MaxLoop Over" << std::endl;
		return 0;
	}


	//修正値を入れるところ
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for( int i = 0 ; i < NumVariable ; i++ ){

		Step[i] = Acc[i] * DiagonalMatrix[i][ NumVariable ];

		if( MaxDisplacement[i] != 0.0 ){
			if( MaxDisplacement[i] < Step[i]  ){
				Step[i] = MaxDisplacement[i];
			}else if( Step[i] < -MaxDisplacement[i]){
				Step[i] = -MaxDisplacement[i];
			}
		}

		//あとで使うかもしれないから残しておく
		//if ((Variable[i][NumVariable] - Step[i])*(Variable[i][NumVariable]) < 0.0){
		//	JacobReset = true;
		//}
		//else if ( 0.9 < fabs(Step[i]/Variable[i][NumVariable]) ){
		//	JacobReset = true;
		//}
		//
		Variable[i][NumVariable] = Variable[i][NumVariable] - Step[i];

		if( MinValue[i] != 0.0 ){
			if( Variable[i][NumVariable] < MinValue[i] ){
				Variable[i][NumVariable] = MinValue[i];
			}
		}
		if( MaxValue[i] != 0.0 ){
			if( Variable[i][NumVariable] > MaxValue[i] ){
				Variable[i][NumVariable] = MaxValue[i];
			}
		}

		for( int k = 0 ; k < NumVariable ; k++ ){
			Variable[i][k] = Variable[i][NumVariable];
		}

		int j;
		j = i;
		//微小移動量の計算，仮定値が(Delta - 1.0)以下の時は常に（delta-1(Delta - 1.0)だけ動かすようにする
		//if (fabs(Variable[i][NumVariable])  < (Delta - 1.0)){
		//		if (0 < Variable[i][NumVariable]){
		//		DeltaValue[i] = (Delta - 1.0)*(Delta - 1.0);
		//		Variable[i][j] = Variable[i][NumVariable] + (Delta - 1.0)*(Delta - 1.0);
		//	}
		//	else{
		//		DeltaValue[i] = -(Delta - 1.0)*(Delta - 1.0);
		//		Variable[i][j] = Variable[i][NumVariable] - (Delta - 1.0)*(Delta - 1.0);
		//	}
		//}
		//微小移動量の計算，仮定値が1以下の時は常に（delta-1）だけ動かすようにする
		if (fabs(Variable[i][NumVariable])  < 1.0){
			if (0 < Variable[i][NumVariable]){
				DeltaValue[i] = (Delta - 1.0);
				Variable[i][j] = Variable[i][NumVariable] + (Delta - 1.0);
			}
			else{
				DeltaValue[i] = -(Delta - 1.0);
				Variable[i][j] = Variable[i][NumVariable] - (Delta - 1.0);
			}
		}
		else{
			DeltaValue[i] = Variable[i][NumVariable] * (Delta - 1.0);
			Variable[i][j] = Variable[i][NumVariable] * Delta;
		}


	}

	return 1;
}

int CNewtonRaphsonMethodPlus::MainLoopReinit(){
	double max;
	int maxi;

	//ヤコビ行列構築をスキップする場合はこのセクションをスキップ
	if( FlagJacobSkip == false ){

		//ヤコビ行列を作る
		//ここは並列化しない方が速い
		for( int i = 0 ; i < NumVariable ; i++ ){
			for( int j = 0 ; j < NumVariable ; j++ ){
				JacobMatrix[i][j] = ( ResultMatrix[i][j] - ResultMatrix[i][NumVariable] ) / DeltaValue[j];
			}
		}
	}

	//結果の行列は毎回書き換える
	for( int i = 0 ; i < NumVariable ; i++ ){
		JacobMatrix[i][NumVariable] = ResultMatrix[i][NumVariable];
	}


	//ガウスの消去法はそのうち関数化する
	//収束していない場合はここでガウスの消去法で修正値を求める．

	//JacobMatrix ==> UpperTriangularMatrix へのコピー
	//ここは並列化しない方が速い
	for (int i = 0; i < NumVariable; i++){
		for( int j = 0 ; j < NumVariable+1 ; j++ ){
			UpperTriangularMatrix[i][j] = JacobMatrix[i][j];
		}
	}


	//前進消去
	for( int k = 0 ; k < NumVariable ; k++ ){

		//部分ピボット選択
		if( UpperTriangularMatrix[k][k] < Min ){
	
			max = 0.0;
			for( int si = k ; si < NumVariable ; si++ ){
				if (1.0 < max){
					break;
				}
				if( max  < fabs( UpperTriangularMatrix[si][k] ) ){
					maxi = si;
					max  = fabs( UpperTriangularMatrix[si][k] );
				}
			}

			//ピボット選択エラー(0割り)
			if( max == 0.0 ){
				std::cout << "error : pivot(" << k << ")" << std::endl;
				std::cout << Name << "のNewtonRaphsonMethodPlusでピボット選択で0割が発生しました．プログラムを終了します" << std::endl;
				system("pause");
				exit(0);
//				return 0;
			}

			//入れ替え
			if( maxi != k ){
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for( int sj = k ; sj < NumVariable+1 ; sj++ ){
					double temp;
					temp = UpperTriangularMatrix[k][sj];
					UpperTriangularMatrix[k][sj] = UpperTriangularMatrix[maxi][sj];
					UpperTriangularMatrix[maxi][sj] = temp;
				}
			}
		}

		//前進消去部分
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for (int i = k + 1; i < NumVariable; i++){
			UpperC[i] = UpperTriangularMatrix[i][k] / UpperTriangularMatrix[k][k];
			if (fabs(UpperC[i]) > SparseMin ){
				for (int j = k; j < NumVariable + 1; j++){
					UpperTriangularMatrix[i][j] = UpperTriangularMatrix[i][j] - (UpperTriangularMatrix[k][j] * UpperC[i]);
				}
			}
		}		
	}

	//UpperTriangularMatrix ==> DiagonalMatrix へのコピー
	//ここは並列化しない方が速い
	for (int i = 0; i < NumVariable; i++){
		for( int j = 0 ; j < NumVariable+1 ; j++ ){
			DiagonalMatrix[i][j] = UpperTriangularMatrix[i][j];
		}
	}

	//後退代入
	for( int k = NumVariable-1 ; k >= 1 ; k-- ){
		for( int i = k-1 ; i >= 0 ; i-- ){
			DiagonalMatrix[i][NumVariable] -= DiagonalMatrix[k][NumVariable] * DiagonalMatrix[i][k] / DiagonalMatrix[k][k];
			DiagonalMatrix[i][k] = 0.0;
		}
	}

	//答え計算
	//ここは並列化しない方が速い
	for (int k = 0; k < NumVariable; k++){
		DiagonalMatrix[k][NumVariable] = DiagonalMatrix[k][NumVariable] / DiagonalMatrix[k][k];
		DiagonalMatrix[k][k] = 1.0;
	}
	

	//ここでカウンターを上げる．
	MainLoopCounter++;
	JacobSkipCounter++;

	//スキップ条件に合致すれば
	if (JacobStep == 1){
		FlagJacobSkip = false;
	}
	else if (JacobReset == true){
		FlagJacobSkip = false;
		JacobReset = false;
		JacobSkipCounter = 0;
	}
	else if (JacobStep <= JacobSkipCounter){
		FlagJacobSkip = false;
		JacobSkipCounter = 0;
	}
	else if (JacobSkipCounter < JacobStep){
		FlagJacobSkip = true;
	}


	if (FlagMatrixOut == true) {
		MatrixOut();
	}

	return 1;
}



//サブループはあんまりやることない
int CNewtonRaphsonMethodPlus::SubLoopInit(){
	if( FlagJacobSkip == true ){
		SubLoopCounter = NumVariable;//jacob行列の構築をスキップする場合はいきなり最後の計算からスタート
	}else{
		SubLoopCounter = 0;
	}
	return 1;
}
int CNewtonRaphsonMethodPlus::SubLoopCheck(){
	//ループカウンタを見て抜け出す
	if( NumVariable < SubLoopCounter ){
		return 0;
	}
	return 1;
}
int CNewtonRaphsonMethodPlus::SubLoopReinit(){
	//サブカウンターを上げる．
	SubLoopCounter++;
	return 1;
}


int CNewtonRaphsonMethodPlus::SetDelta( double _Delta ){
	Delta = _Delta;
	return 1;
}
int CNewtonRaphsonMethodPlus::SetError( double _Error ){
	Error = _Error;
	return 1;
}
int CNewtonRaphsonMethodPlus::SetMaxLoop( int _MaxLoop ){
	MaxLoop = _MaxLoop;
	return 1;
}
int CNewtonRaphsonMethodPlus::SetJacobStep( int _Step ){
	JacobStep = _Step;
	return 1;
}
int CNewtonRaphsonMethodPlus::EndReport( void ){
	time(&TimeNow);

	std::cout << "MaxError = " << MaxError << std::endl;
	std::cout << "Time = " << (difftime(TimeNow,TimeStart)) << std::endl;

	return 1;
}



CNewtonRaphsonMethodPlus::CNewtonRaphsonMethodPlus(){

	NumVariable = 0;
	FlagConbergence = false;
	FlagEndReport = false;
	FlagAccAuto = false;
	FlagMatrixOut = false;
	MaxLoop = 5;

	Variable.resize(1);
	Variable[0].resize(1);

	VariableMemory.resize(1);
	VariableMemory[0].resize(1);

	ResultMatrix.resize(1);
	ResultMatrix[0].resize(1);

	JacobMatrix.resize(1);
	JacobMatrix[0].resize(1);

	UpperTriangularMatrix.resize(1);
	UpperTriangularMatrix[0].resize(1);

	DiagonalMatrix.resize(1);
	DiagonalMatrix[0].resize(1);


	SetDelta( 1.00001 );
	SetError( 0.000001 );

	SetJacobStep( 1 );
	FlagJacobSkip = false;

	Min = 1.0;
	SparseMin = 1e-8;

	SetPrint( false , false );

	ErrorMax0 = 100001.0;
	ErrorMax1 = 100000.0;
	ErrorDif0 = -1.0;
	ErrorDif1 = -1.0;


}


CNewtonRaphsonMethodPlus::~CNewtonRaphsonMethodPlus(void){

}


int CNewtonRaphsonMethodPlus::StatusOut(){
	std::cout << "Variable[ " << Variable.size() << " ][ " << Variable[0].size() << " ]" << std::endl;
	std::cout << "VariableMemory[ " << VariableMemory.size() << " ][ " << VariableMemory[0].size() << " ]" << std::endl;
	std::cout << "JacobMatrix[ " << JacobMatrix.size() << " ][ " << JacobMatrix[0].size() << " ]" << std::endl;
	std::cout << "UpperTriangularMatrix[ " << UpperTriangularMatrix.size() << " ][ " << UpperTriangularMatrix[0].size() << " ]" << std::endl;
	std::cout << "DiagonalMatrix[ " << DiagonalMatrix.size() << " ][ " << DiagonalMatrix[0].size() << " ]" << std::endl;
	return 1;
}

int CNewtonRaphsonMethodPlus::GaussianElimination(){

	return 1;
}



int CNewtonRaphsonMethodPlus::Print(){
	std::cout << " Main Loop = " << MainLoopCounter << " Acc = " << Acc[0] << " MaxError = " << MaxError << " AveError = " << AveError << std::endl;
	return 1;
}

int CNewtonRaphsonMethodPlus::PrintAll(){
	int i;
	std::cout << " Main Loop = " << MainLoopCounter << std::endl;
	for( i = 0 ; i < NumVariable ; i++ ){
		std::cout << " Variable[" << i << "] = " << Variable[i][NumVariable] << " Error[" << i << "] = " << ResultMatrix[i][NumVariable] << std::endl;
	}
	return 1;
}

void CNewtonRaphsonMethodPlus::irekae(double &x,double &y ){double z=x;x=y;y=z;}
void CNewtonRaphsonMethodPlus::irekae(int    &x,int    &y ){int    z=x;x=y;y=z;}

void CNewtonRaphsonMethodPlus::NR_free(void){
	if(Variable.size()>0){
		std::vector<std::vector<double>>().swap(Variable);
	}
	if(VariableMemory.size()>0){
		std::vector<std::vector<double>>().swap(VariableMemory);
	}
	if(RelError.size()>0){
		std::vector<double>().swap(RelError); 
	}
	if(AbsError.size()>0){
		std::vector<double>().swap(AbsError); 
	}
	if(DeltaValue.size()>0){
	std::vector<double>().swap(DeltaValue); 
	}
	if(ResultMatrix.size()>0){
		std::vector<std::vector<double>>().swap(ResultMatrix);
	}
	if(JacobMatrix.size()>0){	
	std::vector<std::vector<double>>().swap(JacobMatrix);
	}
	if(UpperTriangularMatrix.size()>0){
		std::vector<std::vector<double>>().swap(UpperTriangularMatrix);
	}
	if(DiagonalMatrix.size()>0){
		std::vector<std::vector<double>>().swap(DiagonalMatrix);
	}
	if(UpperC.size()>0){
		std::vector<double>().swap(UpperC);
	}
	if(Acc.size()>0){
		std::vector<double>().swap(Acc);
	}
	if(MaxValue.size()>0){	
		std::vector<double>().swap(MaxValue);
	}
	if(MinValue.size()>0){
		std::vector<double>().swap(MinValue);
	}
	if(MaxDisplacement.size()){	
		std::vector<double>().swap(MaxDisplacement);
	}
	if(Step.size()){
		std::vector<double>().swap(Step);
	}
	if(StepControl.size()){
		std::vector<bool>().swap(StepControl);
	}
	if(MinMaxControl.size()){	
	std::vector<bool>().swap(MinMaxControl);
	}
//		std::cout<<"Newtonのデストラクタ"<<std::endl;

}
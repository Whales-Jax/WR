#ifndef __CNewtonRaphsonMethodPlus_H
#define __CNewtonRaphsonMethodPlus_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <vector>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif



#pragma warning( disable: 4996 )
#pragma warning( disable: 4819 )
#pragma warning( disable: 4018 )

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)






//Newton-Raphson法クラス
class CNewtonRaphsonMethodPlus
{
public:





	int Initialize();//初期化関数



	int SetValiable(int,double);//変数をセットする関数
	int SetValiable(int,double,double);//変数をセットする関数
	int SetValiable(int,double,double,double,double);//変数をセットする関数
	int SetValiableName(int, std::string);

	int SetMaxStep(int,double);//maxステップ幅を設定するモジュール
	int SetMinMaxValue(int,double,double);//最大値と最長値を設定するモジュール
	
	int GetValiable(int,double&);//変数を受け取る関数
	int GetValiableDirect(int,int,double&);//変数を受け取る関数

	int SetResult(int,double,double);//計算結果を入れる関数
	int SetResultDirect(int,int,double,double);//計算結果を入れる関数
	int SetResultMatrix( int , int , double );//計算結果を格納するマトリックスに直接値を入れる関数

	int SetError(double);//収束判定をセットする関数
	int SetDelta(double);//微小移動力をセットする関数

	int SetAcc(double);//加速係数をセットする関数
	int SetAccDirect(int,double);//加速係数を変数ごとにセットする関数
	int SetAccAuto(bool _FlagAccAuto, double _AccAutoMin, double _AccAutoMax, double _AccAutoStep);//加速係数を自動的に変更する

	int SetPrint( bool , bool );//errorの結果表示をするかどうかをセットする関数
	int SetEndReport( bool );//ニュートン法終了後にエラー値などの結果を出力するかどうかの設定
	int SetMaxLoop( int );//最大反復回数を入力する関数

	int SetName(std::string);//このニュートンラフソン法を識別するための名前を入力する関数

	int StatusOut();//行列の数などをアウトプットするデバッグ用の関数

	int SetMatrixOut(bool);
	int MatrixOut();//現在のマトリックスをアウトプットする

	int ValueContainerResize( int );

	int MainLoopInit(void);		//メインループ初期化式
	int MainLoopCheck(void);	//メインループ継続条件式
	int MainLoopReinit(void);	//メインループ再初期化式(メインとなる関数)
	int SubLoopInit(void);		//サブループ初期化式
	int SubLoopCheck(void);		//サブループ継続条件式
	int SubLoopReinit(void);	//サブループ再初期化式

	int GaussianElimination(void);

	int Print(void);
	int PrintAll(void);
	int EndReport(void);
	int SetJacobStep( int );

	void irekae(double &x,double &y);	//行列入れ替え用
	void irekae(int    &x,   int &y);

	void NR_free(void);//2014年8月中川原追加　ニュートン法の行列要素の解放
	//変数の数を数えるカウンタ
	//ループカウン


	std::vector<std::vector<double>> Variable;//修正値の配列
	std::vector<std::string> VariableName;//修正値の配列
	std::vector<std::vector<double>> VariableMemory;//修正値の履歴

	std::vector<double> RelError;//相対誤差
	std::vector<double> AbsError;//絶対誤差

	std::vector<double> DeltaValue;//移動量の絶対値


	std::vector<std::vector<double>> ResultMatrix;//計算結果の配列[行][列+1]
	std::vector<std::vector<double>> JacobMatrix;//ヤコビ行列の配列[行][列+1]
	std::vector<std::vector<double>> UpperTriangularMatrix;//上三角行列[行][列+1]
	std::vector<std::vector<double>> DiagonalMatrix;//対格行列[行][列+1]

	std::vector<double> UpperC;//上三角行列を作るときの係数

	std::vector<double> Acc;//加速係数
	std::vector<double> MaxValue;//変数の最大値
	std::vector<double> MinValue;//変数の最小値
	std::vector<double> MaxDisplacement;//変数の最大移動量
	std::vector<double> Step;//ステップ
	std::vector<bool> StepControl;//ステップ制限bool
	std::vector<bool> MinMaxControl;//MinMax制限bool



	int NumVariable;//変数の数=式の数

	int MainLoopCounter;//ニュートン法メインループのカウンタ
	int SubLoopCounter;//ニュートン法サブループのカウンタ

	int JacobSkipCounter;
	bool JacobReset;

	std::string Name;




	int MaxLoop;//最大反復回数

	int JacobStep;//


	double MaxError;
	double AveError;

	double Min;//
	double SparseMin;

	double keisu;
	double keisu2;


	double AccAcc;
	double AccMax;

	//AccAuto関係
	double AccAutoMax;
	double AccAutoMin;
	double AccAutoStep;
	bool FlagAccAuto;
	double ErrorMax0;
	double ErrorMax1;
	double ErrorDif0;
	double ErrorDif1;


	double Error;//収束判定
	double Delta;//数値的に偏微分するときの微小移動量

	bool FlagConbergence;//収束しているかどうかのフラグ.trueで収束している．
	bool FlagPrint;
	bool FlagPrintAll;
	bool FlagEndReport;
	bool FlagJacobSkip;
	bool FlagMatrixOut;//マトリックスをアウトプットするかどうかのフラグ

	time_t TimeStart;
	time_t TimeNow;


	std::fstream MatrixFile;

	CNewtonRaphsonMethodPlus();
	~CNewtonRaphsonMethodPlus();


};
#endif//CNewtonRaphsonMethod

#ifndef __CNewtonRaphsonMethod_H
#define __CNewtonRaphsonMethod_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>

#pragma warning( disable: 4996 )
//#include <boost/numeric/ublas/lu.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>


//ガウスの消去法（完全ピボッティング）を行うための関数．下記参照のこと．分からなくても問題ありません


//多変数ニュートン法メソッドフレームを構築するためのクラス
class CNewtonRaphsonMethod
{
public:
	std::string name;
	//変数
	int i_sys,i_error;		//ループ内で使用するカウンタ
	int count;				//計算回数
	int EscLoop;			//収束したとみなして抜け出す計算回数
	int count_max;			//最大計算回数			
	int count_div;			//半減周期
	int num;				//パラメータの数
	int _num;				//仮パラメータの数
	int ArrayNum;			//配列の数
	int Solver;				//連立方程式の解法 0:完全ピボット選択 1:部分ピボット選択
	int ErrorSolver;
	int ErrorLoopOver;
	bool FlagInitial;		//初期化の確認用
	bool FlagSetup;			//setupの確認用
	bool check;				//収束確認用
	bool check2;			//収束確認これが確認できたら収束
	bool Bprt;				//プリント表示用変数
	bool Bprt_sum;			//プリント表示用変数
	double delta;			//微小移動分の係数 x*delta が微小移動量になる
	double MinimumEPS;		//最小のdelta移動量
	double pibot;			//最小ピボット値
	double err;				//絶対誤差許容地
	double acc;				//緩和係数
	double err_sum;
	double err_ave;
	double err_max;
	double AccAcc;			//加速度勾配の加速度
	double AccMax;			//加速度勾配の最大値
	//配列
	int *ErrorMethod;			//エラー判定 0:相対誤差　1:絶対誤差
	double *Error;
	double *answer;				//現在地点と次の地点との差分
	double *set_value;			//繰り返し計算において一時的に使用する値（いわゆるT_vap_eva_enみたいの）
	double *residual;			//現在地点でのエラー値
	double *MaxDisplacement;	//最大移動量(max displacement)
	double *_MaxDisplacement;	//最大移動量の仮配列
	double *MaxValue;			//最大値(max value)
	double *MinValue;			//最小値(min value)
	double *assumed_value;		//現在地点の値をストック（いわゆるT_vap_eva_en_ntみたいの）
	double *_assumed_value;		//仮定値の仮配列
	double *step_value;			//微小移動地点での値のストック(T_vap_eva_en_nt)
	double **jacob;				//偏微分値（いわゆるJacob[NO][NO]）
	double **ErrorAbs;			//現在地点および微小移動した点でのエラー値(いわゆるError[NO][NO])絶対値
	double **ErrorRel;			//現在地点および微小移動した点でのエラー値(いわゆるError[NO][NO])相対値
	//matAvecX=vecBとしてvecXを求める
	//boost::numeric::ublas::matrix<double> matA;
	//boost::numeric::ublas::vector<double> vecX;
	//boost::numeric::ublas::vector<double> vecB;
	//関数
	CNewtonRaphsonMethod();//仮パラメータの数，一部メモリ確保
	~CNewtonRaphsonMethod(void);				//メモリ解放
	void setup(int n = 1000,double _delta=1.0+1e-6,double _err2=1e-6, double _pibot=-1.0);
	void reset(void);
	void setErrorMethod(int);
	void main_loop_init(void);			//メインループ初期化式
	bool main_loop_check(void);			//メインループ継続条件式
	void main_loop_reinit(void);		//メインループ再初期化式
	void sub_loop_init(void);			//サブループ初期化式
	bool sub_loop_check(void);			//サブループ継続条件式
	void sub_loop_reinit(void);			//サブループ再初期化式
	void initial(int n = 0);			//パラメータ数決定，メモリ確保
	void prt(double value=-1.0);		//エラー表示
	void prt_sum(void);					//全体エラー表示
	void prt2(double value=-1.0);		//エラー表示
	void prt_sum2(void);					//全体エラー表示
	void spcacc(double maxacc);				//acc自動調節

	double getValue(int i);				//修正値を返す関数

	void set_assumed_value(int i,double value,double _max = 0.0);//仮定値代入
	void setValue(int i , double value , double _max = 0.0 , double _MinValue = 0.0 , double _MaxValue = 0.0 );//仮定値代入
	void set_error(double,double x=1.0);//エラーセット
	void set_error2(int,double,double,double x=1.0);//エラーセット2
	void setError(int i ,double Value1 , double Value2 , int _ErrorMethod = 0 , double error = -1 );//エラーセット3
	void setAcc(double);
	void setSolver(int);
	void setEscLoop(int);
	void setDeltaError(double _delta=1.0+1e-6,double _err2=1e-6);


	std::fstream matrixFile;
	void matrixOut(void);
	void matrixFileMake(void);

	//ニュートン法用
private:
	int gaussian_elimination(void);		//完全ピボット選択（遅い
	int gaussian_elimination_2(void);	//部分ピボット選択（速い
//	int lu_factorize_ublas(void);		//lu分解　たぶん速い
	void irekae(double &x,double &y);
	void irekae(int    &x,   int &y);
};
#endif//CNewtonRaphsonMethod

#ifndef __CStructureParameter_h_
#define __CStructureParameter_h_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "CMaster222.h"

class CBasicParameter{
public:
	CBasicParameter(){
		PI = 4.0 * atan(1.0);
		g = 9.8;
	}
	double PI;
	double g;
};

class CStructureParameter 
{

public:
	CStructureParameter(){}
	~CStructureParameter(){}

	//CMAccumulator000
	double V; //!<体積[m3]

	//CMCompressor000
	//double V;//体積
	double rpm;//回転数
	double rps;//回転数
	double rps1;//回転数
	double rps2;//回転数
	double rps_per;//回転数(パーセント表示)
	double C;//膨張弁流量係数
	double C1;//膨張弁流量係数
	double C2;//膨張弁流量係数
	double C3;
	double C4;
	double C5;
	double C6;
	double S;//膨張弁断面積
	double W;//仕事
	double W1;//仕事
	double W2;//仕事
	double W3;//仕事
	double e_V;//体積効率
	double e_S;//断熱効率
	double e_M;//機械効率
	double e_W;//
	double e_INV;//インバーター効率
	double e_V1;//体積効率
	double e_S1;//断熱効率
	double e_M1;//機械効率
	double e_V2;//体積効率
	double e_S2;//断熱効率
	double e_M2;//機械効率
	double e_T1;//温度効率
	double e_T2;//温度効率

	double E;
	double E1;
	double E2;
	double E3;

	double torque1;//軸トルク
	double torque2;//軸トルク
	double torque3;//軸トルク

	double alpha1;
	double alpha2;
	double alpha3;
	double alpha1C;
	double alpha2C;
	double alpha3C;

	double cop;
	double cop1;
	double cop2;
	double cop3;

	double COR;

	int retsu;


	//CMExpansionValve000
	//double C;//膨張弁流量係数
	//double S;//断面積
	double open;//開度[%]
	double maxS;//最大断面積

	double pd_keisu;
	double ht_keisu;
	double ht_keisu_fin;
	double pd_keisu_fin;
	double ht_keisu_oil;
	double pd_keisu_oil;

	//CMHeatExchanger000
	double R1;
	double L1;
	double D1;
	double V1;
	double A1;
	double S1;
	double K1;
	double lambda1;
	double Arate1;
	double theta1;

	double R2;
	double L2;
	double D2;
	double V2;
	double A2;
	double S2;
	double K2;
	double lambda2;
	double Arate2;
	double theta2;

	double R3;
	double L3;
	double D3;
	double V3;
	double A3;
	double S3;
	double K3;
	double lambda3;
	double Arate3;
	double theta3;


	double R4;
	double L4;
	double D4;
	double V4;
	double A4;
	double S4;
	double K4;
	double lambda4;
	double Arate4;
	double theta4;


	double R5;
	double L5;
	double D5;
	double V5;
	double A5;
	double S5;
	double K5;
	double lambda5;
	double Arate5;
	double theta5;


	double X1;
	double Y1;
	double Z1;

	double X2;
	double Y2;
	double Z2;


	double X3;
	double Y3;
	double Z3;

	double X4;
	double Y4;
	double Z4;

	double X5;
	double Y5;
	double Z5;


	double finwariai;
	double D_M;
	double C_M;
	int In;//本数
	int In1;//本数
	int In2;//本数
	double L;//長さ

	int Idiv;
	double Mass_t0;
	double Q;
	double Q1;
	double Q2;
	double Q3;
	double Qtemp;
	double X;
	double KKK;
	double T;
	int DriveMode;//0:冷房　1:暖房
	double  DriveMode2;
	int CDriveMode;//0いつでも実行　1静特性のみ　2動特性のみ
	int Choke;

	double thc;

	int OnOffSwitch;
	
	double unit;
	double IncreasingBorder;
	double DecreasingBorder;
	double WaitTime;


	//フィンのパラメーター
	double fin_x;
	double fin_y;
	double fin_y_pitch;
	double fin_T;
	double fin_pitch;
	double fin_lambda;
	double fin_Ac;//自由通過面積
	double fin_Dec;//代表寸法
	double fin_Vac;//代表風速
	double e_FIN;
	double fin_retsu;//ふぃんの列

	//CMTube000
	//double D1;
	//double V1;
	//double A1;
	//double S1;
	//double K1;
	//double D2;
	//double V2;
	//double A2;
	//double S2;
	//double K2;
	//int In;//本数
	//double L;//長さ
	double outT;
	double outTset;
	double outTw;
	double aveT;
	double DaveT;
	double outG;//蛇口モデル用流量決定定数
	//int Idiv;
	//double D_M;
	//double C_M;
	//double Mass_t0;
	//double Q;
	double deg;//パイプの傾き　度
	double g;
	double H;//高さ
	double Hw;//液の高さ（液位）

	double KA;
	double KA_temp;
	double KA1;
	double KA2;
	double KA3;
	double KA4;

	double G;

	double n1;
	double n2;
	double G_temp1;
	double G_temp2;

	double nopt;
	double ng;
	double H_mfin;

	double beta1;
	double beta2;
	double beta3;


	double error;
	double delta;

	double **Tt;
	double **dTa;
	double **dTt;


	double a[100];
	double b[100];
	double c[100];
	int    n[100];


	//溝付き管　伝熱＆圧損　渡辺

	double groove_fin_Dn;//最大内径[m]
	double groove_fin_n;//フィン数
	double groove_fin_hf;//フィン高さ[m]
	double groove_fin_gamma;//フィンの頂角[°]
	double groove_fin_beta;//フィンのねじれ角[°]
	double alpha_in;
	//double groove_fin_pitch;//フィンピッチ[m]
	//double groove_fin_t;//フィン頂部幅[m]
	
	//渡辺終わり


	CMaster222 pq;
	CMaster222 ef;

	CMaster222 file01;
	CMaster222 file02;
	CMaster222 file03;
	CMaster222 file04;
	CMaster222 file05;

	char filename01[128];
	char filename02[128];
	char filename03[128];
	char filename04[128];
	char filename05[128];

	double SetValue;
	double *SetVariable;
	double *ProcessVariable;
	double *ManipulativeVariable;

	char cSetVariable[128];
	char cProcessVariable[128];
	char cManipulativeVariable[128];

	double *ProcessVariable_A;
	double *ProcessVariable_B;
	double *ProcessVariable_C;

	double *ProcessVariable_Step;
	double *ProcessVariable_PID;
	
	double *Constant;

	char cFolderName[128];
	char cFileName[128];



	int SetModule;
	int ProcessModule;
	int ManipulativeModule;

	char cSetModule[128];
	char cProcessModule[128];
	char cManipulativeModule[128];

	double ManipulativeVariableIni;


	int IndNum;
	int GenNum;
	double SurvivalRate;
	double MutationRate;
	int MinOrMax;
	int SimpleOrMulti;

	double Evaluation;


	int		ObjMinOrMax_00;
	double	*ProcessObjectiveVariable_00;
	char		cProcessObjectiveVariable_00[128];
	char		cProcessObjectiveModule_00[128];
	double	Weight_Param_00;

	int		ObjMinOrMax_01;
	double	*ProcessObjectiveVariable_01;
	char		cProcessObjectiveVariable_01[128];
	char		cProcessObjectiveModule_01[128];
	double	Weight_Param_01;

	int		ObjMinOrMax_02;
	double	*ProcessObjectiveVariable_02;
	char		cProcessObjectiveVariable_02[128];
	char		cProcessObjectiveModule_02[128];
	double	Weight_Param_02;

	char		cProcessConstraintModule_00[128];
	char		cProcessConstraintVariable_00[128];
	double	*ProcessConstraintVariable_00;
	double	Constraint_number_00;
	double	Tuning_alpha_00;

	char		cProcessConstraintModule_01[128];
	char		cProcessConstraintVariable_01[128];
	double	*ProcessConstraintVariable_01;
	double	Constraint_number_01;
	double	Tuning_alpha_01;

	char		cProcessConstraintModule_02[128];
	char		cProcessConstraintVariable_02[128];
	double	*ProcessConstraintVariable_02;
	double	Constraint_number_02;
	double	Tuning_alpha_02;

	char		cProcessConstraintModule_03[128];
	char		cProcessConstraintVariable_03[128];
	double	*ProcessConstraintVariable_03;
	double	Constraint_number_03;
	double	Tuning_alpha_03;


	int     ManipulativeModule_00;
	double *ManipulativeVariable_00;
	char   cManipulativeVariable_00[128];
	char   cManipulativeModule_00[128];
	double ManipulativeVariableMin_00;
	double ManipulativeVariableMax_00;

	int     ManipulativeModule_01;
	double *ManipulativeVariable_01;
	char   cManipulativeVariable_01[128];
	char   cManipulativeModule_01[128];
	double ManipulativeVariableMin_01;
	double ManipulativeVariableMax_01;

	int     ManipulativeModule_02;
	double *ManipulativeVariable_02;
	char   cManipulativeVariable_02[128];
	char   cManipulativeModule_02[128];
	double ManipulativeVariableMin_02;
	double ManipulativeVariableMax_02;

	int     ManipulativeModule_03;
	double *ManipulativeVariable_03;
	char   cManipulativeVariable_03[128];
	char   cManipulativeModule_03[128];
	double ManipulativeVariableMin_03;
	double ManipulativeVariableMax_03;

	int     ManipulativeModule_04;
	double *ManipulativeVariable_04;
	char   cManipulativeVariable_04[128];
	char   cManipulativeModule_04[128];
	double ManipulativeVariableMin_04;
	double ManipulativeVariableMax_04;

	int     ManipulativeModule_05;
	double *ManipulativeVariable_05;
	char   cManipulativeVariable_05[128];
	char   cManipulativeModule_05[128];
	double ManipulativeVariableMin_05;
	double ManipulativeVariableMax_05;

	int     ManipulativeModule_06;
	double *ManipulativeVariable_06;
	char   cManipulativeVariable_06[128];
	char   cManipulativeModule_06[128];
	double ManipulativeVariableMin_06;
	double ManipulativeVariableMax_06;

	int     ManipulativeModule_07;
	double *ManipulativeVariable_07;
	char   cManipulativeVariable_07[128];
	char   cManipulativeModule_07[128];
	double ManipulativeVariableMin_07;
	double ManipulativeVariableMax_07;

	int     ManipulativeModule_08;
	double *ManipulativeVariable_08;
	char   cManipulativeVariable_08[128];
	char   cManipulativeModule_08[128];
	double ManipulativeVariableMin_08;
	double ManipulativeVariableMax_08;

	int     ManipulativeModule_09;
	double *ManipulativeVariable_09;
	char   cManipulativeVariable_09[128];
	char   cManipulativeModule_09[128];
	double ManipulativeVariableMin_09;
	double ManipulativeVariableMax_09;

	int     ManipulativeModule_10;
	double *ManipulativeVariable_10;
	char   cManipulativeVariable_10[128];
	char   cManipulativeModule_10[128];
	double ManipulativeVariableMin_10;
	double ManipulativeVariableMax_10;

	int     ManipulativeModule_11;
	double *ManipulativeVariable_11;
	char   cManipulativeVariable_11[128];
	char   cManipulativeModule_11[128];
	double ManipulativeVariableMin_11;
	double ManipulativeVariableMax_11;

	int    PID_Type;
	double ProportionalGain;
	double IntegralGain;
	double DerivativeGain;
	
	double ManipulativeValueMax;
	double ManipulativeValueMin;

	double ManipulativeValueA;
	double ManipulativeValueB;
	double SetValueAtoB;
	double SetValueBtoA;

	double ValueA;
	double ValueB;
	double ValueC;
	double Value;


	//A+B-Cとかに使う物
	char cProcessVariableA[128];
	char cProcessVariableB[128];
	char cProcessVariableC[128];
	char cProcessModuleA[128];
	char cProcessModuleB[128];
	char cProcessModuleC[128];
	double *ProcessVariableA;
	double *ProcessVariableB;
	double *ProcessVariableC;
	double CalcResult;

	double dummyA;
	double dummyB;
	double dummyC;

	//配列による制御変数
	char cPV[100][128];
	char cPM[100][128];
	char cMV[100][128];
	char cMM[100][128];
	std::vector<double*> PV;		//ProcessVariable
	std::vector<double*> MV;		//ManipulativeVariable
	std::vector<std::vector<double>> VValue;

	//CTimeに必要なもの
	char ManipulativeFile[128];

	double Time;

	bool   check;
	int    flag;
	int    counter;
	double timer;

	//inputとoutputを交換するモジュール
	double *InputToOutputVariable;
	char cInputToOutputModule[128];
	char cInputToOutputVariable[128];

	double *OutputToInputVariable;
	char cOutputToInputModule[128];
	char cOutputToInputVariable[128];

	//sin input と square input
	double amplitude;
	double Aoffset;
	double wavelength;
	double offset;


	double delta_y;
	double delta_x;
	double _delta_x;
	double integral_y;
	double _integral_y;
	double derivative_y;


	double SHF;

	int ImplicitSwitch;//0:Explicit, 1:Implicit, 2:SemiImplicit


	//efm2fluentなどで仕様するものｔ
	fstream FileOut01;
	char ProcessFile[128];
	char FolderPath[128];
	char cProcessID[128];
	char cProcessID_02[128];
	char cController[128];
	char cStatus[128];
	bool ReadFlug;

};

#endif
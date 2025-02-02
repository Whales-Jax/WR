#ifndef __CCellModule_h
#define __CCellModule_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <vector>

using namespace std;

#include "DXProperty_ver06.h"
#include "Moist_Air_Property.h"
#include "CNewtonRaphsonMethod.h"
#include "CStructureParameter.h"
#include "CBisectionMethod.h"
#include "CFluidParameter.h"
#include "CPressureDropAndHeatTransfer.h"
#include "CSuperModule.h"
#include "CMaster222.h"
#include "PropertyWater_ver05.h"
#include "PropertyLithiumBromide_ver05.h"
#include "property_air_ver2_2.h"
#include "PropertyNaClaq_ver01.h"

#define SENGEN_PROTOTYPE( NUM ) void CM##NUM##init(); void CM##NUM##setup(); void CM##NUM##calc(int); void CM##NUM##outletCheck(int); void CM##NUM##renew();

class CCycleSystemModule;

class CCellModule :public CSuperModule {
public:
	CCellModule(){//CCellModuleのコンストラクタ

		ref.reibai = "R32.fld,R125.fld";
		ref.joutai = "IIR";
		ref.massFraction[0] = 0.5;
		ref.massFraction[1] = 0.5;

		Htransfer = 0;
		Pdrop = 0;
		ModuleType = 0;

		bis.setupError( 1e-10 );
		dt = 0.0;

		P.OnOffSwitch = 1;

		keisan_finish = 0;
		portnum = 20;
		portnum_inuse = portnum;
		nVariable = -1;

		Po = new CFluidParameter[portnum+1];
		_Po = new CFluidParameter[portnum+1];
//		MaxStep.resize(portnum);
//		Value.resize(portnum);
//		MinValue.resize(portnum);
//		MaxValue.resize(portnum);
		MaxStep = new double[portnum+1];
		Value = new double*[portnum+1];
		MaxValue = new double[portnum+1];
		MinValue = new double[portnum+1];

		for( int i = 0 ; i <= portnum ; i++ ){
			MaxStep[i] = 0.0;
			MaxValue[i] = 0.0;
			MinValue[i] = 0.0;
		}

		for( int i = 0 ; i <= portnum ; i++ ){
			Po[i].IONum = i;
			Po[i].ConnectModule = -1;//モジュールナンバーは0から振るので，初期値として-1（接続されていない）を入れておく
			Po[i].ConnectPort = 0;//コネクトポートは1から振るので，初期値として0（接続されていない）を入れておく
			Po[i].IOCode = 2;
			Po[i].FluidCode = NULL;
		}

		X_t1 = new CFluidParameter[portnum];
		_X_t1 = new CFluidParameter[portnum];
		X_t0 = new CFluidParameter[portnum];
		_X_t0 = new CFluidParameter[portnum];

		Y_t1 = new CFluidParameter[portnum];
		_Y_t1 = new CFluidParameter[portnum];
		Y_t0 = new CFluidParameter[portnum];
		_Y_t0 = new CFluidParameter[portnum];

		Z_t1 = new CFluidParameter[portnum];
		_Z_t1 = new CFluidParameter[portnum];
		Z_t0 = new CFluidParameter[portnum];
		_Z_t0 = new CFluidParameter[portnum];


		//充填量エラーの初期化
		CycleMAX = 100;
		CycleNUM = 0;
		MassCharge = new double[CycleMAX];
		for( int i = 0 ; i < CycleMAX ; i++ ){
			MassCharge[i] = 0.0;
		}
		ErrorOfContinuity = -1;

		P.ImplicitSwitch = 1;

	}
	~CCellModule(){
	
		delete[] X_t1;
		delete[] _X_t1;
		delete[] X_t0;
		delete[] _X_t0;
	
		delete[] Y_t1;
		delete[] _Y_t1;
		delete[] Y_t0;
		delete[] _Y_t0;

		delete[] Z_t1;
		delete[] _Z_t1;
		delete[] Z_t0;
		delete[] _Z_t0;
	
		delete[] MaxStep;
		delete[] Value;
		delete[] MaxValue;
		delete[] MinValue;

		delete[] Po;
		delete[] _Po;

	}



	//宣言
	int flag1;
	int flag2;

	int nVariable;
	int keisan_finish;
	int ModuleType;

	int xDIV;
	int tDIV;

	int portnum;
	int portnum_inuse;
	int Iportnum;
	int Oportnum;
	CFluidParameter *Po;
	CFluidParameter *_Po;

	//充填量エラー
	int CycleNUM;
	int CycleMAX;
	double *MassCharge;
	int ErrorOfContinuity;

	CCycleSystemModule *SystemCopy;

	CCellModule *CellCopy;
	int divCopy;
	CCellModule *C_CellCopy;
	int C_divCopy;

	bool TopCell;
	bool EndCell;

//	vector <double> MaxStep;
//	vector <double*> Value;
//	vector <double> MinValue;
//	vector <double> MaxValue;
	double *MaxStep;
	double **Value;
	double *MinValue;
	double *MaxValue;


	CFluidParameter R_t1;
	CFluidParameter _R_t1;
	CFluidParameter R_t0;

	CFluidParameter V_t1;
	CFluidParameter _V_t1;
	CFluidParameter V_t0;

	CFluidParameter M_t1;
	CFluidParameter _M_t1;
	CFluidParameter M_t0;

	CFluidParameter A_t1;
	CFluidParameter _A_t1;
	CFluidParameter A_t0;

	CFluidParameter B_t1;
	CFluidParameter _B_t1;
	CFluidParameter B_t0;

	CFluidParameter C_t1;
	CFluidParameter _C_t1;
	CFluidParameter C_t0;

	CFluidParameter *W_t1;
	CFluidParameter *_W_t1;
	CFluidParameter *W_t0;
	CFluidParameter *_W_t0;

	CFluidParameter *X_t1;
	CFluidParameter *_X_t1;
	CFluidParameter *X_t0;
	CFluidParameter *_X_t0;

	CFluidParameter *Y_t1;
	CFluidParameter *_Y_t1;
	CFluidParameter *Y_t0;
	CFluidParameter *_Y_t0;

	CFluidParameter *Z_t1;
	CFluidParameter *_Z_t1;
	CFluidParameter *Z_t0;
	CFluidParameter *_Z_t0;

	DXProperty ref;
	DXProperty airref;
	DXProperty wat;
	Air air;
	PropertyWater wat2;
	PropertyLithiumBromide libr;
	CPressureDropAndHeatTransfer pdht;
//	PropertyAir P_air;
	PropertyNaClaq naclaq;

	CNewtonRaphsonMethod mnm;
	CNewtonRaphsonMethod mnm1;
	CNewtonRaphsonMethod mnm2;
	CBisectionMethod bis;

	CMaster222 gaibu_input;
	char FolderPath[128];


	//この5つの関数はモジュールによって処理が違う
	void init();
	void setup();
	void calc(int);
	void outletCheck(int);
	void renew();
	void SearchPointer( char *_Module , char *_Variable , double* (&Pointer) );

	//この4つの関数はモジュールに限らず処理が同じ だから関数ポインタを作る必要なし
	void getVariable(int,double&,double&,double&,double&);
	void setVariable(int,double);
	void fileout(fstream&);
	void nameout(fstream&);
	void setPointer();
	void setVariablePointer();

	//処理が違うモジュールの関数は関数ポインタで見た目は一括処理
	//モジュールごとにCPPで分割し関数ポインタを使う
	void (CCellModule::*init2)();
	void (CCellModule::*setup2)();
	void (CCellModule::*calc2)(int);
	void (CCellModule::*outletCheck2)(int);
	void (CCellModule::*renew2)();




	//関数の配列たち
	void (CCellModule::*InitPtr[10000])();
	void (CCellModule::*SetupPtr[10000])();
	void (CCellModule::*CalcPtr[10000])(int);
	void (CCellModule::*OCheckPtr[10000])(int);
	void (CCellModule::*RenewPtr[10000])();




	//init()やcalc()関数のプロトタイプ宣言
	

	//0000番台は大規模次世代ヒーポン
	void CM0001init();void CM0001setup();void CM0001calc(int);void CM0001outletCheck(int);void CM0001renew();
	void CM0002init();void CM0002setup();void CM0002calc(int);void CM0002outletCheck(int);void CM0002renew();
	void CM0003init();void CM0003setup();void CM0003calc(int);void CM0003outletCheck(int);void CM0003renew();
	void CM0004init();void CM0004setup();void CM0004calc(int);void CM0004outletCheck(int);void CM0004renew();
	void CM0005init();void CM0005setup();void CM0005calc(int);void CM0005outletCheck(int);void CM0005renew();
	void CM0006init();void CM0006setup();void CM0006calc(int);void CM0006outletCheck(int);void CM0006renew();
	void CM0007init();void CM0007setup();void CM0007calc(int);void CM0007outletCheck(int);void CM0007renew();
	void CM0008init();void CM0008setup();void CM0008calc(int);void CM0008outletCheck(int);void CM0008renew();
	void CM0009init();void CM0009setup();void CM0009calc(int);void CM0009outletCheck(int);void CM0009renew();
	void CM0010init();void CM0010setup();void CM0010calc(int);void CM0010outletCheck(int);void CM0010renew();
	void CM0011init();void CM0011setup();void CM0011calc(int);void CM0011outletCheck(int);void CM0011renew();
	void CM0012init();void CM0012setup();void CM0012calc(int);void CM0012outletCheck(int);void CM0012renew();
	void CM0013init();void CM0013setup();void CM0013calc(int);void CM0013outletCheck(int);void CM0013renew();
	void CM0014init();void CM0014setup();void CM0014calc(int);void CM0014outletCheck(int);void CM0014renew();
	void CM0015init();void CM0015setup();void CM0015calc(int);void CM0015outletCheck(int);void CM0015renew();
	void CM0016init();void CM0016setup();void CM0016calc(int);void CM0016outletCheck(int);void CM0016renew();
	void CM0017init();void CM0017setup();void CM0017calc(int);void CM0017outletCheck(int);void CM0017renew();
	void CM0018init();void CM0018setup();void CM0018calc(int);void CM0018outletCheck(int);void CM0018renew();
	void CM0019init();void CM0019setup();void CM0019calc(int);void CM0019outletCheck(int);void CM0019renew();
	void CM0020init();void CM0020setup();void CM0020calc(int);void CM0020outletCheck(int);void CM0020renew();
	void CM0021init();void CM0021setup();void CM0021calc(int);void CM0021outletCheck(int);void CM0021renew();
	void CM0022init();void CM0022setup();void CM0022calc(int);void CM0022outletCheck(int);void CM0022renew();
	void CM0023init();void CM0023setup();void CM0023calc(int);void CM0023outletCheck(int);void CM0023renew();
	void CM0024init();void CM0024setup();void CM0024calc(int);void CM0024outletCheck(int);void CM0024renew();
	void CM0025init();void CM0025setup();void CM0025calc(int);void CM0025outletCheck(int);void CM0025renew();
	void CM0026init();void CM0026setup();void CM0026calc(int);void CM0026outletCheck(int);void CM0026renew();
	void CM0027init();void CM0027setup();void CM0027calc(int);void CM0027outletCheck(int);void CM0027renew();
	void CM0028init();void CM0028setup();void CM0028calc(int);void CM0028outletCheck(int);void CM0028renew();
	void CM0029init();void CM0029setup();void CM0029calc(int);void CM0029outletCheck(int);void CM0029renew();
	void CM0030init();void CM0030setup();void CM0030calc(int);void CM0030outletCheck(int);void CM0030renew();
	void CM0031init();void CM0031setup();void CM0031calc(int);void CM0031outletCheck(int);void CM0031renew();
	void CM0032init();void CM0032setup();void CM0032calc(int);void CM0032outletCheck(int);void CM0032renew();
	void CM0033init();void CM0033setup();void CM0033calc(int);void CM0033outletCheck(int);void CM0033renew();
	void CM0034init();void CM0034setup();void CM0034calc(int);void CM0034outletCheck(int);void CM0034renew();
	void CM0035init();void CM0035setup();void CM0035calc(int);void CM0035outletCheck(int);void CM0035renew();
	void CM0036init();void CM0036setup();void CM0036calc(int);void CM0036outletCheck(int);void CM0036renew();
	void CM0037init();void CM0037setup();void CM0037calc(int);void CM0037outletCheck(int);void CM0037renew();
	void CM0038init();void CM0038setup();void CM0038calc(int);void CM0038outletCheck(int);void CM0038renew();
	void CM0039init();void CM0039setup();void CM0039calc(int);void CM0039outletCheck(int);void CM0039renew();
	void CM0040init();void CM0040setup();void CM0040calc(int);void CM0040outletCheck(int);void CM0040renew();
	void CM0041init();void CM0041setup();void CM0041calc(int);void CM0041outletCheck(int);void CM0041renew();
	void CM0042init();void CM0042setup();void CM0042calc(int);void CM0042outletCheck(int);void CM0042renew();
	void CM0043init();void CM0043setup();void CM0043calc(int);void CM0043outletCheck(int);void CM0043renew();
	void CM0044init();void CM0044setup();void CM0044calc(int);void CM0044outletCheck(int);void CM0044renew();
	void CM0045init();void CM0045setup();void CM0045calc(int);void CM0045outletCheck(int);void CM0045renew();
	void CM0046init();void CM0046setup();void CM0046calc(int);void CM0046outletCheck(int);void CM0046renew();
	void CM0047init();void CM0047setup();void CM0047calc(int);void CM0047outletCheck(int);void CM0047renew();
	void CM0048init();void CM0048setup();void CM0048calc(int);void CM0048outletCheck(int);void CM0048renew();
	void CM0049init();void CM0049setup();void CM0049calc(int);void CM0049outletCheck(int);void CM0049renew();

	//1000番台は吸収式シリーズ
	SENGEN_PROTOTYPE( 1001 );
	SENGEN_PROTOTYPE( 1002 );
	SENGEN_PROTOTYPE( 1003 );
	SENGEN_PROTOTYPE( 1004 );
	SENGEN_PROTOTYPE( 1005 );
	SENGEN_PROTOTYPE( 1006 );
	SENGEN_PROTOTYPE( 1007 );
	SENGEN_PROTOTYPE( 1008 );
	SENGEN_PROTOTYPE( 1009 );
	SENGEN_PROTOTYPE( 1010 );
	SENGEN_PROTOTYPE( 1011 );
	SENGEN_PROTOTYPE( 1012 );
	SENGEN_PROTOTYPE( 1013 );
	SENGEN_PROTOTYPE( 1014 );
	SENGEN_PROTOTYPE( 1015 );
	SENGEN_PROTOTYPE( 1016 );
	SENGEN_PROTOTYPE( 1017 );
	SENGEN_PROTOTYPE( 1018 );
	SENGEN_PROTOTYPE( 1019 );
	SENGEN_PROTOTYPE( 1020 );
	SENGEN_PROTOTYPE( 1021 );
	SENGEN_PROTOTYPE( 1022 );
	SENGEN_PROTOTYPE( 1023 );
	SENGEN_PROTOTYPE( 1024 );
	SENGEN_PROTOTYPE( 1025 );
	SENGEN_PROTOTYPE( 1026 );
	SENGEN_PROTOTYPE( 1027 );
	SENGEN_PROTOTYPE( 1028 );
	SENGEN_PROTOTYPE( 1029 );
	SENGEN_PROTOTYPE( 1030 );
	SENGEN_PROTOTYPE( 1031 );
	SENGEN_PROTOTYPE( 1032 );
	SENGEN_PROTOTYPE( 1033 );
	SENGEN_PROTOTYPE( 1034 );
	SENGEN_PROTOTYPE( 1035 );
	SENGEN_PROTOTYPE( 1036 );
	SENGEN_PROTOTYPE( 1037 );
	SENGEN_PROTOTYPE( 1038 );
	SENGEN_PROTOTYPE( 1039 );
	SENGEN_PROTOTYPE( 1040 );
	SENGEN_PROTOTYPE( 1041 );
	SENGEN_PROTOTYPE( 1042 );
	SENGEN_PROTOTYPE( 1043 );
	SENGEN_PROTOTYPE( 1044 );
	SENGEN_PROTOTYPE( 1045 );
	SENGEN_PROTOTYPE( 1046 );
	SENGEN_PROTOTYPE( 1047 );
	SENGEN_PROTOTYPE( 1048 );
	SENGEN_PROTOTYPE( 1049 );
	SENGEN_PROTOTYPE( 1050 );
	SENGEN_PROTOTYPE( 1051 );
	SENGEN_PROTOTYPE( 1052 );
	SENGEN_PROTOTYPE( 1053 );
	SENGEN_PROTOTYPE( 1054 );
	SENGEN_PROTOTYPE( 1055 );
	SENGEN_PROTOTYPE( 1056 );
	SENGEN_PROTOTYPE( 1057 );
	SENGEN_PROTOTYPE( 1058 );
	SENGEN_PROTOTYPE( 1059 );
	SENGEN_PROTOTYPE( 1060 );
	SENGEN_PROTOTYPE( 1061 );
	SENGEN_PROTOTYPE( 1062 );
	SENGEN_PROTOTYPE( 1063 );
	SENGEN_PROTOTYPE( 1064 );
	SENGEN_PROTOTYPE( 1065 );
	SENGEN_PROTOTYPE( 1066 );
	SENGEN_PROTOTYPE( 1067 );
	SENGEN_PROTOTYPE( 1068 );
	SENGEN_PROTOTYPE( 1069 );
	SENGEN_PROTOTYPE( 1070 );
	SENGEN_PROTOTYPE( 1071 );
	SENGEN_PROTOTYPE( 1072 );
	SENGEN_PROTOTYPE( 1073 );
	SENGEN_PROTOTYPE( 1074 );
	SENGEN_PROTOTYPE( 1075 );
	SENGEN_PROTOTYPE( 1076 );
	SENGEN_PROTOTYPE( 1077 );
	SENGEN_PROTOTYPE( 1078 );
	SENGEN_PROTOTYPE( 1079 );
	SENGEN_PROTOTYPE( 1080 );
	SENGEN_PROTOTYPE( 1081 );
	SENGEN_PROTOTYPE( 1082 );
	SENGEN_PROTOTYPE( 1083 );
	SENGEN_PROTOTYPE( 1084 );
	SENGEN_PROTOTYPE( 1085 );
	SENGEN_PROTOTYPE( 1086 );
	SENGEN_PROTOTYPE( 1087 );
	SENGEN_PROTOTYPE( 1088 );
	SENGEN_PROTOTYPE( 1089 );
	SENGEN_PROTOTYPE( 1090 );
	SENGEN_PROTOTYPE( 1091 );
	SENGEN_PROTOTYPE( 1092 );
	SENGEN_PROTOTYPE( 1093 );
	SENGEN_PROTOTYPE( 1094 );
	SENGEN_PROTOTYPE( 1095 );
	SENGEN_PROTOTYPE( 1096 );
	SENGEN_PROTOTYPE( 1097 );
	SENGEN_PROTOTYPE( 1098 );
	SENGEN_PROTOTYPE( 1099 );
	SENGEN_PROTOTYPE( 1100 );
	SENGEN_PROTOTYPE( 1101 );
	SENGEN_PROTOTYPE( 1102 );
	SENGEN_PROTOTYPE( 1103 );
	SENGEN_PROTOTYPE( 1104 );
	SENGEN_PROTOTYPE( 1105 );
	SENGEN_PROTOTYPE( 1106 );
	SENGEN_PROTOTYPE( 1107 );
	SENGEN_PROTOTYPE( 1108 );
	SENGEN_PROTOTYPE( 1109 );
	SENGEN_PROTOTYPE( 1110 );
	SENGEN_PROTOTYPE( 1111 );
	SENGEN_PROTOTYPE( 1112 );
	SENGEN_PROTOTYPE( 1113 );
	SENGEN_PROTOTYPE( 1114 );
	SENGEN_PROTOTYPE( 1115 );
	SENGEN_PROTOTYPE( 1116 );
	SENGEN_PROTOTYPE( 1117 );
	SENGEN_PROTOTYPE( 1118 );
	SENGEN_PROTOTYPE( 1119 );
	SENGEN_PROTOTYPE( 1120 );
	SENGEN_PROTOTYPE( 1121 );
	SENGEN_PROTOTYPE( 1122 );
	SENGEN_PROTOTYPE( 1123 );
	SENGEN_PROTOTYPE( 1124 );
	SENGEN_PROTOTYPE( 1125 );
	SENGEN_PROTOTYPE( 1126 );
	SENGEN_PROTOTYPE( 1127 );
	SENGEN_PROTOTYPE( 1128 );
	SENGEN_PROTOTYPE( 1129 );
	SENGEN_PROTOTYPE( 1130 );
	SENGEN_PROTOTYPE( 1131 );
	SENGEN_PROTOTYPE( 1132 );
	SENGEN_PROTOTYPE( 1133 );
	SENGEN_PROTOTYPE( 1134 );
	SENGEN_PROTOTYPE( 1135 );
	SENGEN_PROTOTYPE( 1136 );
	SENGEN_PROTOTYPE( 1137 );
	SENGEN_PROTOTYPE( 1138 );
	SENGEN_PROTOTYPE( 1139 );
	SENGEN_PROTOTYPE( 1140 );
	SENGEN_PROTOTYPE( 1141 );
	SENGEN_PROTOTYPE( 1142 );
	SENGEN_PROTOTYPE( 1143 );
	SENGEN_PROTOTYPE( 1144 );
	SENGEN_PROTOTYPE( 1145 );
	SENGEN_PROTOTYPE( 1146 );
	SENGEN_PROTOTYPE( 1147 );
	SENGEN_PROTOTYPE( 1148 );
	SENGEN_PROTOTYPE( 1149 );
	SENGEN_PROTOTYPE( 1150 );
	SENGEN_PROTOTYPE( 1151 );
	SENGEN_PROTOTYPE( 1152 );
	SENGEN_PROTOTYPE( 1153 );
	SENGEN_PROTOTYPE( 1154 );
	SENGEN_PROTOTYPE( 1155 );
	SENGEN_PROTOTYPE( 1156 );
	SENGEN_PROTOTYPE( 1157 );
	SENGEN_PROTOTYPE( 1158 );
	SENGEN_PROTOTYPE( 1159 );
	SENGEN_PROTOTYPE( 1160 );
	SENGEN_PROTOTYPE( 1161 );
	SENGEN_PROTOTYPE( 1162 );
	SENGEN_PROTOTYPE( 1163 );
	SENGEN_PROTOTYPE( 1164 );
	SENGEN_PROTOTYPE( 1165 );
	SENGEN_PROTOTYPE( 1166 );
	SENGEN_PROTOTYPE( 1167 );
	SENGEN_PROTOTYPE( 1168 );
	SENGEN_PROTOTYPE( 1169 );
	SENGEN_PROTOTYPE( 1170 );
	SENGEN_PROTOTYPE( 1171 );
	SENGEN_PROTOTYPE( 1172 );
	SENGEN_PROTOTYPE( 1173 );
	SENGEN_PROTOTYPE( 1174 );
	SENGEN_PROTOTYPE( 1175 );
	SENGEN_PROTOTYPE( 1176 );
	SENGEN_PROTOTYPE( 1177 );
	SENGEN_PROTOTYPE( 1178 );
	SENGEN_PROTOTYPE( 1179 );
	SENGEN_PROTOTYPE( 1180 );
	SENGEN_PROTOTYPE( 1181 );
	SENGEN_PROTOTYPE( 1182 );
	SENGEN_PROTOTYPE( 1183 );
	SENGEN_PROTOTYPE( 1184 );
	SENGEN_PROTOTYPE( 1185 );
	SENGEN_PROTOTYPE( 1186 );
	SENGEN_PROTOTYPE( 1187 );
	SENGEN_PROTOTYPE( 1188 );
	SENGEN_PROTOTYPE( 1189 );
	SENGEN_PROTOTYPE( 1190 );
	SENGEN_PROTOTYPE( 1191 );
	SENGEN_PROTOTYPE( 1192 );
	SENGEN_PROTOTYPE( 1193 );
	SENGEN_PROTOTYPE( 1194 );
	SENGEN_PROTOTYPE( 1195 );
	SENGEN_PROTOTYPE( 1196 );
	SENGEN_PROTOTYPE( 1197 );
	SENGEN_PROTOTYPE( 1198 );
	SENGEN_PROTOTYPE( 1199 );
	SENGEN_PROTOTYPE( 1200 );


	//2000番台は圧縮式シリーズ
	SENGEN_PROTOTYPE( 2001 );
	SENGEN_PROTOTYPE( 2002 );
	SENGEN_PROTOTYPE( 2003 );
	SENGEN_PROTOTYPE( 2004 );
	SENGEN_PROTOTYPE( 2005 );
	SENGEN_PROTOTYPE( 2006 );
	SENGEN_PROTOTYPE( 2007 );
	SENGEN_PROTOTYPE( 2008 );
	SENGEN_PROTOTYPE( 2009 );
	SENGEN_PROTOTYPE( 2010 );
	SENGEN_PROTOTYPE( 2011 );
	SENGEN_PROTOTYPE( 2012 );
	SENGEN_PROTOTYPE( 2013 );
	SENGEN_PROTOTYPE( 2014 );
	SENGEN_PROTOTYPE( 2015 );
	SENGEN_PROTOTYPE( 2016 );
	SENGEN_PROTOTYPE( 2017 );
	SENGEN_PROTOTYPE( 2018 );
	SENGEN_PROTOTYPE( 2019 );
	SENGEN_PROTOTYPE( 2020 );
	SENGEN_PROTOTYPE( 2021 );
	SENGEN_PROTOTYPE( 2022 );
	SENGEN_PROTOTYPE( 2023 );
	SENGEN_PROTOTYPE( 2024 );
	SENGEN_PROTOTYPE( 2025 );
	SENGEN_PROTOTYPE( 2026 );
	SENGEN_PROTOTYPE( 2027 );
	SENGEN_PROTOTYPE( 2028 );
	SENGEN_PROTOTYPE( 2029 );
	SENGEN_PROTOTYPE( 2030 );
	SENGEN_PROTOTYPE( 2031 );
	SENGEN_PROTOTYPE( 2032 );
	SENGEN_PROTOTYPE( 2033 );
	SENGEN_PROTOTYPE( 2034 );
	SENGEN_PROTOTYPE( 2035 );
	SENGEN_PROTOTYPE( 2036 );
	SENGEN_PROTOTYPE( 2037 );
	SENGEN_PROTOTYPE( 2038 );
	SENGEN_PROTOTYPE( 2039 );
	SENGEN_PROTOTYPE( 2040 );
	SENGEN_PROTOTYPE( 2041 );
	SENGEN_PROTOTYPE( 2042 );
	SENGEN_PROTOTYPE( 2043 );
	SENGEN_PROTOTYPE( 2044 );
	SENGEN_PROTOTYPE( 2045 );
	SENGEN_PROTOTYPE( 2046 );
	SENGEN_PROTOTYPE( 2047 );
	SENGEN_PROTOTYPE( 2048 );
	SENGEN_PROTOTYPE( 2049 );
	SENGEN_PROTOTYPE( 2050 );
	SENGEN_PROTOTYPE( 2051 );
	SENGEN_PROTOTYPE( 2052 );
	SENGEN_PROTOTYPE( 2053 );
	SENGEN_PROTOTYPE( 2054 );
	SENGEN_PROTOTYPE( 2055 );
	SENGEN_PROTOTYPE( 2056 );
	SENGEN_PROTOTYPE( 2057 );
	SENGEN_PROTOTYPE( 2058 );
	SENGEN_PROTOTYPE( 2059 );
	SENGEN_PROTOTYPE( 2060 );
	SENGEN_PROTOTYPE( 2061 );
	SENGEN_PROTOTYPE( 2062 );
	SENGEN_PROTOTYPE( 2063 );
	SENGEN_PROTOTYPE( 2064 );
	SENGEN_PROTOTYPE( 2065 );
	SENGEN_PROTOTYPE( 2066 );
	SENGEN_PROTOTYPE( 2067 );
	SENGEN_PROTOTYPE( 2068 );
	SENGEN_PROTOTYPE( 2069 );
	SENGEN_PROTOTYPE( 2070 );
	SENGEN_PROTOTYPE( 2071 );
	SENGEN_PROTOTYPE( 2072 );
	SENGEN_PROTOTYPE( 2073 );
	SENGEN_PROTOTYPE( 2074 );
	SENGEN_PROTOTYPE( 2075 );
	SENGEN_PROTOTYPE( 2076 );
	SENGEN_PROTOTYPE( 2077 );
	SENGEN_PROTOTYPE( 2078 );
	SENGEN_PROTOTYPE( 2079 );
	SENGEN_PROTOTYPE( 2080 );
	SENGEN_PROTOTYPE( 2081 );
	SENGEN_PROTOTYPE( 2082 );
	SENGEN_PROTOTYPE( 2083 );
	SENGEN_PROTOTYPE( 2084 );
	SENGEN_PROTOTYPE( 2085 );
	SENGEN_PROTOTYPE( 2086 );
	SENGEN_PROTOTYPE( 2087 );
	SENGEN_PROTOTYPE( 2088 );
	SENGEN_PROTOTYPE( 2089 );
	SENGEN_PROTOTYPE( 2090 );
	SENGEN_PROTOTYPE( 2091 );
	SENGEN_PROTOTYPE( 2092 );
	SENGEN_PROTOTYPE( 2093 );
	SENGEN_PROTOTYPE( 2094 );
	SENGEN_PROTOTYPE( 2095 );
	SENGEN_PROTOTYPE( 2096 );
	SENGEN_PROTOTYPE( 2097 );
	SENGEN_PROTOTYPE( 2098 );
	SENGEN_PROTOTYPE( 2099 );
	SENGEN_PROTOTYPE( 2100 );


	//2100番シリーズは圧縮式
	SENGEN_PROTOTYPE( 2101 );
	SENGEN_PROTOTYPE( 2102 );
	SENGEN_PROTOTYPE( 2103 );
	SENGEN_PROTOTYPE( 2104 );
	SENGEN_PROTOTYPE( 2105 );
	SENGEN_PROTOTYPE( 2106 );
	SENGEN_PROTOTYPE( 2107 );
	SENGEN_PROTOTYPE( 2108 );
	SENGEN_PROTOTYPE( 2109 );
	SENGEN_PROTOTYPE( 2110 );
	SENGEN_PROTOTYPE( 2111 );
	SENGEN_PROTOTYPE( 2112 );
	SENGEN_PROTOTYPE( 2113 );
	SENGEN_PROTOTYPE( 2114 );
	SENGEN_PROTOTYPE( 2115 );
	SENGEN_PROTOTYPE( 2116 );
	SENGEN_PROTOTYPE( 2117 );
	SENGEN_PROTOTYPE( 2118 );
	SENGEN_PROTOTYPE( 2119 );
	SENGEN_PROTOTYPE( 2120 );
	SENGEN_PROTOTYPE( 2121 );
	SENGEN_PROTOTYPE( 2122 );
	SENGEN_PROTOTYPE( 2123 );
	SENGEN_PROTOTYPE( 2124 );
	SENGEN_PROTOTYPE( 2125 );
	SENGEN_PROTOTYPE( 2126 );
	SENGEN_PROTOTYPE( 2127 );
	SENGEN_PROTOTYPE( 2128 );
	SENGEN_PROTOTYPE( 2129 );
	SENGEN_PROTOTYPE( 2130 );
	SENGEN_PROTOTYPE( 2131 );
	SENGEN_PROTOTYPE( 2132 );
	SENGEN_PROTOTYPE( 2133 );
	SENGEN_PROTOTYPE( 2134 );
	SENGEN_PROTOTYPE( 2135 );
	SENGEN_PROTOTYPE( 2136 );
	SENGEN_PROTOTYPE( 2137 );
	SENGEN_PROTOTYPE( 2138 );
	SENGEN_PROTOTYPE( 2139 );
	SENGEN_PROTOTYPE( 2140 );
	SENGEN_PROTOTYPE( 2141 );
	SENGEN_PROTOTYPE( 2142 );
	SENGEN_PROTOTYPE( 2143 );
	SENGEN_PROTOTYPE( 2144 );
	SENGEN_PROTOTYPE( 2145 );
	SENGEN_PROTOTYPE( 2146 );
	SENGEN_PROTOTYPE( 2147 );
	SENGEN_PROTOTYPE( 2148 );
	SENGEN_PROTOTYPE( 2149 );
	SENGEN_PROTOTYPE( 2150 );
	SENGEN_PROTOTYPE( 2151 );
	SENGEN_PROTOTYPE( 2152 );
	SENGEN_PROTOTYPE( 2153 );
	SENGEN_PROTOTYPE( 2154 );
	SENGEN_PROTOTYPE( 2155 );
	SENGEN_PROTOTYPE( 2156 );
	SENGEN_PROTOTYPE( 2157 );
	SENGEN_PROTOTYPE( 2158 );
	SENGEN_PROTOTYPE( 2159 );
	SENGEN_PROTOTYPE( 2160 );
	SENGEN_PROTOTYPE( 2161 );
	SENGEN_PROTOTYPE( 2162 );
	SENGEN_PROTOTYPE( 2163 );
	SENGEN_PROTOTYPE( 2164 );
	SENGEN_PROTOTYPE( 2165 );
	SENGEN_PROTOTYPE( 2166 );
	SENGEN_PROTOTYPE( 2167 );
	SENGEN_PROTOTYPE( 2168 );
	SENGEN_PROTOTYPE( 2169 );



	//3000番台はNEDO吸収シリーズ
	SENGEN_PROTOTYPE( 3000 );
	SENGEN_PROTOTYPE( 3001 );
	SENGEN_PROTOTYPE( 3002 );
	SENGEN_PROTOTYPE( 3003 );
	SENGEN_PROTOTYPE( 3004 );
	SENGEN_PROTOTYPE( 3005 );
	SENGEN_PROTOTYPE( 3006 );
	SENGEN_PROTOTYPE( 3007 );
	SENGEN_PROTOTYPE( 3008 );
	SENGEN_PROTOTYPE( 3009 );
	SENGEN_PROTOTYPE( 3010 );
	SENGEN_PROTOTYPE( 3011 );
	SENGEN_PROTOTYPE( 3012 );
	SENGEN_PROTOTYPE( 3013 );
	SENGEN_PROTOTYPE( 3014 );
	SENGEN_PROTOTYPE( 3015 );
	SENGEN_PROTOTYPE( 3016 );
	SENGEN_PROTOTYPE( 3017 );
	SENGEN_PROTOTYPE( 3018 );
	SENGEN_PROTOTYPE( 3019 );
	SENGEN_PROTOTYPE( 3020 );
	SENGEN_PROTOTYPE( 3021 );
	SENGEN_PROTOTYPE( 3022 );
	SENGEN_PROTOTYPE( 3023 );
	SENGEN_PROTOTYPE( 3024 );
	SENGEN_PROTOTYPE( 3025 );
	SENGEN_PROTOTYPE( 3026 );
	SENGEN_PROTOTYPE( 3027 );
	SENGEN_PROTOTYPE( 3028 );
	SENGEN_PROTOTYPE( 3029 );
	SENGEN_PROTOTYPE( 3030 );
	SENGEN_PROTOTYPE( 3031 );
	SENGEN_PROTOTYPE( 3032 );
	SENGEN_PROTOTYPE( 3033 );
	SENGEN_PROTOTYPE( 3034 );
	SENGEN_PROTOTYPE( 3035 );
	SENGEN_PROTOTYPE( 3036 );
	SENGEN_PROTOTYPE( 3037 );
	SENGEN_PROTOTYPE( 3038 );
	SENGEN_PROTOTYPE( 3039 );
	SENGEN_PROTOTYPE( 3040 );
	SENGEN_PROTOTYPE( 3041 );
	SENGEN_PROTOTYPE( 3042 );
	SENGEN_PROTOTYPE( 3043 );
	SENGEN_PROTOTYPE( 3044 );
	SENGEN_PROTOTYPE( 3045 );
	SENGEN_PROTOTYPE( 3046 );
	SENGEN_PROTOTYPE( 3047 );
	SENGEN_PROTOTYPE( 3048 );
	SENGEN_PROTOTYPE( 3049 );
	SENGEN_PROTOTYPE( 3050 );
	SENGEN_PROTOTYPE( 3051 );
	SENGEN_PROTOTYPE( 3052 );
	SENGEN_PROTOTYPE( 3053 );
	SENGEN_PROTOTYPE( 3054 );
	SENGEN_PROTOTYPE( 3055 );
	SENGEN_PROTOTYPE( 3056 );
	SENGEN_PROTOTYPE( 3057 );
	SENGEN_PROTOTYPE( 3058 );
	SENGEN_PROTOTYPE( 3059 );
	SENGEN_PROTOTYPE( 3060 );
	SENGEN_PROTOTYPE( 3061 );
	SENGEN_PROTOTYPE( 3062 );
	SENGEN_PROTOTYPE( 3063 );
	SENGEN_PROTOTYPE( 3064 );
	SENGEN_PROTOTYPE( 3065 );
	SENGEN_PROTOTYPE( 3066 );
	SENGEN_PROTOTYPE( 3067 );
	SENGEN_PROTOTYPE( 3068 );
	SENGEN_PROTOTYPE( 3069 );
	SENGEN_PROTOTYPE( 3070 );
	SENGEN_PROTOTYPE( 3071 );
	SENGEN_PROTOTYPE( 3072 );
	SENGEN_PROTOTYPE( 3073 );
	SENGEN_PROTOTYPE( 3074 );
	SENGEN_PROTOTYPE( 3075 );
	SENGEN_PROTOTYPE( 3076 );
	SENGEN_PROTOTYPE( 3077 );
	SENGEN_PROTOTYPE( 3078 );
	SENGEN_PROTOTYPE( 3079 );
	SENGEN_PROTOTYPE( 3080 );
	SENGEN_PROTOTYPE( 3081 );
	SENGEN_PROTOTYPE( 3082 );
	SENGEN_PROTOTYPE( 3083 );
	SENGEN_PROTOTYPE( 3084 );
	SENGEN_PROTOTYPE( 3085 );
	SENGEN_PROTOTYPE( 3086 );
	SENGEN_PROTOTYPE( 3087 );
	SENGEN_PROTOTYPE( 3088 );
	SENGEN_PROTOTYPE( 3089 );
	SENGEN_PROTOTYPE( 3090 );
	SENGEN_PROTOTYPE( 3091 );
	SENGEN_PROTOTYPE( 3092 );
	SENGEN_PROTOTYPE( 3093 );
	SENGEN_PROTOTYPE( 3094 );
	SENGEN_PROTOTYPE( 3095 );
	SENGEN_PROTOTYPE( 3096 );
	SENGEN_PROTOTYPE( 3097 );
	SENGEN_PROTOTYPE( 3098 );
	SENGEN_PROTOTYPE( 3099 );
	SENGEN_PROTOTYPE( 3100 );


	//中尾研究室用4000番台
	SENGEN_PROTOTYPE( 4001 );
	SENGEN_PROTOTYPE( 4002 );
	SENGEN_PROTOTYPE( 4003 );
	SENGEN_PROTOTYPE( 4004 );
	SENGEN_PROTOTYPE( 4005 );
	SENGEN_PROTOTYPE( 4006 );
	SENGEN_PROTOTYPE( 4007 );
	SENGEN_PROTOTYPE( 4008 );
	SENGEN_PROTOTYPE( 4009 );
	SENGEN_PROTOTYPE( 4010 );
	SENGEN_PROTOTYPE( 4011 );
	SENGEN_PROTOTYPE( 4012 );
	SENGEN_PROTOTYPE( 4013 );
	SENGEN_PROTOTYPE( 4014 );
	SENGEN_PROTOTYPE( 4015 );
	SENGEN_PROTOTYPE( 4016 );
	SENGEN_PROTOTYPE( 4017 );
	SENGEN_PROTOTYPE( 4018 );
	SENGEN_PROTOTYPE( 4019 );
	SENGEN_PROTOTYPE( 4020 );
	SENGEN_PROTOTYPE( 4021 );
	SENGEN_PROTOTYPE( 4022 );
	SENGEN_PROTOTYPE( 4023 );
	SENGEN_PROTOTYPE( 4024 );
	SENGEN_PROTOTYPE( 4025 );
	SENGEN_PROTOTYPE( 4026 );
	SENGEN_PROTOTYPE( 4027 );
	SENGEN_PROTOTYPE( 4028 );
	SENGEN_PROTOTYPE( 4029 );
	SENGEN_PROTOTYPE( 4030 );
	SENGEN_PROTOTYPE( 4031 );
	SENGEN_PROTOTYPE( 4032 );
	SENGEN_PROTOTYPE( 4033 );
	SENGEN_PROTOTYPE( 4034 );
	SENGEN_PROTOTYPE( 4035 );
	SENGEN_PROTOTYPE( 4036 );
	SENGEN_PROTOTYPE( 4037 );
	SENGEN_PROTOTYPE( 4038 );
	SENGEN_PROTOTYPE( 4039 );
	SENGEN_PROTOTYPE( 4040 );
	SENGEN_PROTOTYPE( 4041 );
	SENGEN_PROTOTYPE( 4042 );
	SENGEN_PROTOTYPE( 4043 );
	SENGEN_PROTOTYPE( 4044 );
	SENGEN_PROTOTYPE( 4045 );
	SENGEN_PROTOTYPE( 4046 );
	SENGEN_PROTOTYPE( 4047 );
	SENGEN_PROTOTYPE( 4048 );
	SENGEN_PROTOTYPE( 4049 );
	SENGEN_PROTOTYPE( 4050 );
	SENGEN_PROTOTYPE( 4051 );
	SENGEN_PROTOTYPE( 4052 );
	SENGEN_PROTOTYPE( 4053 );
	SENGEN_PROTOTYPE( 4054 );
	SENGEN_PROTOTYPE( 4055 );
	SENGEN_PROTOTYPE( 4056 );
	SENGEN_PROTOTYPE( 4057 );
	SENGEN_PROTOTYPE( 4058 );
	SENGEN_PROTOTYPE( 4059 );
	SENGEN_PROTOTYPE( 4060 );
	SENGEN_PROTOTYPE( 4061 );
	SENGEN_PROTOTYPE( 4062 );
	SENGEN_PROTOTYPE( 4063 );
	SENGEN_PROTOTYPE( 4064 );
	SENGEN_PROTOTYPE( 4065 );
	SENGEN_PROTOTYPE( 4066 );
	SENGEN_PROTOTYPE( 4067 );
	SENGEN_PROTOTYPE( 4068 );
	SENGEN_PROTOTYPE( 4069 );
	SENGEN_PROTOTYPE( 4070 );
	SENGEN_PROTOTYPE( 4071 );
	SENGEN_PROTOTYPE( 4072 );
	SENGEN_PROTOTYPE( 4073 );
	SENGEN_PROTOTYPE( 4074 );
	SENGEN_PROTOTYPE( 4075 );
	SENGEN_PROTOTYPE( 4076 );
	SENGEN_PROTOTYPE( 4077 );
	SENGEN_PROTOTYPE( 4078 );
	SENGEN_PROTOTYPE( 4079 );
	SENGEN_PROTOTYPE( 4080 );
	SENGEN_PROTOTYPE( 4081 );
	SENGEN_PROTOTYPE( 4082 );
	SENGEN_PROTOTYPE( 4083 );
	SENGEN_PROTOTYPE( 4084 );
	SENGEN_PROTOTYPE( 4085 );
	SENGEN_PROTOTYPE( 4086 );
	SENGEN_PROTOTYPE( 4087 );
	SENGEN_PROTOTYPE( 4088 );
	SENGEN_PROTOTYPE( 4089 );
	SENGEN_PROTOTYPE( 4090 );
	SENGEN_PROTOTYPE( 4091 );
	SENGEN_PROTOTYPE( 4092 );
	SENGEN_PROTOTYPE( 4093 );
	SENGEN_PROTOTYPE( 4094 );
	SENGEN_PROTOTYPE( 4095 );
	SENGEN_PROTOTYPE( 4096 );
	SENGEN_PROTOTYPE( 4097 );
	SENGEN_PROTOTYPE( 4098 );
	SENGEN_PROTOTYPE( 4099 );
	SENGEN_PROTOTYPE( 4100 );
	
	SENGEN_PROTOTYPE( 5010 );
	SENGEN_PROTOTYPE( 5011 );


	//5000番ってなんだっけ？サンプルだっけ？
	void CM5001init();void CM5001setup();void CM5001calc(int);void CM5001outletCheck(int);void CM5001renew();
	void CM5002init();void CM5002setup();void CM5002calc(int);void CM5002outletCheck(int);void CM5002renew();
	void CM5003init();void CM5003setup();void CM5003calc(int);void CM5003outletCheck(int);void CM5003renew();
	void CM5004init();void CM5004setup();void CM5004calc(int);void CM5004outletCheck(int);void CM5004renew();
	void CM5005init();void CM5005setup();void CM5005calc(int);void CM5005outletCheck(int);void CM5005renew();
	void CM5006init();void CM5006setup();void CM5006calc(int);void CM5006outletCheck(int);void CM5006renew();
	void CM5007init();void CM5007setup();void CM5007calc(int);void CM5007outletCheck(int);void CM5007renew();
	void CM5008init();void CM5008setup();void CM5008calc(int);void CM5008outletCheck(int);void CM5008renew();
	void CM5009init();void CM5009setup();void CM5009calc(int);void CM5009outletCheck(int);void CM5009renew();
//	void CM5010init();void CM5010setup();void CM5010calc(int);void CM5010outletCheck(int);void CM5010renew();
//	void CM5011init();void CM5011setup();void CM5011calc(int);void CM5011outletCheck(int);void CM5011renew();
	void CM5012init();void CM5012setup();void CM5012calc(int);void CM5012outletCheck(int);void CM5012renew();
	void CM5013init();void CM5013setup();void CM5013calc(int);void CM5013outletCheck(int);void CM5013renew();
	void CM5014init();void CM5014setup();void CM5014calc(int);void CM5014outletCheck(int);void CM5014renew();
	void CM5015init();void CM5015setup();void CM5015calc(int);void CM5015outletCheck(int);void CM5015renew();
	void CM5016init();void CM5016setup();void CM5016calc(int);void CM5016outletCheck(int);void CM5016renew();
	void CM5017init();void CM5017setup();void CM5017calc(int);void CM5017outletCheck(int);void CM5017renew();
	void CM5018init();void CM5018setup();void CM5018calc(int);void CM5018outletCheck(int);void CM5018renew();
	void CM5019init();void CM5019setup();void CM5019calc(int);void CM5019outletCheck(int);void CM5019renew();
	void CM5020init();void CM5020setup();void CM5020calc(int);void CM5020outletCheck(int);void CM5020renew();
	void CM5021init();void CM5021setup();void CM5021calc(int);void CM5021outletCheck(int);void CM5021renew();
	void CM5022init();void CM5022setup();void CM5022calc(int);void CM5022outletCheck(int);void CM5022renew();
	void CM5023init();void CM5023setup();void CM5023calc(int);void CM5023outletCheck(int);void CM5023renew();
	void CM5024init();void CM5024setup();void CM5024calc(int);void CM5024outletCheck(int);void CM5024renew();
	void CM5025init();void CM5025setup();void CM5025calc(int);void CM5025outletCheck(int);void CM5025renew();
	void CM5026init();void CM5026setup();void CM5026calc(int);void CM5026outletCheck(int);void CM5026renew();
	void CM5027init();void CM5027setup();void CM5027calc(int);void CM5027outletCheck(int);void CM5027renew();
	void CM5028init();void CM5028setup();void CM5028calc(int);void CM5028outletCheck(int);void CM5028renew();
	void CM5029init();void CM5029setup();void CM5029calc(int);void CM5029outletCheck(int);void CM5029renew();
	void CM5030init();void CM5030setup();void CM5030calc(int);void CM5030outletCheck(int);void CM5030renew();
	void CM5031init();void CM5031setup();void CM5031calc(int);void CM5031outletCheck(int);void CM5031renew();
	void CM5032init();void CM5032setup();void CM5032calc(int);void CM5032outletCheck(int);void CM5032renew();
	void CM5033init();void CM5033setup();void CM5033calc(int);void CM5033outletCheck(int);void CM5033renew();
	void CM5034init();void CM5034setup();void CM5034calc(int);void CM5034outletCheck(int);void CM5034renew();
	void CM5035init();void CM5035setup();void CM5035calc(int);void CM5035outletCheck(int);void CM5035renew();
	void CM5036init();void CM5036setup();void CM5036calc(int);void CM5036outletCheck(int);void CM5036renew();
	void CM5037init();void CM5037setup();void CM5037calc(int);void CM5037outletCheck(int);void CM5037renew();
	void CM5038init();void CM5038setup();void CM5038calc(int);void CM5038outletCheck(int);void CM5038renew();
	void CM5039init();void CM5039setup();void CM5039calc(int);void CM5039outletCheck(int);void CM5039renew();
	void CM5040init();void CM5040setup();void CM5040calc(int);void CM5040outletCheck(int);void CM5040renew();
	void CM5041init();void CM5041setup();void CM5041calc(int);void CM5041outletCheck(int);void CM5041renew();
	void CM5042init();void CM5042setup();void CM5042calc(int);void CM5042outletCheck(int);void CM5042renew();
	void CM5043init();void CM5043setup();void CM5043calc(int);void CM5043outletCheck(int);void CM5043renew();
	void CM5044init();void CM5044setup();void CM5044calc(int);void CM5044outletCheck(int);void CM5044renew();
	void CM5045init();void CM5045setup();void CM5045calc(int);void CM5045outletCheck(int);void CM5045renew();
	void CM5046init();void CM5046setup();void CM5046calc(int);void CM5046outletCheck(int);void CM5046renew();
	void CM5047init();void CM5047setup();void CM5047calc(int);void CM5047outletCheck(int);void CM5047renew();
	void CM5048init();void CM5048setup();void CM5048calc(int);void CM5048outletCheck(int);void CM5048renew();
	void CM5049init();void CM5049setup();void CM5049calc(int);void CM5049outletCheck(int);void CM5049renew();
	void CM5050init();void CM5050setup();void CM5050calc(int);void CM5050outletCheck(int);void CM5050renew();



	//9000番以上は制御モジュール

	SENGEN_PROTOTYPE( 9001 );
	SENGEN_PROTOTYPE( 9002 );
	SENGEN_PROTOTYPE( 9003 );
	SENGEN_PROTOTYPE( 9004 );
	SENGEN_PROTOTYPE( 9005 );
	SENGEN_PROTOTYPE( 9006 );
	SENGEN_PROTOTYPE( 9007 );
	SENGEN_PROTOTYPE( 9008 );
	SENGEN_PROTOTYPE( 9009 );
	SENGEN_PROTOTYPE( 9010 );
	SENGEN_PROTOTYPE( 9011 );
	SENGEN_PROTOTYPE( 9012 );
	SENGEN_PROTOTYPE( 9013 );
	SENGEN_PROTOTYPE( 9014 );
	SENGEN_PROTOTYPE( 9015 );
	SENGEN_PROTOTYPE( 9016 );
	SENGEN_PROTOTYPE( 9017 );
	SENGEN_PROTOTYPE( 9018 );
	SENGEN_PROTOTYPE( 9019 );
	SENGEN_PROTOTYPE( 9020 );
	SENGEN_PROTOTYPE( 9021 );
	SENGEN_PROTOTYPE( 9022 );
	SENGEN_PROTOTYPE( 9023 );
	SENGEN_PROTOTYPE( 9024 );
	SENGEN_PROTOTYPE( 9025 );
	SENGEN_PROTOTYPE( 9026 );
	SENGEN_PROTOTYPE( 9027 );
	SENGEN_PROTOTYPE( 9028 );
	SENGEN_PROTOTYPE( 9029 );
	SENGEN_PROTOTYPE( 9030 );
	SENGEN_PROTOTYPE( 9031 );
	SENGEN_PROTOTYPE( 9032 );
	SENGEN_PROTOTYPE( 9033 );
	SENGEN_PROTOTYPE( 9034 );
	SENGEN_PROTOTYPE( 9035 );
	SENGEN_PROTOTYPE( 9036 );
	SENGEN_PROTOTYPE( 9037 );
	SENGEN_PROTOTYPE( 9038 );
	SENGEN_PROTOTYPE( 9039 );
	SENGEN_PROTOTYPE( 9040 );
	SENGEN_PROTOTYPE( 9041 );
	SENGEN_PROTOTYPE( 9042 );
	SENGEN_PROTOTYPE( 9043 );
	SENGEN_PROTOTYPE( 9044 );
	SENGEN_PROTOTYPE( 9045 );
	SENGEN_PROTOTYPE( 9046 );
	SENGEN_PROTOTYPE( 9047 );
	SENGEN_PROTOTYPE( 9048 );
	SENGEN_PROTOTYPE( 9049 );
	SENGEN_PROTOTYPE( 9050 );
	SENGEN_PROTOTYPE( 9051 );
	SENGEN_PROTOTYPE( 9052 );
	SENGEN_PROTOTYPE( 9053 );
	SENGEN_PROTOTYPE( 9054 );
	SENGEN_PROTOTYPE( 9055 );
	SENGEN_PROTOTYPE( 9056 );
	SENGEN_PROTOTYPE( 9057 );
	SENGEN_PROTOTYPE( 9058 );
	SENGEN_PROTOTYPE( 9059 );
	SENGEN_PROTOTYPE( 9060 );
	SENGEN_PROTOTYPE( 9061 );
	SENGEN_PROTOTYPE( 9062 );
	SENGEN_PROTOTYPE( 9063 );
	SENGEN_PROTOTYPE( 9064 );
	SENGEN_PROTOTYPE( 9065 );
	SENGEN_PROTOTYPE( 9066 );
	SENGEN_PROTOTYPE( 9067 );
	SENGEN_PROTOTYPE( 9068 );
	SENGEN_PROTOTYPE( 9069 );
	SENGEN_PROTOTYPE( 9070 );
	SENGEN_PROTOTYPE( 9071 );
	SENGEN_PROTOTYPE( 9072 );
	SENGEN_PROTOTYPE( 9073 );
	SENGEN_PROTOTYPE( 9074 );
	SENGEN_PROTOTYPE( 9075 );
	SENGEN_PROTOTYPE( 9076 );
	SENGEN_PROTOTYPE( 9077 );
	SENGEN_PROTOTYPE( 9078 );
	SENGEN_PROTOTYPE( 9079 );
	SENGEN_PROTOTYPE( 9080 );
	SENGEN_PROTOTYPE( 9081 );
	SENGEN_PROTOTYPE( 9082 );
	SENGEN_PROTOTYPE( 9083 );
	SENGEN_PROTOTYPE( 9084 );
	SENGEN_PROTOTYPE( 9085 );
	SENGEN_PROTOTYPE( 9086 );
	SENGEN_PROTOTYPE( 9087 );
	SENGEN_PROTOTYPE( 9088 );
	SENGEN_PROTOTYPE( 9089 );
	SENGEN_PROTOTYPE( 9090 );
	SENGEN_PROTOTYPE( 9091 );
	SENGEN_PROTOTYPE( 9092 );
	SENGEN_PROTOTYPE( 9093 );
	SENGEN_PROTOTYPE( 9094 );
	SENGEN_PROTOTYPE( 9095 );
	SENGEN_PROTOTYPE( 9096 );
	SENGEN_PROTOTYPE( 9097 );
	SENGEN_PROTOTYPE( 9098 );
	SENGEN_PROTOTYPE( 9099 );
	SENGEN_PROTOTYPE( 9100 );

};



#endif
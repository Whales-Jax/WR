/********************************************************************************
 * DXProperty_ver04.h                                                           *
 *                                                                              *
 * 2008/04 開発チーム                                                           *
 * Propertyに２成分用のセットアップ関数追加．by小林                             *
 * Propertyに複成分用のセットアップ関数追加．ver.3.0 by平松						*
 ********************************************************************************/

#ifndef __DXPROPERTY_VER04_H__
#define __DXPROPERTY_VER04_H__

#include <windows.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include "CCSVFile.h"
using namespace std;

/*********************************************************************************
 * Fluid クラス -- 任意の状態の物性を格納するクラス
 *
 * メンバ変数
 *     T       -- 温度 deg.C
 *     P       -- 圧力 kPa
 *     rho     -- 密度 kg/m^3
 *     h       -- 比エンタルピ kJ/kg
 *     s       -- 比エントロピ kJ/(kg*K)
 *     x       -- 乾き度 -
 *     cp      -- 定圧比熱 kJ/(kg*K)
 *     thc     -- 熱伝導率 W/(m*K)
 *     visc    -- 粘性 uPa*s
 *     st      -- 表面張力 N/m
 * メンバ変数
 *     Fluid   -- コンストラクタ
 *********************************************************************************/
class DXFluid
{
	public:
		double T;
		double P;
		double rho;
		double h;
		double s;
		double x;
		double cp;
		double cv;
		double thc;
		double visc;
		double st;
		double u;

		// コンストラクタ
		DXFluid()
		{
			T    = 0.0;
			P    = 0.0;
			rho  = 0.0;
			h    = 0.0;
			s    = 0.0;
			x    = 0.0;
			cp   = 0.0;
			cv   = 0.0;
			thc  = 0.0;
			visc = 0.0;
			st   = 0.0;
		}

		~DXFluid(){}
};
/*********************************************************************************
 * Property クラス -- refpropを利用して冷媒の物性値を計算するクラス
 *
 * メンバ関数
 *     DllSetup         -- コンストラクタ
 *     ~DllSetup        -- デストラクタ
 *     freeDll          -- DLL解放関数
 *     setup            -- セットアップ関数
 *     setref           -- 基準状態設定関数
 *     xmole            -- モル分率・分子量関数
 *     xmass            -- 質量分率・分子量関数
 *     criticalPoint    -- 臨界状態関数
 *     state_tp         -- 単相状態関数(温度・圧力入力)
 *     state_ph         -- 任意状態関数(圧力・比エンタルピ入力)
 *     sat_t            -- 飽和状態関数(温度入力)
 *     sat_p            -- 飽和状態関数(圧力入力)
 *********************************************************************************/
// DLL呼び出し関数の型宣言
typedef void (__stdcall *SetupDllType) (long&,char*,char*,char*,long&,char*,long,long,long,long);
typedef void (__stdcall *SetrefDllType) (char*,long&,double*,double&,double&,double&,double&,long&,char*,long,long);
typedef void (__stdcall *XmassDllType) (double*,double*,double&);
typedef void (__stdcall *XmoleDllType) (double*,double*,double&);
typedef void (__stdcall *EntDllType) (double&,double&,double*,double&);
typedef void (__stdcall *CritpDllType) (double*,double&,double&,double&,long&,char*,long);
typedef void (__stdcall *SatDllType) (double&,double*,long&,double&,double&,double&,double*,double*,long&,char*,long);
typedef void (__stdcall *PhflshDllType) (double&,double&,double*,double&,double&,double&,double&,double*,double*,double&,double&,double&,double&,double&,double&,long&,char*,long);
typedef void (__stdcall *PtflshDllType) (double&,double&,double*,double&,double&,double&,double*,double*,double&,double&,double&,double&,double&,double&,double&,long&,char*,long);
typedef void (__stdcall *SurtenDllType) (double&,double&,double&,double*,double*,double&,long&,char*,long);
typedef void (__stdcall *SurftDllType)(double &,double &,double *,double &,long &,char*,long);
typedef void (__stdcall *ThermdllTYPE)  (double&,double&,double*,double&,double&,double&,double&,double&,double&,double&,double&);
typedef void (__stdcall *TpflshDllType) (double&,double&,double*,double&,double&,double&,double*,double*,double&, double&,double&,double&,double&,double&,double&,long&,char*,long);
typedef void (__stdcall *TprhoDllType)(double&,double&,double*,long&,long&,double&,long&,char*,long);
typedef void (__stdcall *TrnprpDllType)(double&,double&,double*,double&,double&,long&,char*,long);

class DXProperty
{
	private:
		// DLLのハンドル
		HINSTANCE hRefpropDll;
		// DLL内の関数のエントリ
		SetupDllType  setupDll;
		SetrefDllType setrefDll;
		XmassDllType  xmassDll;
		XmoleDllType  xmoleDll;

		CritpDllType  critpDll;
		SatDllType    sattDll;
		SatDllType    satpDll;
		PhflshDllType phflshDll;
		PtflshDllType ptflshDll;
		SurtenDllType surtenDll;
		SurftDllType  surftDll;
		ThermdllTYPE  thermDll;
		TpflshDllType tpflshDll;
		TprhoDllType  tprhoDll;
		TrnprpDllType trnprpDll;

		// 成分数
		long numberOfComponents_;
		// モル分率
		double *moleFraction;

	public:
		// 質量分率・分子量
		double *massFraction;
		double molecularWeight;
		ofstream refdata;

		DXFluid **data_table;	//!<state_ph用データテーブル
		DXFluid **data_table_p;	//!<sat_p用データテーブル
		DXFluid **data_table_tp;//!<state_tp用データテーブル
		DXFluid **data_table_t;	//!<sat_t用データテーブル
		DXFluid **data_table_du;//!<state_vu用データテーブル
		int rep_h;
		int rep_p;
		int rep_t;
		int rep_d;
		int rep_u;
		bool _DLLsetup;

		int calccount;

		double min_P;
		double min_h;
		double min_T;
		double max_P;
		double max_h;
		double max_T;

		string PATH;
		ostringstream name;

		string reibai;
		string joutai;
		DXFluid Rc,Rl,Rv;


		// コンストラクタ
		explicit DXProperty();
		// デストラクタ
		virtual ~DXProperty();
		// DLL解放関数
		void freeDll();
		void LoadDLL(std::string DllFileName);
		// セットアップ関数
		void setup();
		// モル分率・分子量関数
		void xmole(double *massFraction, double *moleFraction, double &molecularWeight);
		// 質量分率・分子量関数
		void xmass(double *moleFraction, double *massFraction, double &molecularWeight);
		// 分子量関数
		double molecularWeight_f();
		// 臨界状態関数
		void criticalPoint(DXFluid &crit);
		// 表面張力関数
		double surfaceTension_t(double temperature);

		// 単相状態関数
		void state_tp(double temperature, double pressure);
		// 任意状態関数
		void state_ph(double pressure, double enthalpy);
		// 飽和状態関数(温度入力)
		void sat_t(double temperature);
		// 飽和状態関数(圧力入力)
		void sat_p(double pressure);




		void CopyTable( DXProperty &ref );
		void set_step( int h , int p , int t);
		void init_table(void);
		void load_table(std::string fname);
		void load_table2(std::string fname);
		void make_table(double start_h , double end_h , double start_p , double end_p ,double start_t , double end_t );
		bool state_ph2(double pressure, double enthalpy);
		bool sat_p2(double pressure);
		bool state_tp2(double temperature , double pressure );
		bool sat_t2(double temperature);

		bool state_ph3(double pressure, double enthalpy);
		bool sat_p3(double pressure);
		bool state_tp3(double temperature , double pressure );
		bool sat_t3(double temperature);


		//state_vu
		void init_table_du(void);
		void load_table_du(std::string fname);
		void make_table_du(double start_d , double end_d , double start_u , double end_u );
		bool state_du2(double density , double energy );

};

#endif
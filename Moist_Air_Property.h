/***********************************************************
 * Air_Property.h
 *
 * 2008/07 平松雅弘
 ***********************************************************/
/*
	空気物性の計算は『M.CONDE ENGINEERING, 2007』を参照すること
	乾き空気と水蒸気の熱伝導率と粘性係数の式は出所がわからず精度が不明なので，
	別のものを使用．
	乾き空気は
		Journal of Physical and chemical reference data, 1985, 14,(4), P947-970
		『Viscosity and Thermal Conductivity of Dry Air in the Gaseous Phase』
		K. Kadoya, N. Matsunaga, A. Nagashima
	を参照
	水蒸気は
		『property_water_ver2.h』	（1980 日本機械学会蒸気表）
	を参照すること
	水の物性に関しては新しい国際基準になっていないので，徹底してやりたい人は最新
	の日本機械学会蒸気表を見れば式がわかります．また，水の物性は0.0度以下では使
	えないので，計算範囲が同様のrefpropを使うことも検討中です．でも，熱伝導率と
	粘性係数を計算するためだけにファイルが増えるのも微妙なので，今のところは実装
	する気はないです．

	注）property_water_ver2.hを使用しているのでincludeするのを忘れずに！！！
*/
#include <cmath>

#ifndef __MOIST_AIR_PROPERTY_H__
#define __MOIST_AIR_PROPERTY_H__

/*********************************************************************************
 * Air クラス -- 任意の状態の物性を格納するクラス
 *
 * メンバ変数
 *     T       -- 温度 deg.C
 *     P       -- 圧力 kPa
 *     rho     -- 密度 kg/m^3
 *     h       -- 比エンタルピ kJ/kg(DA)
 *     s       -- 比エントロピ kJ/(kg(DA)*K)
 *     x       -- 絶対湿度 kg/kg(DA)
 *	   phi	   -- 相対湿度 -
 *     cp      -- 定圧比熱 kJ/(kg(DA)*K)
 *     thc     -- 熱伝導率 kW/(m*K)
 *     visc    -- 粘性 Pa*s
 *     Psat    -- 飽和水蒸気圧力 kPa
 *     G       -- 流量 kg(DA)/s
 *     dP      -- 圧損 kPa
 * メンバ関数
 *     Air     -- コンストラクタ
 *     Air_Ptx -- 状態関数（圧力・温度・絶対湿度入力）
 *     Air_Phx -- 状態関数（圧力・比エンタルピー・絶対湿度入力）
 ********************************************************************************/
class Air
{
	public:
		double T;
		double P;
		double rho;
		double h;
		double x;
		double cp;
		double thc;
		double visc;
		double Psat;
		double G;
		double dP;

		double _P;
		double _h;
		double _x;

		bool okka;
		int LoopCount;

		// コンストラクタ
		Air()
		{
			T	 = 0.0;
			P	 = 101.325; // 大気圧
			rho	 = 0.0;
			h	 = 0.0;
			x	 = 0.0;
			cp	 = 0.0;
			thc	 = 0.0;
			visc = 0.0;
			Psat = 0.0;
			G    = 0.0;
			dP   = 0.0;
		}

		// 状態関数（圧力・温度・絶対湿度入力）
		void Air_Ptx ( double pressure, double temperature, double humidity );
		
		// 状態関数（圧力・比エンタルピー・絶対湿度入力）
		void Air_Phx ( double pressure, double enthalpy, double humidity );

		// 下記簡易計算用
		// 状態関数（温度・絶対湿度入力）
		void Air_tx ( double temperature, double humidity );

		// 状態関数（比エンタルピー・絶対湿度入力）
		void Air_hx ( double enthalpy, double humidity );
};

// 以下，湿り空気物性関数

// 湿り空気比容積関数	-100 <= T <= 200
double Air_Ptxv ( double P, double T, double x );
// 湿り空気比容積関数	-60 <= T <= 200
double Air_Ptphiv ( double P, double T, double phi );

// 湿り空気密度関数	-100 <= T <= 200
double Air_Ptxrho ( double P, double T, double x );
// 湿り空気密度関数	-60 <= T <= 200
double Air_Ptphirho ( double P, double T, double phi );

// 湿り空気比エンタルピー関数	-100 <= T <= 200
double Air_Ptxh ( double P, double T, double x );
// 湿り空気比エンタルピー関数	-60 <= T <= 200
double Air_Ptphih ( double P, double T, double phi );

// 湿り空気比エントロピー関数	-100 <= T <= 200
double Air_Ptxs ( double P, double T, double x );
// 湿り空気比エントロピー関数	-60 <= T <= 200
double Air_Ptphis ( double P, double T, double phi );

// 湿り空気相対湿度関数	-60 <= T <= 365
double Air_Ptxphi ( double P, double T, double x );
// 湿り空気相対湿度関数	-100 <= T <= 200
double Air_Pttwbphi ( double P, double T, double Twb );

// 湿り空気絶対湿度関数	-60 <= T <= 365
double Air_Ptphix ( double P, double T, double phi );
// 湿り空気絶対湿度関数	-100 <= T <= 200
double Air_Phtx ( double P, double h, double T );
// 湿り空気絶対湿度関数	-60 <= T <= 200
double Air_Pttwbx ( double P, double T, double Twb );

// 湿り空気温度関数	-100 <= T <= 200
double Air_Phxt ( double P, double h, double x );
// 湿り空気温度関数	-60 <= T <= 200
double Air_Phphit ( double P, double h, double phi );

// 飽和水蒸気圧関数	-60 <= T <= 365
double Air_tPs ( double T );

// 湿り空気定圧比熱関数	-100 <= T <= 200
double Air_Ptxcp ( double P, double T, double x );

// 湿り空気熱伝導率関数	0.01 <= T <= 350
double Air_Ptxthc ( double P, double T, double x );

// 湿り空気粘性係数関数	0.01 <= T <= 350
double Air_Ptxvisc ( double P, double T, double x );

// 湿り空気動粘性係数関数	0.01 <= T <= 200
double Air_Ptxkine ( double P, double T, double x );

// 湿り空気熱拡散率関数	0.01 <= T <= 200
double Air_Ptxthd ( double P, double T, double x );

// 湿り空気水蒸気の拡散係数関数	-20 <= T <= 300
double Air_PtDh2o ( double P, double T );

// プラントル数関数	0.01 <= T <= 200
double Air_PtxPr ( double P, double T, double x );

// シュミット数関数	0.01 <= T <= 200
double Air_PtxSc ( double P, double T, double x );

// 露点温度関数	x >= 0.0
double Air_Pxtdp ( double P, double x );

// 凝固点温度関数	P <= 200.0[MPa]
double Air_Ptmel ( double P );



// 下記の関数は上記の関数内で使用する関数

// 湿球温度関数	-60 <= T <= 200
double Air_Ptxtwb ( double P, double T, double x );
// 湿球温度関数	-100 <= T <= 200
double Air_Ptphitwb ( double P, double T, double phi );

// 飽和凝縮相比エンタルピー関数	-100 <= T <= 200
double Air_Pthc ( double P, double T );

// Bm, Cm計算用関数
void Bm_Cm ( double T, double Xa, double Xw, double &Bm, double &Cm );

// 修正係数f算出関数
double f_ ( double P, double T );

// 乾き空気熱伝導率関数	-193 <= T <= 1727
double Dair_Ptthc ( double P, double T );

// 乾き空気粘性係数関数	-193 <= T <= 1727
double Dair_Ptvisc ( double P, double T );



// 以下，大気圧（101.325kPa）での簡易空気関数

//比体積関数
inline double Air_txv ( double T, double x )
{
	return 4.554 * 1e-3 * ( 273.15 + T ) * ( 0.622 + x );
}

// 密度関数
inline double Air_txrho ( double T, double x )
{
	return 1.0 / Air_txv ( T, x );
}

// 比エンタルピー関数
inline double Air_txh ( double T, double x )
{
	return ( 1005.22 + 0.02615 * T ) * T + x * ( 2500800 + 1868 * T );
}

// 空気温度関数
inline double Air_hxt ( double h, double x )
{
	double a,b,c,d;
	double T;

	a = 1005.22;
	b = 0.02615;
	c = 2500800.0;
	d = 1868.0;

	// 温度
	T = ( -( d * x + a ) + pow( pow( x * d + a, 2.0 ) - 4.0 * b * ( x * c - h ), 0.5 ) ) / ( 2.0 * b );

	return ( T );
}

// 熱伝導率
inline double Air_tthc ( double T )
{
	return 0.02624 + 7.58e-5 * ( T - 27.0 );
}

// 粘性係数
inline double Air_tvisc ( double T )
{
	return 1.859e-5 + 4.32e-8 * ( T - 27.0 );
}


// デシカント計算用
// 水蒸気の蒸発潜熱
inline double Wat_tlh(double t)//latentheat
{
	return 2500800 + 1868 * t;//2358500.0 - 2460.0 * ( t - 60.0 );
}

#endif
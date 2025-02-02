/***********************************************************
 * Air_Property.cpp
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

// 07/08 一部変更 by 平松

#include <iostream>
#include <cmath>
#include "Moist_Air_Property.h"
#include "property_water_ver2.h"

const double Error = 1e-9; // 収束判定
const double Delta = 1.0 + 1e-3; // 刻み幅
const int Loop = 1000; // 繰り返し計算回数
int LoopCount;

const double Ma = 28.9645; // 空気の分子量 g/mol
const double Mw = 18.01528; // 水の分子量 g/mol
const double R = 8.31441; // 気体定数 J/(mol/K)

/********************************************************************************
 * Air::Air_Ptx -- 任意状態関数(圧力・温度・絶対湿度入力)
 *                 任意の圧力・温度・絶対湿度における状態の算出
 *
 * 引数(入力値)
 *     pressure     -- 圧力 kPa
 *     temperature  -- 温度 deg.C
 *     humidity     -- 絶対湿度 kg/kg(DA)
 ********************************************************************************/
void Air::Air_Ptx ( double pressure, double temperature, double humidity )
{
	this->P = pressure;
	this->T = temperature;
	this->x = humidity;
	this->rho = Air_Ptxrho ( pressure, temperature, humidity );
	this->h = Air_Ptxh ( pressure, temperature, humidity );
	//this->cp = Air_Ptxcp ( pressure, temperature, humidity );
	this->thc = Air_Ptxthc ( pressure, temperature, humidity );
	this->visc = Air_Ptxvisc ( pressure, temperature, humidity );
	this->Psat = Air_tPs ( temperature );
}

/********************************************************************************
 * Air::Air_Phx -- 任意状態関数(圧力・比エントロピー・絶対湿度入力)
 *                 任意の圧力・比エントロピー・絶対湿度における状態の算出
 *
 * 引数(入力値)
 *     pressure    -- 圧力 kPa
 *     enthalpy    -- 比エンタルピー kJ/kg(DA)
 *     humidity    -- 絶対湿度 kg/kg(DA)
 ********************************************************************************/
void Air::Air_Phx ( double pressure, double enthalpy, double humidity )
{
	if( _P == pressure &&
		_h == enthalpy &&
		_x == humidity){
			return;
	}

	okka = true;

	double temperature;

	this->P = pressure;
	this->T = temperature = Air_Phxt ( pressure, enthalpy, humidity );
	this->x = humidity;
	this->rho = Air_Ptxrho ( pressure, temperature, humidity );
	this->h = enthalpy;
	//this->cp = Air_Ptxcp ( pressure, temperature, humidity );
	this->thc = Air_Ptxthc ( pressure, temperature, humidity );
	this->visc = Air_Ptxvisc ( pressure, temperature, humidity );
	this->Psat = Air_tPs ( temperature );

	_P = pressure;
	_h = enthalpy;
	_x = humidity;

}

/********************************************************************************
 * Air::Air_tx -- 任意状態関数(温度・絶対湿度入力)
 *                 任意の温度・絶対湿度における状態の算出
 *
 * 引数(入力値)
 *     temperature  -- 温度 deg.C
 *     humidity     -- 絶対湿度 kg/kg(DA)
 ********************************************************************************/
void Air::Air_tx ( double temperature, double humidity )
{
	this->T = temperature;
	this->x = humidity;
	this->rho = Air_txrho ( temperature, humidity );
	this->h = Air_txh ( temperature, humidity );
	this->thc = Air_tthc ( temperature );
	this->visc = Air_tvisc ( temperature );
	this->Psat = Air_tPs ( temperature );

}

/********************************************************************************
 * Air::Air_hx -- 任意状態関数(比エントロピー・絶対湿度入力)
 *                 任意の比エントロピー・絶対湿度における状態の算出
 *
 * 引数(入力値)
 *     enthalpy    -- 比エンタルピー J/kg(DA)
 *     humidity     -- 絶対湿度 kg/kg(DA)
 ********************************************************************************/
void Air::Air_hx ( double enthalpy, double humidity )
{
	double temperature;
	this->T = temperature = Air_hxt ( enthalpy, humidity );
	this->x = humidity;
	this->rho = Air_txrho ( temperature, humidity );
	this->h = enthalpy;
	this->thc = Air_tthc ( temperature );
	this->visc = Air_tvisc ( temperature );
	this->Psat = Air_tPs ( temperature );

}

/***********************************************
* 湿り空気比容積関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：比容積[m^3/kg(DA)]
************************************************/
double Air_Ptxv ( double P, double T, double x )
{
	double v, v_;
	double na, nw;
	double Xa, Xw;
	double Bm, Cm;

	// 単位変換
	P *= 1e3; // kPa -> Pa
	T += 273.15; // deg.C -> K

	// モル
	na = 1.0 / Ma;
	nw = x / Mw;

	// モル分率
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// Bm, Cm計算
	Bm_Cm ( T, Xa, Xw, Bm, Cm );
	
	// 比容積計算
	// 仮定
	v = R * T / P;

	LoopCount = 0;
	while ( LoopCount < Loop )
	{
		v_ = R * T / P * ( 1.0 + Bm * 1e-6 / v + Cm * 1e-12 / ( v * v ) );
		
		// 収束判定
		if ( fabs ( 1.0 - v_ / v ) < Error ) break;

		// 値の更新
		v = v_;

		LoopCount++;
	}

	// 単位換算
	v = v_ / ( Xa * Ma * 1e-3 ); // m^3/mol -> m^3/kg(DA)

	return v;
}

/***********************************************
* 湿り空気比容積関数	-60 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       相対湿度[-]
* 出力：比容積[m^3/kg(DA)]
************************************************/
double Air_Ptphiv ( double P, double T, double phi )
{
	double x;
	double v;

	// 絶対湿度
	x = Air_Ptphix ( P, T, phi );

	// 比容積
	v = Air_Ptxv ( P, T, x );

	return ( v );
}

/***********************************************
* 湿り空気密度関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：密度[kg(DA)/m^3]
************************************************/
double Air_Ptxrho ( double P, double T, double x )
{
	double rho;

	// 密度
	rho = 1.0 / Air_Ptxv ( P, T, x );

	return ( rho );
}

/***********************************************
* 湿り空気密度関数	-60 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       相対湿度[-]
* 出力：密度[kg(DA)/m^3]
************************************************/
double Air_Ptphirho ( double P, double T, double phi )
{
	double rho;

	// 密度
	rho = 1.0 / Air_Ptphiv ( P, T, phi );

	return ( rho );
}

/***********************************************
* 湿り空気比エンタルピー関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：比エンタルピー[kJ/kg(DA)]
************************************************/
double Air_Ptxh ( double P, double T, double x )
{
	int i;

	double na, nw;
	double Xa, Xw;
	double Bm, Cm;
	double ha, hw, hm;
	double v;
	double G[6], H[6];
	double GT, HT;
	double T_delta, Bm_delta, Cm_delta;

	// 初期化
	GT = 0.0;
	HT = 0.0;

	G[0] = 0.63290874 * 1e1;
	G[1] = 0.28709015 * 1e2;
	G[2] = 0.26431805 * 1e-2;
	G[3] = -0.10405863 * 1e-4;
	G[4] = 0.18660410 * 1e-7;
	G[5] = -0.9784331 * 1e-11;

	H[0] = -0.5008 * 1e-2;
	H[1] = 0.32491829 * 1e2;
	H[2] = 0.65576345 * 1e-2;
	H[3] = -0.26442147 * 1e-4;
	H[4] = 0.517517889 * 1e-7;
	H[5] = -0.31541624 * 1e-10;

	ha = -7914.1982; // J/mol
	hw = 35994.17; // J/mol

	// モル
	na = 1.0 / Ma;
	nw = x / Mw;

	// モル分率
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// 比容積v
	v = Air_Ptxv ( P, T, x );

	// 単位換算
	v *= Xa * Ma * 1e-3 * 1e6; // m^3/kg(DA) -> cm^3/mol

	// 単位変換
	T += 273.15; // deg.C -> K

	// GT, HT
	for ( i = 0 ; i <= 5 ; i++ )
	{
		GT += G[i] * pow ( T, i );
		HT += H[i] * pow ( T, i );
	}

	// Bm, Cm
	Bm_Cm ( T, Xa, Xw, Bm, Cm );

	// Bm_delta, Cm_delta
	T_delta = T * Delta;
	Bm_Cm ( T_delta, Xa, Xw, Bm_delta, Cm_delta );

	// 比エンタルピー
	hm = Xa * ( GT + ha ) + Xw * ( HT + hw ) + R * T * ( ( Bm - T * ( Bm_delta - Bm ) / ( T_delta - T ) ) / v 
		+ ( Cm - 0.5 * T * ( Cm_delta - Cm ) / ( T_delta - T ) ) / ( v * v ) );

	// 単位換算
	hm /= ( Xa * Ma ); // J/mol -> kJ/kg(DA)

	return ( hm );
}

/***********************************************
* 湿り空気比エンタルピー関数	-60 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       相対湿度[-]
* 出力：比エンタルピー[kJ/kg(DA)]
************************************************/
double Air_Ptphih ( double P, double T, double phi )
{
	double x;
	double h;

	// 絶対湿度
	x = Air_Ptphix ( P, T, phi );

	// 比エンタルピー
	h = Air_Ptxh ( P, T, x );

	return ( h );
}

/***********************************************
* 湿り空気比エントロピー関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：比エントロピー[kJ/(kg(DA)*K)]
************************************************/
double Air_Ptxs ( double P, double T, double x )
{
	int i;

	double na, nw;
	double Xa, Xw;
	double Bm, Cm;
	double sa, sw, sm;
	double v;
	double R0;
	double M[6], N[7];
	double MT, NT;
	double T_delta, Bm_delta, Cm_delta;

	// 初期化
	MT = 0.0;
	NT = 0.0;

	M[0] = 0.34373874 * 1e2;
	M[1] = 0.52863609 * 1e-2;
	M[2] = -0.15608795 * 1e-4;
	M[3] = 0.24880547 * 1e-7;
	M[4] = -0.12230416 * 1e-10;
	M[5] = 0.28709015 * 1e2;

	N[0] = 0.2196603 * 1e1;
	N[1] = 0.19743819 * 1e-1;
	N[2] = -0.70128225 * 1e-4;
	N[3] = 0.14866252 * 1e-6;
	N[4] = -0.14524437 * 1e-9;
	N[5] = 0.55663583 * 1e-13;
	N[6] = 0.32284652 * 1e2;

	sa = -196.125465; // J/(mol･K)
	sw = -63.31449; // J/(mol･K)

	// モル
	na = 1.0 / Ma;
	nw = x / Mw;

	// モル分率
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// 比容積v
	v = Air_Ptxv ( P, T, x );

	// 単位換算
	v *= Xa * Ma * 1e-3 * 1e6; // m^3/kg(DA) -> cm^3/mol

	// 単位変換
	P *= 1e3; // kPa -> Pa
	T += 273.15; // deg.C -> K
	R0 = R * 1e6;

	// MT, NT
	for ( i = 0 ; i <= 4 ; i++ )
	{
		MT += M[i] * pow ( T, i );
		NT += N[i] * pow ( T, i );
	}
	NT += N[i] * pow ( T, i );

	// Bm, Cm
	Bm_Cm ( T, Xa, Xw, Bm, Cm );

	// Bm_delta, Cm_delta
	T_delta = T * Delta;
	Bm_Cm ( T_delta, Xa, Xw, Bm_delta, Cm_delta );

	// 比エントロピー
	sm = Xa * ( MT + M[5] * log ( T ) + sa ) + Xw * ( NT + N[6] * log ( T ) + sw ) - R * log ( P / 101325.0 ) + Xa * R * log ( ( P * v ) / ( Xa * R0 * T ) )
		+ Xw * R * log ( ( P * v ) / ( Xw * R0 * T ) ) - R * ( ( Bm + T * ( Bm_delta - Bm ) / ( T_delta - T ) ) / v
		+ 0.5 * ( Cm + T * ( Cm_delta - Cm ) / ( T_delta - T ) ) / ( v * v ) );

	// 単位換算
	sm /= ( Xa * Ma ); // J/(mol･K) -> kJ/(kg(DA)･K)

	return ( sm );
}

/***********************************************
* 湿り空気比エントロピー関数	-60 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       相対湿度[-]
* 出力：比エントロピー[kJ/(kg(DA)*K)]
************************************************/
double Air_Ptphis ( double P, double T, double phi )
{
	double x;
	double s;

	// 絶対湿度
	x = Air_Ptphix ( P, T, phi );

	// 比エントロピー
	s = Air_Ptxs ( P, T, x );

	return ( s );
}

/***********************************************
* 湿り空気相対湿度関数	-60 <= T <= 365
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：相対湿度[-]
************************************************/
double Air_Ptxphi ( double P, double T, double x )
{
	double phi;
	double Pvs;
	double f;
	double na, nw;
	double Xw;

	// モル
	na = 1.0 / Ma;
	nw = x / Mw;

	// モル分率
	Xw = nw / ( na + nw );

	// 飽和水蒸気圧
	Pvs = Air_tPs ( T );

	// f
	f = f_ ( P, T );

	// 相対湿度
	phi = Xw * P / ( f * Pvs );

	return ( phi );
}

/***********************************************
* 湿り空気相対湿度関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       湿球温度[deg.C]
* 出力：相対湿度[-]
************************************************/
double Air_Pttwbphi ( double P, double T, double Twb )
{
	double phi;
	double Twb_;
	double dT;

	int i, j;
	double X[2], Y[2];

	// 温度差
	dT = T - Twb;

	if ( dT < 0.0 )
	{
		std::cout << std::endl;
		std::cout << "Twb exceeds T" << std::endl;
		exit(1);
	}

	// 仮定値
	if ( dT > 10.0 ) phi = 0.30;
	else if ( dT > 9.0 ) phi = 0.40;
	else if ( dT > 8.0 ) phi = 0.45;
	else if ( dT > 7.0 ) phi = 0.50;
	else if ( dT > 6.0 ) phi = 0.55;
	else if ( dT > 5.0 ) phi = 0.65;
	else if ( dT > 4.0 ) phi = 0.70;
	else if ( dT > 3.0 ) phi = 0.75;
	else if ( dT > 2.0 ) phi = 0.85;
	else if ( dT > 1.0 ) phi = 0.90;
	else phi = 1.0;
	
	for ( i = 1 ; i <= Loop ; i++ )
	{
		for ( j = 0, phi *= Delta ; j < 2 ; j++ )
		{
			if ( j == 1 ) phi /= Delta;

			// 湿球温度
			Twb_ = Air_Ptphitwb ( P, T, phi );

			X[j] = phi;
			Y[j] = 1.0 - Twb_ / Twb;
		} // end of j loop
		
		// 収束判定
		if ( fabs ( Y[1] ) < Error ) break;
		else phi = X[1] - Y[1] * ( X[0] - X[1] ) / ( Y[0] - Y[1] );

		// 繰り返し計算回数
		if ( i == Loop )
		{
			std::cout << std::endl;
			std::cout << "Error -> Air_Pttwbphi" << std::endl;
			exit(1);
		}
	} // end of i loop

	return ( phi );
}

/***********************************************
* 湿り空気絶対湿度関数	-60 <= T <= 365
* 入力：圧力[kPa]
*       温度[deg.C]
*       相対湿度[-]
* 出力：絶対湿度[kg/kg(DA)]
************************************************/
double Air_Ptphix ( double P, double T, double phi )
{
	double x;
	double Pvs;
	double f;
	double na;
	double Xw;

	// モル
	na = 1.0 / Ma;//空気1kgあたりのモル数 空気1molで約28kgくらい（ソースの上方に、物性値として入力）

	// 飽和水蒸気圧  M.CONDE ENGINEERING, 2007	の式[10]より
	Pvs = Air_tPs ( T );//

	// f
	f = f_ ( P, T );

	// モル分率 M.CONDE ENGINEERING, 2007	の式[9]より
	Xw = phi * f * Pvs / P;

	if ( Xw > 1.0 )
	{
		std::cout << std::endl;
		std::cout << "Pvs exceeds P -> Air_Ptphix" << std::endl;
		exit(1);
	}

	// 絶対湿度
	x = na * Xw * Mw / ( 1.0 - Xw );

	return ( x );
}

/***********************************************
* 湿り空気絶対湿度関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       比エンタルピー[kJ/kg(DA)]
*       温度[deg.C]
* 出力：絶対湿度[kg/kg(DA)]
************************************************/
double Air_Phtx ( double P, double h, double T )
{
	double x;
	double T_;
	int i, j;
	double X[2], Y[2];

	// 仮定
	x = ( h - 1.006 * T ) / ( 1.805 * T + 2501.0 );

	for ( i = 1 ; i <= Loop ; i++ )
	{
		for ( j = 0, x *= Delta ; j < 2 ; j++ )
		{
			if ( j == 1 )	x /= Delta;

			// 温度
			T_ = Air_Phxt ( P, h, x );

			X[j] = x;
			Y[j] = T - T_;
		} // end of j loop

		// 収束判定
		if ( fabs ( Y[1] ) < Error ) break;
		else x = X[1] - Y[1] * ( X[0] - X[1] ) / ( Y[0] - Y[1] );

		// 繰り返し計算回数
		if ( i == Loop )
		{
			std::cout << std::endl;
			std::cout << "Error -> Air_Phxt1" << std::endl;
			exit(1);
		}
	} // end of i loop

	return ( x );
}

/***********************************************
* 湿り空気絶対湿度関数	-60 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       湿球温度[deg.C]
* 出力：絶対湿度[kg/kg(DA)]
************************************************/
double Air_Pttwbx ( double P, double T, double Twb )
{
	double x;
	double phi;

	// 相対湿度
	phi = Air_Pttwbphi ( P, T, Twb );

	// 絶対湿度
	x = Air_Ptphix ( P, T, phi );

	return ( x );
}

/***********************************************
* 湿り空気温度関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       比エンタルピー[kJ/kg(DA)]
*       絶対湿度[kg/kg(DA)]
* 出力：温度[deg.C]
************************************************/
double Air_Phxt ( double P, double h, double x )
{
	double T;
	double h_;

	int i, j;
	double X[2], Y[2];

	// 仮定値
	T = ( h - 2501.0 * x ) / ( 1.006 + 1.805 * x );

	if( - 0.001 < T && T < 0.001 ){
		T = 0.001;
	}//仮定値が0になったときに計算できないバグ修正　大野

	for ( i = 1 ; i <= Loop ; i++ )
	{
		for ( j = 0, T *= Delta ; j < 2 ; j++ )
		{
			if ( j == 1 )	T /= Delta;

			// 比エンタルピー
			h_ = Air_Ptxh ( P, T, x );

			X[j] = T;
			Y[j] = h - h_;
		} // end of j loop
	
		// 収束判定
		if ( fabs ( Y[1] ) < Error ) break;
		else T = X[1] - Y[1] * ( X[0] - X[1] ) / ( Y[0] - Y[1] );

		// 繰り返し計算回数
		if ( i == Loop )
		{
			std::cout << std::endl;
			std::cout << "Error -> Air_Phxt2" << std::endl;
			exit(1);
		}
	} // end of i loop

	return ( T );
}

/***********************************************
* 湿り空気温度関数	-60 <= T <= 200
* 入力：圧力[kPa]
*       比エンタルピー[kJ/kg(DA)]
*       相対湿度[-]
* 出力：温度[deg.C]
************************************************/
double Air_Phphit ( double P, double h, double phi )
{
	double T;
	double h_;

	int i, j;
	double X[2], Y[2];

	// 仮定値
	if ( h < -18.537 ) T = -20.0;
	else if ( h < -6.0866 ) T = -10.0;
	else if ( h < -0.3758 ) T = -5.0;
	else if ( h < 5.9637 ) T = -2.0;
	else if ( h < 9.4353 ) T = 0.01;
	else if ( h < 29.261 ) T = 10.0;
	else if ( h < 57.358 ) T = 20.0;
	else if ( h < 99.602 ) T = 30.0;
	else if ( h < 165.87 ) T = 40.0;
	else if ( h < 273.71 ) T = 50.0;
	else if ( h < 457.54 ) T = 60.0;
	else if ( h < 796.33 ) T = 70.0;
	else if ( h < 1524.8 ) T = 80.0;
	else if ( h < 3812.0 ) T = 90.0;
	else if ( h < 6908.2 ) T = 96.0;
	else T = 98.0;

	for ( i = 1 ; i <= Loop ; i++ )
	{
		for ( j = 0, T *= Delta ; j < 2 ; j++ )
		{
			if ( j == 1 ) T /= Delta;

			// 比エンタルピー
			h_ = Air_Ptphih ( P, T, phi );

			X[j] = T;
			Y[j] = h - h_;
		} // end of j loop
		
		// 収束判定
		if ( fabs ( Y[1] ) < Error ) break;
		else T = X[1] - Y[1] * ( X[0] - X[1] ) / ( Y[0] - Y[1] ) * 0.8;

		// 繰り返し計算回数
		if ( i == Loop )
		{
			std::cout << std::endl;
			std::cout << "Error -> Air_Phxt3" << std::endl;
			exit(1);
		}
	}

	return ( T );
}

/***********************************************
* 飽和水蒸気圧関数	-60 <= T <= 365
* 入力：温度[deg.C]圧力[kPa]
* 出力：圧力[kPa]
************************************************/
double Air_tPs ( double T )
{
	double Pvs;//M.CONDE ENGINEERING, 2007の6ページ目のMole Fraction and Humidity Ratioより Psv:飽和水蒸気圧
	double t, s;
	double Pcr, Tcr;
	double A, A1, A2, A3, A4, A5, A6;

	Pcr = 22064.0; // kPa M.CONDE ENGINEERING, 2007の6ページ目のMole Fraction and Humidity Ratioの式より Pcr,H2O
	Tcr = 647.14; // K  M.CONDE ENGINEERING, 2007の6ページ目のMole Fraction and Humidity Ratioの式より Tcr,H2O

	// 単位換算
	T += 273.15;

	// θ，τ
	s = T / Tcr;//θ=T/Tcr,H2O
	t = 1.0 - s;//τ=1-θ

	if ( T >= 213.15 && T < 273.16 )
	{
		A1 = 2.442663;
		A2 = -11.413077;
		A3 = -15.109346;
		A4 = 1.119193;
		A5 = 18.159568;
		A6 = -6.138264;
	}
	else if ( T >= 273.16 && T <= 647.14 )
	{
		A1 = -7.858230;
		A2 = 1.839910;
		A3 = -11.781100;
		A4 = 22.670500;
		A5 = -15.939300;
		A6 = 1.775160;
	}
	else 
	{
		std::cout << std::endl;
		std::cout << "Over the T range -> Air_tPs" << std::endl;
		std::cout << T << std::endl;
		exit(1);
	}

	// 飽和水蒸気圧
	A = A1 * t + A2 * pow ( t, 1.5 ) + A3 * pow ( t, 3.0 ) + A4 * pow ( t, 3.5 ) + A5 * pow ( t, 4.0 ) + A6 * pow ( t, 7.5 );//M.CONDE ENGINEERING, 2007の6ページ目の式[10]よりln(Psv/Pcr,H2O)=1/θ*(A1*τ+A2*τ^1.5+A3*τ^3+A4*τ^3.5+A5*τ^4+A6*τ^7.5)  
	Pvs = Pcr * exp ( A / s );

	return ( Pvs );
}

/***********************************************
* 湿り空気定圧比熱関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：定圧比熱[kJ/(kg(DA)*K)]
************************************************/
double Air_Ptxcp ( double P, double T, double x )
{
	double cp;
	double h1, h2;
	const double dT = 0.05;

	// h1, h2
	h1 = Air_Ptxh ( P, T - dT, x );
	h2 = Air_Ptxh ( P, T + dT, x );

	// 定圧比熱
	cp = ( h2 - h1 ) / ( 2 * dT );

	return ( cp );
}

/***********************************************
* 湿り空気熱伝導率関数	0.01 <= T <= 350
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：熱伝導率[kW/(m*K)]
************************************************/
double Air_Ptxthc ( double P, double T, double x )
{
	double thc;
	double Pa, Pw;
	double f;
	double thc_a, thc_w;
	double Gaw, Gwa;
	double na, nw;
	double Xa, Xw;
	double visc_a, visc_w;

	// モル
	na = 1.0 / Ma;
	nw = x / Mw;

	// モル分率
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// f
	f = f_ ( P, T );

	// 分圧
	Pw = P * Xw / f;
	Pa = P - Pw;

	// 乾き空気粘性係数，熱伝達率
	visc_a = Dair_Ptvisc ( Pa, T );
	thc_a = Dair_Ptthc ( Pa, T );

	if ( x == 0.0 )	thc = thc_a;
	else
	{
		// 水蒸気粘性係数，熱伝達率
		visc_w = sh_myuv ( Pw, T );
		thc_w = sh_lamv ( Pw, T );

		// Gaw, Gwa
		Gaw = sqrt ( 2.0 ) / 4.0 * sqrt ( Mw / ( Ma + Mw ) ) * pow ( 1.0 + sqrt ( visc_a / visc_w ) * pow ( Mw / Ma, 0.25 ), 2.0 );
		Gwa = sqrt ( 2.0 ) / 4.0 * sqrt ( Ma / ( Ma + Mw ) ) * pow ( 1.0 + sqrt ( visc_w / visc_a ) * pow ( Ma / Mw, 0.25 ), 2.0 );

		// 熱伝導率
		thc = thc_a / ( 1.0 + Gaw * Xw / Xa ) + thc_w / ( 1.0 + Gwa * Xa / Xw );
	}

	// 単位換算
	thc *= 1e-3; // W/(m*K) -> kW/(m*K)

	return ( thc );
}

/***********************************************
* 湿り空気粘性係数関数	0.01 <= T <= 350
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：粘性係数[Pa*s]
************************************************/
double Air_Ptxvisc ( double P, double T, double x )
{
	double visc;
	double Pa, Pw;
	double f;
	double Gaw, Gwa;
	double na, nw;
	double Xa, Xw;
	double visc_a, visc_w;

	// モル
	na = 1.0 / Ma;
	nw = x / Mw;

	// モル分率
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// f
	f = f_ ( P, T );

	// 分圧
	Pw = P * Xw / f;
	Pa = P - Pw;

	// 乾き空気粘性係数
	visc_a = Dair_Ptvisc ( Pa, T );

	if ( x == 0.0 )	visc = visc_a;
	else
	{
		// 水蒸気粘性係数，熱伝達率
		visc_w = sh_myuv ( Pw, T );

		// Gaw, Gwa
		Gaw = sqrt ( 2.0 ) / 4.0 * sqrt ( Mw / ( Ma + Mw ) ) * pow ( 1.0 + sqrt ( visc_a / visc_w ) * pow ( Mw / Ma, 0.25 ), 2.0 );
		Gwa = sqrt ( 2.0 ) / 4.0 * sqrt ( Ma / ( Ma + Mw ) ) * pow ( 1.0 + sqrt ( visc_w / visc_a ) * pow ( Ma / Mw, 0.25 ), 2.0 );

		// 粘性係数
		visc = visc_a / ( 1.0 + Gaw * Xw / Xa ) + visc_w / ( 1.0 + Gwa * Xa / Xw );
	}

	return ( visc );
}

/***********************************************
* 湿り空気動粘性係数関数	0.01 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：動粘性係数[m^2/s]
************************************************/
double Air_Ptxkine ( double P, double T, double x )
{
	double kine;
	double visc, rho;

	// 密度
	rho = Air_Ptxrho ( P, T, x );

	// 粘性係数
	visc = Air_Ptxvisc ( P, T, x );

	// 動粘性係数
	kine = visc / rho;

	return ( kine );
}

/***********************************************
* 湿り空気熱拡散係数関数	0.01 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：熱拡散係数[m^2/s]
************************************************/
double Air_Ptxthd ( double P, double T, double x )
{
	double thd;
	double rho, cp, thc;

	// 密度
	rho = Air_Ptxrho ( P, T, x );

	// 比熱
	cp = Air_Ptxcp ( P, T, x );

	// 熱伝導率
	thc = Air_Ptxthc ( P, T, x );

	thd = thc / ( rho * cp );

	return ( thd );
}

/***********************************************
* 湿り空気水蒸気の拡散係数関数	-20 <= T <= 300
* 入力：圧力[kPa]
*       温度[deg.C]
* 出力：水蒸気の拡散係数[m^2/s]
************************************************/
double Air_PtDh2o ( double P, double T )
{
	double D;

	// 単位換算
	T += 273.15; // deg.C -> K
	P *= 1e3; // kPa -> Pa

	// 場合分け
	if ( T >= 253.15 && T <= 353.15 )	D = 104.91143 * 1e-6 * pow ( T, 1.774 ) / P;

	else if ( T > 353.15 && T <= 573.15 )	D = 805.2375 * 1e-6 / P * pow ( T, 2.5 ) / ( T + 190 );

	else
	{
		std::cout << std::endl;
		std::cout << "over the temperature range -> Air_PtDh2o" << std::endl;
		std::cout << T << std::endl;
		exit (1);
	}

	return ( D );
}

/***********************************************
* 湿り空気プロントル数関数	0.01 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：プラントル数[-]
************************************************/
double Air_PtxPr ( double P, double T, double x )
{
	double Pr;
	double thd, kine;

	// 熱拡散係数
	thd = Air_Ptxthd ( P, T, x );

	// 動粘性係数
	kine = Air_Ptxkine ( P, T, x );

	// プラントル数
	Pr = kine / thd;

	return ( Pr );
}

/***********************************************
* 湿り空気シュミット数関数	0.01 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：シュミット数[-]
************************************************/
double Air_PtxSc ( double P, double T, double x )
{
	double Sc;
	double D, kine;

	// 拡散係数
	D = Air_PtDh2o ( P, T );

	// 動粘性係数
	kine = Air_Ptxkine ( P, T, x );

	// シュミット数
	Sc = kine / D;

	return ( Sc );
}

/***********************************************
* 湿り空気露点温度関数	x >= 0.0
* 入力：圧力[kPa]
*       絶対湿度[kg/kg(DA)]
* 出力：温度[deg.C]
************************************************/
double Air_Pxtdp ( double P, double x )
{
	double x_;
	double Tdp;
	double TL, TR, dT;

	// 初期刻み幅
	dT = 1.0;

	// 露点温度仮定値
	if(x<0.0000067) Tdp=-60.0;
	else if(x<0.0000243) Tdp=-50.0;
	else if(x<0.0000793) Tdp=-40.0;
	else if(x<0.0002346) Tdp=-30.0;
	else if(x<0.0006373) Tdp=-20.0;
	else if(x<0.0016061) Tdp=-10.0;
	else if(x<0.0037896) Tdp=0.01;
	else if(x<0.0076609) Tdp=10.0;
	else if(x<0.0147571) Tdp=20.0;
	else if(x<0.0273278) Tdp=30.0;
	else if(x<0.0491377) Tdp=40.0;
	else if(x<0.0868538) Tdp=50.0;
	else if(x<0.1535323) Tdp=60.0;
	else if(x<0.2791495) Tdp=70.0;
	else if(x<0.5528998) Tdp=80.0;
	else if(x<1.4291908) Tdp=90.0;
	else if(x<3.1943411) Tdp=95.0;
	else Tdp=99.0;

	// 露点温度
	if(x==0.0)	Tdp = -100.0;
	else
	{
		// 仮定露点絶対湿度
		x_ = Air_Ptphix ( P, Tdp, 1.0 );

		LoopCount = 0;
		while ( x_ >= x || LoopCount < Loop )
		{
			// 露点温度仮定
			Tdp -= dT;

			// 仮定露点絶対湿度
			x_ = Air_Ptphix ( P, Tdp, 1.0 );

			LoopCount++;
		}

		TL = Tdp;
		TR = Tdp + dT;

		LoopCount = 0;
		while ( LoopCount < Loop )
		{
			// 露点温度仮定
			Tdp = ( TL + TR ) / 2.0;

			// 仮定露点絶対湿度
			x_ = Air_Ptphix ( P, Tdp, 1.0 );

			// 収束判定
			if ( fabs ( 1.0 - x_ / x ) < Error ) break;

			// 値の更新
			if ( x_ > x ) TR = Tdp;
			else TL = Tdp;

			LoopCount++;
		}
	}

	return ( Tdp );
}

/***********************************************
* 湿り空気凝固点温度関数	P <= 200.0[MPa]
* 入力：圧力[kPa]
* 出力：温度[deg.C]
************************************************/
double Air_Ptmel ( double P )
{
	double Tm;
	double Ptp, Ttp;

	Ptp = 0.006112; // bar
	Ttp = 273.16; // K

	// 単位換算
	P /= 1e8; // kPa -> bar

	// 凝固点温度
	Tm = Ttp * pow ( - ( P - Ptp ) / 3950 + 1.0, 1.0 / 9.0 );

	// 単位換算
	Tm -= 273.15; // K -> deg.C

	return ( Tm );
}



// 下記の関数は上記の関数内で使用する関数

/***********************************************
* 湿り空気湿球温度関数	-60 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       絶対湿度[kg/kg(DA)]
* 出力：温度[deg.C]
************************************************/
double Air_Ptxtwb ( double P, double T, double x )
{
	double Twb;
	double phi;

	// 相対湿度
	phi = Air_Ptxphi ( P, T, x );

	// 湿球温度
	Twb = Air_Ptphitwb ( P, T, phi );

	return ( Twb );
}

/***********************************************
* 湿り空気湿球温度関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       相対湿度[-]
* 出力：温度[deg.C]
************************************************/
double Air_Ptphitwb ( double P, double T, double phi )
{
	double Twb;
	double x, xwb;
	double h, hwb, hc, h_;

	double TL, TR, dT;

	// 初期刻み幅
	dT = 1.0;

	// 湿球温度
	if ( phi == 1.0 )	Twb = T;
	else
	{
		// 湿球温度仮定
		Twb = T;

		// 絶対湿度
		x = Air_Ptphix ( P, T, phi );

		// 比エンタルピー
		h = Air_Ptxh ( P, T, x );

		// 飽和絶対湿度
		xwb = Air_Ptphix ( P, Twb, 1.0 );

		// 飽和比エンタルピー
		hwb = Air_Ptxh ( P, Twb, xwb );

		// 凝縮相比エンタルピー
		hc = Air_Pthc ( P, Twb );

		LoopCount = 0;
		while ( h < hwb - ( xwb - x ) * hc  )
		{
			// 湿球温度仮定
			Twb -= dT;

			// 飽和絶対湿度
			xwb = Air_Ptphix ( P, Twb, 1.0 );

			// 飽和比エンタルピー
			hwb = Air_Ptxh ( P, Twb, xwb );

			// 凝縮相比エンタルピー
			hc = Air_Pthc ( P, Twb );

			LoopCount++;
			if( LoopCount > Loop ){
				std::cout << "Air_Ptphitwb loop over1";
			}
		}

		TL = Twb - dT;
		TR = Twb + dT;

		LoopCount = 0;
		while ( LoopCount < Loop )
		{
			// 湿球温度仮定
			Twb = ( TL + TR ) / 2.0;

			// 飽和絶対湿度
			xwb = Air_Ptphix ( P, Twb, 1.0 );

			// 飽和比エンタルピー
			hwb = Air_Ptxh ( P, Twb, xwb );

			// 凝縮相比エンタルピー
			hc = Air_Pthc ( P, Twb );

			// 比エンタルピー
			h_ = hwb - ( xwb - x ) * hc;

			// 収束判定
			if ( fabs ( 1.0 - h_ / h ) < Error && fabs ( 1.0 - TL / Twb ) < Error ) break;

			// 値の更新
			if ( h_ > h ) TR = Twb;
			else TL = Twb;

			LoopCount++;
		}

	}

	return ( Twb );
}

/***********************************************
* 湿り空気湿球温度関数	-100 <= T <= 200
* 入力：圧力[kPa]
*       温度[deg.C]
*       相対湿度[-]
* 出力：温度[deg.C]
************************************************/
double Air_Pthc ( double P, double T )
{
	double h;
	double Psv;
	double T_delta, Psv_delta;
	double Tm;
	double v;
	double z[5], l[7], o[6], g[8];
	double zT, lT, oT, gT1, gT2; 
	double a, b;

	int i;

	// 初期化
	zT = 0.0;
	lT = 0.0;
	oT = 0.0;
	gT1 = gT2 = 0.0;

	z[0] = -0.647595 * 1e3;
	z[1] = 0.274292;
	z[2] = 0.2910583 * 1e-2;
	z[3] = 0.1083437 * 1e-5;
	z[4] = 0.107 * 1e-5;

	l[0] = -0.11411380 * 1e4;
	l[1] = 0.41930463 * 1e1;
	l[2] = -0.8134865 * 1e-4;
	l[3] = 0.1451133 * 1e-6;
	l[4] = -0.1005230 * 1e-9;
	l[5] = -0.563473;
	l[6] = -0.036;

	o[0] = -0.1141837121 * 1e4;
	o[1] = 0.4194325677 * 1e1;
	o[2] = -0.6908894163 * 1e-4;
	o[3] = 0.105555302 * 1e-6;
	o[4] = -0.7111382234 * 1e-10;
	o[5] = 0.6059 * 1e-6;

	g[0] = -0.2403360201 * 1e4;
	g[1] = -0.140758895 * 1e1;
	g[2] = 0.1068287657;
	g[3] = -0.2914492351 * 1e-3;
	g[4] = 0.373497936 * 1e-6;
	g[5] = -0.21203787 * 1e-9;
	g[6] = -0.3424442728 * 1e1;
	g[7] = 0.1619785 * 1e-1;

	// 飽和水蒸気圧
	Psv = Air_tPs ( T );

	// 凝固点温度
	Tm = Air_Ptmel ( P );

	// 凝縮相（固体）
	if ( T < Tm )
	{
		// 単位換算
		T += 273.15;

		for ( i = 0 ; i <= 3 ; i++ )	zT += z[i] * pow ( T, i );

		// 比エンタルピー
		h = zT + z[4] * Psv;
	}
	// 凝縮相（液体）
	else
	{
		// 飽和水蒸気圧
		T_delta = T * Delta;
		Psv_delta = Air_tPs ( T_delta );

		// 単位換算
		T += 273.15; // deg.C -> K

		// 場合分け
		if ( T <= 373.15 )
		{
			for ( i = 0 ; i <= 4 ; i++ ) lT += l[i] * pow ( T, i );

			// α
			a = lT + l[5] * pow ( 10.0, l[6] * ( T - 273.15 ) );
		}
		else if ( T > 373.15 && T <= 403.128 )
		{
			for ( i = 0 ; i <= 4 ; i++ ) oT += o[i] * pow ( T, i );

			// α
			a = oT;
		}
		else if ( T > 403.128 && T <= 473.15 )
		{
			for ( i = 0 ; i <= 4 ; i++ ) oT += o[i] * pow ( T, i );

			// α
			a = oT - o[5] * pow ( T - 403.128, 3.1 );
		}

		// 比容積
		for ( i = 0 ; i <= 5 ; i++ )	gT1 += g[i] * pow ( T, i );
		for ( i = 6 ; i <= 7 ; i++ )	gT2 += g[i] * pow ( T, i - 6 );
		v = 18015.28 * gT2 / gT1;

		// 単位換算
		v *= 1e-6; // cm^3/kg -> m^3/kg

		// β
		b = T * v * ( Psv_delta - Psv ) / ( T_delta - ( T - 273.15 ) );

		// 比エンタルピー
		h = a - 0.01214 + b;
	}
	
	return ( h );
}

/***********************************************
* Bm, Cm計算用関数
* 入力：温度[K]
*       モル比[-]
* 出力：Bm[cm^3/mol]
*       Cm[cm^6/mol^2]
************************************************/
void Bm_Cm ( double T, double Xa, double Xw, double &Bm, double &Cm )
{
	int i;

	double Baa, Baw, Bww;
	double Caaa, Caaw, Caww, Cwww;
	double B1[4], C1[3];
	double B2[3], C2[3];
	double D[5], E[5], F[5];

	// 初期化
	Baa = Baw = Bww = 0.0;
	Caaa = Caaw = Caww = Cwww = 0.0;

	//EOS for Air 実在気体に補正するため，ビリアル方程式　　Z=PVa/(RT)=1+Baa/Va+Caaa/(Va^2)+……に合うようにする．Vw:1molあたりの体積(mol体積)．ビリアル係数は第三ビリアル係数以下は省略している．
	//B1[i]:M.CONDE ENGINEERING, 2007の2ページの式[1]に入れる係数  Baa:[cm^3/mol]=Σ(i=0～3)B1[i]*T^(-i)
	B1[0] = 0.349568 * 1e2;
	B1[1] = -0.668772 * 1e4;
	B1[2] = -0.210141 * 1e7;
	B1[3] = 0.924746 * 1e8;
	//C1[i]:M.CONDE ENGINEERING, 2007の2ページの式[1]に入れる係数 Caaa:[cm^6/(mol^2)]=Σ(i=0～3)C1[i]*T^(-i)  C1[3]はM.CONDE ENGINEERING, 2007の2ページの表には書かれていない……
	C1[0] = 0.125975 * 1e4;
	C1[1] = -0.190905 * 1e6;
	C1[2] = 0.632467 * 1e8;

	//EOS for Water Vapour 実在気体に補正するため，ビリアル方程式　　Z=PVw/(RT)=1+Bww/Vw+Cwww/(Vw^2)+……に合うようにする．Vw:1molあたりの体積(mol体積)．ビリアル係数は第三ビリアル係数以下は省略している．
	//B2[i]:M.CONDE ENGINEERING, 2007の2ページの式[2]に入れる係数  Bww:[cm^3/mol]=RT*(B2[0]+B2[1]*e^(B2[2]/T))
	B2[0] = 0.70 * 1e-8;
	B2[1] = -0.147184 * 1e-8;
	B2[2] = 1734.29;
	//C2[i]:M.CONDE ENGINEERING, 2007の2ページの式[2]に入れる係数  Cwww:[cm^6/(mol^2)]=((RT)^2)*((C2[0]+C2[1]*e^(C2[2]/T))^2)*Bww^2
	C2[0] = 0.104 * 1e-14;
	C2[1] = -0.335297 * 1e-17;
	C2[2] = 3645.09;

	//EOS for the Mixture(humid air) 実在気体に補正するため，ビリアル方程式　　Z=P*平均Vm/(RT)=1+Bm/平均Vm+Cm/(平均Vm^2)+……に合うようにする．平均Vm:1molあたりの体積(mol体積)．ビリアル係数は第三ビリアル係数以下は省略している．
	//Bm[cm^3/mol]=(Xa^2)*Baa+2*Xa*Xw*Baw+(Xw^2)*Bww  M.CONDE ENGINEERING, 2007の3ページの式[3]
	//Baw[cm^3/mol]=Σ(i=0～4)D[i]*(T^(-i))  M.CONDE ENGINEERING, 2007の3ページの式[5]
	//D[i]:M.CONDE ENGINEERING, 2007の2ページの式[2]に入れる係数  Bww:[cm^3/mol]=RT*(B2[0]+B2[1]*e^(B2[2]/T))
	D[0] = 0.32366097 * 1e2;
	D[1] = -0.141138 * 1e5;
	D[2] = -0.1244535 * 1e7;
	D[3] = 0.0;
	D[4] = -0.2348789 * 1e10;

	//Cm[cm^6/(mol^2)]=Xa^3*Caaa+3*(Xa^2)*Xw*Caaw+3*Xa*(Xw^2)*Caww+(Xw^3)*Cwww  M.CONDE ENGINEERING, 2007の3ページの式[3]
	//Caaw[(cm^6)/(mol^2)]=Σ(i=0～4)E[i]*(T^-i)
	E[0] = 0.482737 * 1e3;
	E[1] = 0.105678 * 1e6;
	E[2] = -0.656394 * 1e8;
	E[3] = 0.299444 * 1e11;
	E[4] = -0.319317 * 1e13;
	//Caww[(cm^6)/(mol^2)]=F4*Σ(i=0～4)E[i]*(T^-i)
	F[0] = -0.10728876 * 1e2;
	F[1] = 0.347802 * 1e4;
	F[2] = -0.383383 * 1e6;
	F[3] = 0.33406 * 1e8;
	F[4] = -1.0 * 1e-6;

	// Bww, Cwww
	Bww = R * 1e6 * T * ( B2[0] + B2[1] * exp ( ( B2[2] / T ) ) );
	Cwww = pow ( R * 1e6 * T, 2 ) * pow ( C2[0] + C2[1] * exp ( C2[2] / T ), 2 ) * ( Bww * Bww );

	// Caaa
	for ( i = 0 ; i <= 2 ; i++ )	Caaa += C1[i] * pow ( T, -i );

	// Baa, Caww
	for ( i = 0 ; i <= 3 ; i++ )
	{
		Baa += B1[i] * pow ( T, -i );
		Caww += F[4] * exp ( ( F[i] * pow ( T, -i ) ) );
	}

	// Baw, Caaw
	for ( i = 0 ; i <= 4 ; i++ )
	{
		Baw += D[i] * pow ( T, -i );
		Caaw += E[i] * pow ( T, -i );
	}

	// Bm, Cm
	Bm = Xa * Xa * Baa + 2.0 * Xa * Xw * Baw + Xw * Xw * Bww;
	Cm = Xa * Xa * Xa * Caaa + 3.0 * Xa * Xa * Xw * Caaw + 3.0 * Xa * Xw * Xw * Caww + Xw * Xw * Xw * Cwww;
}

/***********************************************
* 修正係数f算出関数	 M.CONDE ENGINEERING, 2007の7ページ式[12]
* 入力：圧力[kPa]
*       温度[deg.C]
* 出力：修正係数[-]
************************************************/
double f_ ( double P, double T )
{
	int i;

	double f, f_, f1, f2, f3;
	double R0;

	double Xas;
	double Pvs;

	double v;
	double Tm;
	double g[8], s[3];
	double gT1, gT2, sT;

	double K, Ka, Ko, Kn, Kmax_o, Kmax_n;
	double So, Sn;
	double xSo, xSn;
	double x[5];

	double Baa, Baw, Bww;
	double Caaa, Caaw, Caww, Cwww;
	double B1[4], C1[3];
	double B2[3], C2[3];
	double D[5], E[5], F[5];

	// 初期化
	gT1 = gT2 = 0.0;
	sT = 0.0;
	xSo = 0.0;
	xSn = 0.0;
	Baa = Baw = Bww = 0.0;
	Caaa = Caaw = Caww = Cwww = 0.0;
	//g[i]:liquid water 凝縮過程の液相の比容積[cm^3/mol]の計算に用いる係数 M.CONDE ENGINEERING, 2007の7ページ目の下の表より
	g[0] = -0.2403360201 * 1e4;
	g[1] = -0.140758895 * 1e1;
	g[2] = 0.1068287657;
	g[3] = -0.2914492351 * 1e-3;
	g[4] = 0.373497936 * 1e-6;
	g[5] = -0.21203787 * 1e-9;
	g[6] = -0.3424442728 * 1e1;
	g[7] = 0.1619785 * 1e-1;
	//s[i]:solid water 凝縮過程の固相の比容積[cm^3/mol]の計算に用いる係数 M.CONDE ENGINEERING, 2007の7ページ目の下の表より
	s[0] = 0.1070003 * 1e-2;
	s[1] = -0.249936 * 1e-7;
	s[2] = 0.371611 * 1e-9;
	//Kaの値には空気の主要な構成要素である酸素と窒素のKを用いる．酸素のKoと窒素のKnの計算に用いる二つの物性値の入力 Ka=Ko*Kn/(0.22*Ko+0.78*Kn)*101325.16  M.CONDE ENGINEERING, 2007の9ページ目の式[18]
	Kmax_o = 7.08 * 1e4;//M.CONDE ENGINEERING, 2007の9ページ目の表よりKmax[atm/mol fraction](oxygen)=7.08*10^4
	Kmax_n = 12.39 * 1e4;//M.CONDE ENGINEERING, 2007の9ページ目の表よりKmax[atm/mol fraction](nitrogen)=12.39*10^4
	//x[i]: Himmelblauが解析し構築した式に用いる係数 K=Kmax*10^(Σ(i=0～4)x[i]*θ^(-i) )　M.CONDE ENGINEERING, 2007の9ページ目の式[19]
	x[0] = -1.142;
	x[1] = 2.846;
	x[2] = -2.486;
	x[3] = 0.9761;
	x[4] = -0.2001;

	B1[0] = 0.349568 * 1e2;
	B1[1] = -0.668772 * 1e4;
	B1[2] = -0.210141 * 1e7;
	B1[3] = 0.924746 * 1e8;

	C1[0] = 0.125975 * 1e4;
	C1[1] = -0.190905 * 1e6;
	C1[2] = 0.632467 * 1e8;

	B2[0] = 0.70 * 1e-8;
	B2[1] = -0.147184 * 1e-8;
	B2[2] = 1734.29;

	C2[0] = 0.104 * 1e-14;
	C2[1] = -0.335297 * 1e-17;
	C2[2] = 3645.09;

	D[0] = 0.32366097 * 1e2;
	D[1] = -0.141138 * 1e5;
	D[2] = -0.1244535 * 1e7;
	D[3] = 0.0;
	D[4] = -0.2348789 * 1e10;

	E[0] = 0.482737 * 1e3;
	E[1] = 0.105678 * 1e6;
	E[2] = -0.656394 * 1e8;
	E[3] = 0.299444 * 1e11;
	E[4] = -0.319317 * 1e13;

	F[0] = -0.10728876 * 1e2;
	F[1] = 0.347802 * 1e4;
	F[2] = -0.383383 * 1e6;
	F[3] = 0.33406 * 1e8;
	F[4] = -1.0 * 1e-6;

	// 飽和蒸気圧
	Pvs = Air_tPs ( T );

	// 凝固点温度
	Tm = Air_Ptmel ( P );

	// 単位換算
	P *= 1e3; // kPa -> Pa
	Pvs *= 1e3; // kPa -> Pa
	T += 273.15; // deg.C -> K
	Tm += 273.15; // deg.C -> K
	R0 = R * 1e6;

	// 凝縮相（液体）M.CONDE ENGINEERING, 2007の7ページ式[13]
	if ( T > Tm )
	{
		for ( i = 0 ; i <= 5 ; i++ )	gT1 += g[i] * pow ( T, i );
		for ( i = 6 ; i <= 7 ; i++ )	gT2 += g[i] * pow ( T, i - 6 );

		// 比容積[cm^3/mol]
		v = 18015.28 * gT2 / gT1;
	}
	// 凝縮相（固体）M.CONDE ENGINEERING, 2007の7ページ式[14]
	else
	{
		for ( i = 0 ; i <= 2 ; i++ )	sT += s[i] * pow ( T, i );

		// 比容積[cm^3/mol]
		v = 18015.28 * sT;
	}


	// 0.0 deg.C以上
	if ( T >= 273.15 )
	{
		// Ko
		So = 844.3 / T - 1.305;

		for ( i = 0 ; i <= 4 ; i++ )	xSo += x[i] * pow ( So, -i );

		Ko = Kmax_o * pow ( 10.0, xSo );

		// Kn
		Sn = 797.2 / T - 1.232;

		for ( i = 0 ; i <= 4 ; i++ )	xSn += x[i] * pow ( Sn, -i );

		Kn = Kmax_n * pow ( 10.0, xSn );

		// Ka
		Ka = ( Ko * Kn ) / ( 0.22 * Ko + 0.78 * Kn ) * 101325.16;

		// K
		K = 1.0 / Ka;
	}
	// 0.0 deg.C以下
	else	K = 0.0;

	// Bww, Cwww
	Bww = R * 1e6 * T * ( B2[0] + B2[1] * exp ( ( B2[2] / T ) ) );
	Cwww = pow ( R * 1e6 * T, 2 ) * pow ( C2[0] + C2[1] * exp ( C2[2] / T ), 2 ) * ( Bww * Bww );

	// Caaa
	for ( i = 0 ; i <= 2 ; i++ )	Caaa += C1[i] * pow ( T, -i );

	// Baa, Caww
	for ( i = 0 ; i <= 3 ; i++ )
	{
		Baa += B1[i] * pow ( T, -i );
		Caww += F[4] * exp ( ( F[i] * pow ( T, -i ) ) );
	}

	// Baw, Caaw
	for ( i = 0 ; i <= 4 ; i++ )
	{
		Baw += D[i] * pow ( T, -i );
		Caaw += E[i] * pow ( T, -i );
	}

	// f
	// 仮定
	f = 1.0;

	LoopCount = 0;
	while ( LoopCount < Loop )
	{
		// Xas
		Xas = ( P - f * Pvs ) / P;

		// f
		f1 = ( ( 1.0 + K * Pvs ) * ( P - Pvs ) - 0.5 * K * ( P * P - Pvs * Pvs ) ) * v / ( R0 * T ) + log ( 1.0 - Xas * P * K )
			+ Xas * Xas * P / ( R0 * T ) * Baa - 2.0 * Xas * Xas * P / ( R0 * T ) * Baw - ( P - Pvs - Xas * Xas * P ) / (R0 * T ) * Bww;

		f2 = pow ( Xas, 3 ) * P * P / pow ( R0 * T, 2 ) * Caaa + 3.0 * Xas * Xas * ( 1.0 - 2.0 * Xas ) * P * P / ( 2.0 * pow ( R0 * T, 2 ) ) * Caaw
			- 3.0 * Xas * Xas * ( 1.0 - Xas ) * P * P / pow ( R0 * T, 2 ) * Caww - Cwww * ( ( 1.0 + 2.0 * Xas ) * pow ( 1.0 - Xas, 2 ) * P * P - Pvs * Pvs )
			/ ( 2.0 * pow ( R0 * T, 2 ) );

		f3 = - Baa * Bww * Xas * Xas * ( 1.0 - 3.0 * Xas ) * ( 1.0 - Xas ) * P * P / pow ( R0 * T, 2 ) -2.0 * pow ( Xas, 3 ) * Baa * Baw * ( 2.0 - 3.0 * Xas )
			* P * P / pow ( R0 * T, 2 ) + 6.0 * Xas * Xas * Bww * Baw * pow ( 1.0 - Xas, 2 ) * P * P / pow ( R0 * T, 2 ) - 3.0 * pow ( Xas, 4 ) * P * P
			* Baa * Baa / ( 2.0 * pow ( R0 * T, 2 ) ) - 2.0 * Xas * Xas * Baw * Baw * ( 1.0 - Xas ) * ( 1.0 - 3.0 * Xas ) * P * P / pow ( R0 * T, 2 )
			- Bww * Bww * ( Pvs * Pvs - ( 1.0 + 3.0 * Xas ) * pow ( 1.0 - Xas, 3 ) * P * P ) / ( 2.0 * pow ( R0 * T, 2 ) );

		f_ = exp ( f1 + f2 + f3 );

		// 収束判定
		if ( fabs ( 1 - f_ / f ) < Error ) break;

		// 値の更新
		f = f_;

		LoopCount++;
	}
	return ( f_ );
}

/***********************************************
* 乾き空気熱伝導率関数		193 <= T <= 1727
* 入力：圧力[kPa]
*       温度[deg.C]
* 出力：熱伝導率[kW/(m*K)]
************************************************/
double Dair_Ptthc ( double P, double T )
{
	double thc;
	double thc0, thc1;
	double rho;
	double Tc, rhoc;
	double Tr, rhor;
	double l;
	double C1, C05;
	double C[5], D[5];

	int i;

	// 初期化
	thc1 = 0.0;

	C1 = 0.239503;
	C05 = 0.649768 * 1e-2;
	C[0] = 1.00000;
	C[1] = -0.192615 * 1e1;
	C[2] = 0.200383 * 1e1;
	C[3] = -0.107553 * 1e1;
	C[4] = 0.229414;

	D[0] = 0.402287;
	D[1] = 0.3566030;
	D[2] = -0.163159;
	D[3] = 0.138059;
	D[4] = -0.201725 * 1e-1;

	Tc = 132.5;
	rhoc = 314.3;
	l = 0.259778 * 1e-1;

	// 密度
	rho = Air_Ptxrho ( P, T, 0.0 );

	// Tr, rhor
	Tr = ( T + 273.15 ) / Tc;
	rhor = rho / rhoc;

	// thc0
	thc0 = C1 * Tr + C05 * pow ( Tr, 0.5 );
	for ( i = 0 ; i <= 4 ; i++ )	thc0 += C[i] * pow ( Tr, -i );

	// thc1
	for ( i = 0 ; i <= 4 ; i++ )	thc1 += D[i] * pow ( rhor, i + 1 );

	// 熱伝導率
	thc = l * ( thc0 + thc1 );

	return ( thc );
}

/***********************************************
* 乾き空気粘性係数関数		193 <= T <= 1727
* 入力：圧力[kPa]
*       温度[deg.C]
* 出力：粘性係数[Pa*s]
************************************************/
double Dair_Ptvisc ( double P, double T )
{
	double visc;
	double visc0, visc1;
	double rho;
	double Tc, rhoc;
	double Tr, rhor;
	double H;
	double A1, A05;
	double A[5], B[5];

	int i;

	// 初期化
	visc1 = 0.0;

	A1 = 0.128517;
	A05 = 0.260661 * 1e1;
	A[0] = -1.00000;
	A[1] = -0.709661;
	A[2] = 0.662534;
	A[3] = -0.197846;
	A[4] = 0.770147 * 1e-2;

	B[0] = 0.465601;
	B[1] = 0.126469 * 1e1;
	B[2] = -0.511425;
	B[3] = 0.274600;

	Tc = 132.5;
	rhoc = 314.3;
	H = 0.616090 * 1e-5;

	// 密度
	rho = Air_Ptxrho ( P, T, 0.0 );

	// Tr, rhor
	Tr = ( T + 273.15 ) / Tc;
	rhor = rho / rhoc;

	// visc0
	visc0 = A1 * Tr + A05 * pow ( Tr, 0.5 );
	for ( i = 0 ; i <= 4 ; i++ )	visc0 += A[i] * pow ( Tr, -i );

	// visc1
	for ( i = 0 ; i <= 3 ; i++ )	visc1 += B[i] * pow ( rhor, i + 1 );

	// 熱伝導率
	visc = H * ( visc0 + visc1 );

	return ( visc );
}
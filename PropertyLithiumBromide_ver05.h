#ifndef _PropertyLithiumBromide_ver05_h_
#define _PropertyLithiumBromide_ver05_h_

//////////////////////////////////////////////////////////////////////////////////
//																				//
//  LiBr 物性値関数(SI単位)														//
//																				//
//------------------------------------------------------------------------------//
//	関数名	   |説明				|単位		|入力パラメータ			|出展	//
//------------------------------------------------------------------------------//
//	sc_h_XT	   |溶液の比エンタルピ	|[kJ/kg]	|(濃度，温度)			|*1		//
//		       (使用範囲:40%～70%,15℃～165℃)									//
//	sc_T_Xh	   |溶液の温度			|[℃]		|(濃度，比エンタルピ)	|*1		//
//		       (使用範囲:40%～70%,15℃～165℃)									//
//	sc_T_XTsat |溶液の温度			|[℃]		|(濃度，飽和温度)		|*1		//
//		       (使用範囲:40%～70%,5℃～175℃)									//
//	sc_X_TTsat |溶液の濃度	　　＊陰解法|		|(温度,飽和温度)		|*1'	//
//																				//
//	sc_Tsat_XT |溶液の飽和温度	＊陰解法|		|(濃度, 温度)			|*1'	//
//																				//
//	sc_rho_XT  |溶液の密度			|[kg/m^3]	|(濃度，温度)			|*4		//
//																				//
//	sc_X_Trho  |溶液の濃度 ＊陰解法 |[kg/m^3]	|(温度，密度)			|*4'	//
//																				//
//	sc_visc_XT |溶液の粘性係数		|[Pa･s]		|(濃度，温度)			|*2		//
//																				//
//	sc_visc_XT2|溶液の粘性係数		|[Pa･s]		|(濃度，温度)			|*3		//
//																				//
// 	sc_thc_XT  |溶液の熱伝導率		|[W/m･K]	|(濃度，温度)			|*6		//
//																				//
// 	sc_thc_XT2 |溶液の熱伝導率		|[W/m･K]	|(濃度，温度)			|*3		//
//																				//
//	sc_cp_XT   |溶液の定圧比熱		|[kJ/kg･K]	|(濃度，温度)			|*2		//
//																				//
//	sc_cp_XT2  |溶液の定圧比熱		|[kJ/kg･K]	|(濃度，温度)			|*3		//
//																				//
//	sc_st_XT   |溶液の表面張力		|[N/m]		|(濃度，温度)			|*5-1	//
//																				//
//	sc_st_XT2  |溶液の表面張力		|[N/m]		|(濃度，温度)			|*3		//	おかしい
//																				//
//	sc_d_X     |溶液の物質拡散係数	|[m^2/s]	|(濃度)					|*3		//
//																				//
//	sc_d_XT    |溶液の物質拡散係数	|[m^2/s]	|(濃度，温度)			|*4		//
//																				//
//////////////////////////////////////////////////////////////////////////////////

// *1	2005 ASHRAE Handbook - Fundamentals (SI)
// *2	新版・第6版  冷凍空調便覧Ⅰ巻, (2010)
// *3	日本物性学会編, 新編 熱物性ハンドブック, 養賢堂, (2008)  出展記述なし
//		cp　枷場・植村, 冷凍35,397(1960) 815だと思われる． 
// *4	東京ガスの資料
// *5	Patterson, M.R. and Perez-Blanco, H., "Numerical fits of the properties of lithium-bromide water solutions", ASHRAE annual meeting, June 1988.
// *5-1	植村, 吸収式冷凍機用冷媒-吸収剤系の物性, 冷凍, Vol.52, No.600, 1975.
// *6	PROPERTIES OF LITHIUM BROMIDE-WATER SOLUTIONS AT HIGH TEMPERATURES AND CONDENSATIONS - PART I. Thermal Conductivity", ASHRAE Trans. 1990 Vol.96
//ASHRAE Trans. 1990 Vol.96 PROPERTIES OF LITHIUM BROMIDE-WATER SOLUTIONS AT HIGH TEMPERATURES AND CONDENSATIONS - PART II.  Density and Viscosity
//ASHRAE Trans. 1992 Vol.98 PROPERTIES OF LITHIUM BROMIDE-WATER SOLUTIONS AT HIGH TEMPERATURES AND CONDENSATIONS - PART III. Specific Heat


#include <cmath>

class PropertyLithiumBromide{
public:

	PropertyLithiumBromide();
	~PropertyLithiumBromide();

	double sc_h_XT    ( double _X , double _T    );
	double sc_T_Xh    ( double _X , double _h    );
	double sc_T_XTsat ( double _X , double _Tsat );
	double sc_X_TTsat ( double _T , double _Tsat );
	double sc_Tsat_XT ( double _X , double _T    );
	double sc_rho_XT  ( double _X , double _T    );
	double sc_X_Trho  ( double _T , double _rho  );
	double sc_visc_XT ( double _X , double _T    );
	double sc_visc_XT2( double _X , double _T    );
	double sc_thc_XT  ( double _X , double _T    );
	double sc_thc_XT2 ( double _X , double _T    );
	double sc_cp_XT   ( double _X , double _T    );
	double sc_cp_XT2  ( double _X , double _T    );
	double sc_st_XT   ( double _X , double _T    );
	double sc_st_XT2  ( double _X , double _T    );
	double sc_d_X	  ( double _X );
	double sc_d_XT    ( double _X , double _T    );

	double A[5];
	double B[5];
	double C[5];
	double D[4];
	double E[4];
	double F[2];
	double G[2];

	double An;
	double Bn;
	double Cn;
	double Dn;
	double En;

	double a;
	double b;

	double h;
	double T;
	double Tsat;
	double X;
	double rho;
	double thc;
	double visc;
	double cp;
	double st;
	double d;

	double t10;
	double a10;
	double b10;
	double c10;

	double dec;
	double eps;

	double E0;
	double E1;
	double E2;

};

#endif
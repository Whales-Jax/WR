/*
 *
 *  水の物性値 (1980SI 日本機械学会蒸気表)
 *
 *
 *                                     2003年7月28日　堀内雄次
 *
 *  ●使用範囲 0.01 ～ 350 ℃
 *
 *  ●引数・・・圧縮液,過熱蒸気：（圧力[kPa],温度[℃]）
 *              飽和液,飽和蒸気：（温度[℃]）
 *
 *  ●各物性関数
 *  ――――――┬――――┬――――┬――――┬――――┬――――┬――
 *      物性    │  単位  │ 圧縮液 │ 飽和液 │飽和蒸気│過熱蒸気│備考
 *  ――――――┼――――┼――――┼――――┼――――┼――――┼――
 *  比エントロピ│ kJ/kgK │sc_sl   │sat_sl  │sat_sv  │sh_sv   │
 *  比エクセルギ│ kJ/kg  │sc_el   │sat_el  │sat_ev  │sh_ev   │※1
 *  比エンタルピ│ kJ/kg  │sc_hl   │sat_hl  │sat_hv  │sh_hv   │
 *  比容積      │ m^3/kg │sc_vl   │sat_vl  │sat_vv  │sh_vv   │
 *  密度        │ kg/m^3 │sc_roul │sat_roul│sat_rouv│sh_rouv │※2
 *  定圧比熱    │ kJ/kgK │sc_cpl  │sat_cpl │sat_cpv │sh_cpv  │
 *  粘性係数    │ Pas    │sc_myul │sat_myul│sat_myuv│sh_myuv │
 *  熱伝導率    │ W/mK   │sc_laml │sat_laml│sat_lamv│sh_lamv │
 *  動粘性係数  │ m^2/s  │sc_nyul │sat_nyul│sat_nyuv│sh_nyuv │※2
 *  プラントル数│ -      │sc_Prl  │sat_Prl │sat_Prv │sh_Prv  │※2
 *  ――――――┴――――┴――――┴――――┴――――┴――――┴――
 *  ※1 エクセルギで使用する基準状態の圧力，温度を指定したい場合は基準圧力[kPa]，基準温度[℃]の順で引数に追加
 *      してください．デフォルトでは基準圧力101.325[kPa],基準温度20[℃]になっています．
 *      ex) sc_el(p,t)      //この場合，基準圧力101.325[kPa],基準温度20[℃]で計算
 *          sc_el(p,t,p0,t0)//この場合，基準圧力p0[kPa],基準温度t0[℃]で計算
 *  ※2 例えば密度ならば(密度)=1/(比容積)などのように他の関数の単純な組み合わせにより算出
 *  ※3 飽和液は圧縮液の第一引数に飽和圧を，飽和蒸気は過熱蒸気の第一引数に飽和圧を入れて算出している従属関数
 *
 *
 *  ●飽和温度と飽和圧力の関係
 *  ――――┬―――┬―――┬―――――
 *    物性  │ 単位 │関数名│  入力
 *  ――――┼―――┼―――┼―――――
 *  飽和圧力│ kPa  │ t_p  │温度[deg]
 *  飽和温度│ deg  │ p_t  │圧力[kPa]
 *  ――――┴―――┴―――┴―――――
 *  ※p_tはt_pにNewton法を適用した従属関数
 *
 *
 *  ●飽和温度と飽和圧力の関係-2
 *  ――――┬―――┬―――┬―――――
 *    物性  │ 単位 │関数名│  入力
 *  ――――┼―――┼―――┼―――――
 *  飽和圧力│ kPa  │ p_rg │温度[deg]
 *  飽和温度│ deg  │ t_rg │圧力[kPa]
 *  ――――┴―――┴―――┴―――――
 *  ※p_rgとt_rgはt_pとp_tと全く同様の関数を使用
 *
 *
 *  【補足説明】
 *
 *  1.関数引用元及び使用範囲
 *    記載の関数は「1980SI日本機械学会蒸気表」の一部を抜粋したもので，使用範囲は0.01～350℃である．全文を抜粋す
 *  ることで以下のことが可能になる．
 *  ・使用範囲を0.01～800℃に拡張
 *  ・ギブス関数,ヘルムホルツ関数,表面張力,静的誘電率,イオン積,pH,ラプラス係数の追加
 *
 *  2.wat_win.hとの差異
 *    property_water.hは以下のことを行うべくゼロから新規に作成されている．
 *  ・wat_win.hの間違い箇所の修正
 *    wat_win.hでは蒸気表に「df/dx=多項式」と記載されている場合，多項式を使用せず(f(x+delta_x)-f(x))/deltaにより
 *    df/dxを算出している．この方法では高温になるにつれ精度が落ちるため高温域では比容積などの決して負の値になら
 *    ないものが負になるなどの致命的な症状が起こる．これらの症状は過熱蒸気及び飽和蒸気の物性値全てに起こる．prop
 *    erty_water.hでは多項式を使用して算出している．
 *  ・エントロピを算出する関数の追加
 *
 *  【更新履歴】
 *  ・比エクセルギ関数の追加．コメントの訂正  2003-07-28 by 堀内雄次
 *  ・圧力の入力，出力を全てkPaに変換 2003-07-28 by 秋間雄二
 *  ・エントロピ，エンタルピの単位を全てkJ/kgに変換 2003-07-28 by 秋間雄二
 *  ・LiBr_SIを使用している者用にp_rg,t_rgを追加（各々t_pとp_tを呼出しているだけ）2003-07-28 by 秋間雄二
 *
 */

#ifndef __PROPERTY_WATER_H
#define __PROPERTY_WATER_H
#include <math.h>

/*━━━━━━━━━━━━━━━━━━━━━━━ 実用国際状態式(1967) ━━━━━━━━━━━━━━━━━━━━━━━*/

/*************************************
p109,116
臨界点，三重点及び誘導した定数の数値等
*************************************/
namespace
{
	double Tc = 647.3;      //K
	double Pc = 22120000.0; //J/m^3
	double vc = 0.00317;    //m^3/kg
	double I1 = 4.260321148;//-
	double T0 = 273.15;     //K
	double L0 = 1.574373327E1;
	double L1 =-3.417061978E1;
	double L2 = 1.931380707E1;
	double B0 = 1.683599274E+01;
	double alpha0 = 0.0;
	double alpha1 = 0.0;
}

/*************************************
p114
B関数の定数
*************************************/
#define USING_CONST_VALUE_NZLX                                                            \
int n[9],z[9][4],l[9],x[9][3];                                                            \
n[1] = 2; z[1][1] = 13; z[1][2] =  3;                                                     \
n[2] = 3; z[2][1] = 18; z[2][2] =  2; z[2][3] =  1;                                       \
n[3] = 2; z[3][1] = 18; z[3][2] = 10;                                                     \
n[4] = 2; z[4][1] = 25; z[4][2] = 14;                                                     \
n[5] = 3; z[5][1] = 32; z[5][2] = 28; z[5][3] = 24;                                       \
n[6] = 2; z[6][1] = 12; z[6][2] = 11;               l[6] = 1; x[6][1] = 14;               \
n[7] = 2; z[7][1] = 24; z[7][2] = 18;               l[7] = 1; x[7][1] = 19;               \
n[8] = 2; z[8][1] = 24; z[8][2] = 14;               l[8] = 2; x[8][1] = 54; x[8][2] = 27;

/*************************************
p115 7.1.1
部分領域1の定数
*************************************/
#define USING_CONST_VALUE_A                                                               \
double A[23],a[13];                                                                       \
A[0]  =  6.824687741E+03; A[12] = -2.616571843E-02; a[2]  = 5.362162162E-04;              \
A[1]  = -5.422063673E+02; A[13] =  1.522411790E-03; a[3]  = 1.720000000E+00;              \
A[2]  = -2.096666205E+04; A[14] =  2.284279054E-02; a[4]  = 7.342278489E-02;              \
A[3]  =  3.941286787E+04; A[15] =  2.421647003E+02; a[5]  = 4.975858870E-02;              \
A[4]  = -6.733277739E+04; A[16] =  1.269716088E-10; a[6]  = 6.537154300E-01;              \
A[5]  =  9.902381028E+04; A[17] =  2.074838328E-07; a[7]  = 1.150000000E-06;              \
A[6]  = -1.093911774E+05; A[18] =  2.174020350E-08; a[8]  = 1.510800000E-05;              \
A[7]  =  8.590841667E+04; A[19] =  1.105710498E-09; a[9]  = 1.418800000E-01;              \
A[8]  = -4.511168742E+04; A[20] =  1.293441934E+01; a[10] = 7.002753165E+00;              \
A[9]  =  1.418138926E+04; A[21] =  1.308119072E-05; a[11] = 2.995284926E-04;              \
A[10] = -2.017271113E+03; A[22] =  6.047626338E-14; a[12] = 2.040000000E-01;              \
A[11] =  7.982692717E+00; a[1]  =  8.438375405E-01;

/*************************************
p115 7.1.2
部分領域2の定数
*************************************/
#define USING_CONST_VALUE_B                                                               \
double bb,B[10][7],b[9][9];                                                               \
                            B[3][2] =  1.069036614E-01; B[9][0] =  1.936587558E+02;       \
B[0][1] =  2.856067796E+01; B[4][1] = -5.975336707E-01; B[9][1] = -1.388522425E+03;       \
B[0][2] = -5.438923329E+01; B[4][2] = -8.847535804E-02; B[9][2] =  4.126607219E+03;       \
B[0][3] =  4.330662834E-01; B[5][1] =  5.958051609E-01; B[9][3] = -6.508211677E+03;       \
B[0][4] = -6.547711697E-01; B[5][2] = -5.159303373E-01; B[9][4] =  5.745984054E+03;       \
B[0][5] =  8.565182058E-02; B[5][3] =  2.075021122E-01; B[9][5] = -2.693088365E+03;       \
B[1][1] =  6.670375918E-02; B[6][1] =  1.190610271E-01; B[9][6] =  5.235718623E+02;       \
B[1][2] =  1.388983801E-00; B[6][2] = -9.867174132E-02; bb      =  7.633333333E-01;       \
B[2][1] =  8.390104328E-02; B[7][1] =  1.683998803E-01; b[6][1] =  4.006073948E-01;       \
B[2][2] =  2.614670893E-02; B[7][2] = -5.809438001E-02; b[7][1] =  8.636081627E-02;       \
B[2][3] = -3.373439453E-02; B[8][1] =  6.552390126E-03; b[8][1] = -8.532322921E-01;       \
B[3][1] =  4.520918904E-01; B[8][2] =  5.710218649E-04; b[8][2] =  3.460208861E-01;

/*************************************
p116 7.1.5
飽和線の定数
*************************************/
#define USING_CONST_VALUE_K                                                               \
double k[10];                                                                             \
k[1] = -7.691234564E+00; k[4] =  6.423285504E+01; k[7] = 2.097506760E+01;                 \
k[2] = -2.608023696E+01; k[5] = -1.189646225E+02; k[8] = 1.000000000E+09;                 \
k[3] = -1.681706546E+02; k[6] =  4.167117320E+00; k[9] = 6.000000000E+00;

/**************************
p113 5 K関数
入力:換算温度[-]
出力:換算飽和圧力[-]
**************************/
double k_func(double t)
{
	USING_CONST_VALUE_K
	double x1,x2,x3,x4,p;

	x1=0.0;
	for(int i=1;i<=5;i++)
	{
		x1+=k[i]*pow(1.0-t,i);
	}
	x2=1.0+k[6]*(1.0-t)+k[7]*pow(1.0-t,2.0);
	x3=1.0-t;
	x4=k[8]*pow(1.0-t,2.0)+k[9];

	p=exp(1.0/t*x1/x2-x3/x4);

	return p;
}
/**************************
p116 8.2 L関数
入力:換算温度[-]
出力:部分領域2と3の境界の値
**************************/
double L_func(double t)
{
	double p;
	p=L0+L1*t+L2*t*t;
	return p;
}
double Ld_func(double t)
{
	double p;
	p=L1+2.0*L2*t;
	return p;
}
/****************************
p117 χ1
入力:換算圧力[-],換算温度[-]
出力:換算比容積[-]
****************************/
double chi1(double p,double t)
{
	USING_CONST_VALUE_A
	double Y,Z,x1,x2,x3,x4,x5,v;

	Y=1.0-a[1]*pow(t,2.0)-a[2]*pow(t,-6.0);
	Z=Y+sqrt(a[3]*pow(Y,2.0)-2.0*a[4]*t+2.0*a[5]*p);

	x1=A[11]*a[5]*pow(Z,-5.0/17.0);
	x2=A[12]+A[13]*t+A[14]*pow(t,2.0)+A[15]*pow(a[6]-t,10.0)+A[16]*pow(a[7]+pow(t,19.0),-1.0);
	x3=-pow(a[8]+pow(t,11.0),-1.0)*(A[17]+2.0*A[18]*p+3.0*A[19]*p*p);
	x4=-A[20]*pow(t,18.0)*(a[9]+pow(t,2.0))*(-3.0*pow(a[10]+p,-4.0)+a[11]);
	x5=3.0*A[21]*(a[12]-t)*p*p+4.0*A[22]*pow(t,-20.0)*pow(p,3.0);

	v=x1+x2+x3+x4+x5;

	return v;
}
/****************************
p117 σ1
入力:換算圧力[-],換算温度[-]
出力:換算エントロピ[-]
****************************/
double sigma1(double p,double t)
{
	USING_CONST_VALUE_A
	int nyu;
	double Y,Z,Yd,x1,x2,x3,x4,x5,x6,x7,x8,x9,s;

	Y=1.0-a[1]*pow(t,2.0)-a[2]*pow(t,-6.0);
	Z=Y+sqrt(a[3]*pow(Y,2.0)-2.0*a[4]*t+2.0*a[5]*p);
	Yd=-2.0*a[1]*t+6.0*a[2]*pow(t,-7.0);

	x1=alpha1;
	x2=A[0]*log(t);
	x3=0.0;
	for(nyu=2;nyu<=10;nyu++)
	{
		x3+=(nyu-1.0)*A[nyu]*pow(t,nyu-2.0);
	}
	x4=A[11]*((5.0/12.0*Z-(a[3]-1.0)*Y)*Yd+a[4])*pow(Z,-5.0/17.0);
	x5=(-A[13]-2.0*A[14]*t+10.0*A[15]*pow(a[6]-t,9.0)+19.0*A[16]*pow(a[7]+pow(t,19.0),-2.0)*pow(t,18.0))*p;
	x6=11.0*pow(a[8]+pow(t,11.0),-2.0)*pow(t,10.0)*(A[17]*p+A[18]*p*p+A[19]*p*p*p);
	x7=A[20]*pow(t,17.0)*(18.0*a[9]+20.0*t*t)*(pow(a[10]+p,-3.0)+a[11]*p);
	x8=A[21]*p*p*p;
	x9=20.0*A[22]*pow(t,-21.0)*pow(p,4.0);

	s=-x1+x2-x3+x4+x5-x6+x7+x8+x9;

	return s;
}
/****************************
p117 ε1
入力:換算圧力[-],換算温度[-]
出力:換算エンタルピ[-]
****************************/
double epsilon1(double p,double t)
{
	USING_CONST_VALUE_A
	int nyu;
	double Y,Z,Yd,x1,x2,x3,x4,x5,x6,x7,x8,h;

	Y =1.0-a[1]*pow(t,2.0)-a[2]*pow(t,-6.0);
	Z =Y+sqrt(a[3]*pow(Y,2.0)-2.0*a[4]*t+2.0*a[5]*p);
	Yd=-2.0*a[1]*t+6.0*a[2]*pow(t,-7.0);

	x1=alpha0;
	x2=A[0]*t;
	x3=0.0;
	for(nyu=1;nyu<=10;nyu++)
	{
		x3+=(nyu-2.0)*A[nyu]*pow(t,nyu-1.0);
	}
	x4=A[11]*(Z*(17.0*(Z/29.0-Y/12.0)+5.0*t*Yd/12.0)+a[4]*t-(a[3]-1.0)*t*Y*Yd)*pow(Z,-5.0/17.0);
	x5=(A[12]-A[14]*pow(t,2.0)+A[15]*(9.0*t+a[6])*pow(a[6]-t,9.0)+A[16]*(20.0*pow(t,19.0)+a[7])*pow(a[7]+pow(t,19.0),-2.0))*p;
	x6=(12.0*pow(t,11.0)+a[8])*pow(a[8]+pow(t,11.0),-2.0)*(A[17]*p+A[18]*p*p+A[19]*p*p*p);
	x7=A[20]*pow(t,18.0)*(17.0*a[9]+19.0*t*t)*(pow(a[10]+p,-3.0)+a[11]*p);
	x8=A[21]*a[12]*p*p*p+21.0*A[22]*pow(t,-20.0)*pow(p,4.0);

	h=x1+x2-x3+x4+x5-x6+x7+x8;

	return(h);
}
/****************************
p117 χ2
入力:換算圧力[-],換算温度[-]
出力:換算比容積[-]
****************************/
double chi2(double p,double t)
{
	USING_CONST_VALUE_B
	USING_CONST_VALUE_NZLX
	int myuu,nyu,lambda;
	double v,y1,y2,y3,y4,sigma,sigma_bunbo,sigma_bunshi;
	double X=exp(bb*(1.0-t));

	y1=I1*t/p;
	y2=0.0;
	for(myuu=1;myuu<=5;myuu++)
	{
		sigma=0.0;
		for(nyu=1;nyu<=n[myuu];nyu++)
		{
			sigma+=B[myuu][nyu]*pow(X,z[myuu][nyu]);
		}
		y2+=myuu*pow(p,myuu-1.0)*sigma;
	}
	y3=0.0;
	for(myuu=6;myuu<=8;myuu++)
	{
		sigma_bunshi=sigma_bunbo=0.0;
		for(nyu=1;nyu<=n[myuu];nyu++)
		{
			sigma_bunshi+=B[myuu][nyu]*pow(X,z[myuu][nyu]);
		}
		for(lambda=1;lambda<=l[myuu];lambda++)
		{
			sigma_bunbo+=b[myuu][lambda]*pow(X,x[myuu][lambda]);
		}
		y3+=(myuu-2.0)*pow(p,1.0-myuu)*sigma_bunshi/pow((pow(p,2.0-myuu)+sigma_bunbo),2.0);
	}
	sigma=0.0;
	for(nyu=0;nyu<=6;nyu++)
	{
		sigma+=B[9][nyu]*pow(X,nyu);
	}
	y4=11.0*pow(p/L_func(t),10.0)*sigma;

	v=y1-y2-y3+y4;

	return v;
}
/****************************
p118 σ2
入力:換算圧力[-],換算温度[-]
出力:換算エントロピ[-]
****************************/
double sigma2(double p,double t)
{
	USING_CONST_VALUE_B
	USING_CONST_VALUE_NZLX
	int myuu,nyu,lambda;
	double s,y1,y2,y3,y4,y5,y6,y7,sigma,sigma_bunbo,sigma_bunshi,bunshi;
	double X=exp(bb*(1.0-t));

	y1=alpha1;
	y2=I1*log(p);
	y3=B0*log(t);
	y4=0.0;
	for(nyu=1;nyu<=5;nyu++)
	{
		y4+=(nyu-1.0)*B[0][nyu]*pow(t,nyu-2.0);
	}
	y5=0.0;
	for(myuu=1;myuu<=5;myuu++)
	{
		sigma=0.0;
		for(nyu=1;nyu<=n[myuu];nyu++)
		{
			sigma+=z[myuu][nyu]*B[myuu][nyu]*pow(X,z[myuu][nyu]);
		}
		y5+=pow(p,myuu)*sigma;
	}
	y5*=bb;
	y6=0.0;
	for(myuu=6;myuu<=8;myuu++)
	{
		bunshi=0.0;
		for(nyu=1;nyu<=n[myuu];nyu++)
		{
			sigma_bunshi=sigma_bunbo=0.0;
			for(lambda=1;lambda<=l[myuu];lambda++)
			{
				sigma_bunshi+=x[myuu][lambda]*b[myuu][lambda]*pow(X,x[myuu][lambda]);
			}
			for(lambda=1;lambda<=l[myuu];lambda++)
			{
				sigma_bunbo+=b[myuu][lambda]*pow(X,x[myuu][lambda]);
			}
			sigma=z[myuu][nyu]-sigma_bunshi/(pow(p,2.0-myuu)+sigma_bunbo);
			bunshi+=B[myuu][nyu]*pow(X,z[myuu][nyu])*sigma;
		}
		sigma_bunbo=0.0;
		for(lambda=1;lambda<=l[myuu];lambda++)
		{
			sigma_bunbo+=b[myuu][lambda]*pow(X,x[myuu][lambda]);
		}
		y6+=bunshi/(pow(p,2.0-myuu)+sigma_bunbo);
	}
	y6*=bb;
	y7=0.0;
	for(nyu=0;nyu<=6;nyu++)
	{
		y7+=(10.0*Ld_func(t)/L_func(t)+nyu*bb)*B[9][nyu]*pow(X,nyu);
	}
	y7*=p*pow(p/L_func(t),10.0);

	s=-y1-y2+y3-y4-y5-y6+y7;

	return s;
}
/****************************
p117 ε2
入力:換算圧力[-],換算温度[-]
出力:換算エンタルピ[-]
****************************/
double epsilon2(double p,double t)
{
	USING_CONST_VALUE_B
	USING_CONST_VALUE_NZLX
	int myuu,nyu,lambda;
	double h,y1,y2,y3,y4,y5,y6,sigma,sigma_bunbo,sigma_bunshi,bunshi;
	double X=exp(bb*(1.0-t));

	y1=alpha0;
	y2=B0*t;
	y3=0.0;
	for(nyu=1;nyu<=5;nyu++)
	{
		y3+=B[0][nyu]*(nyu-2.0)*pow(t,nyu-1.0);
	}
	y4=0.0;
	for(myuu=1;myuu<=5;myuu++)
	{
		sigma=0.0;
		for(nyu=1;nyu<=n[myuu];nyu++)
		{
			sigma+=B[myuu][nyu]*(1+z[myuu][nyu]*bb*t)*pow(X,z[myuu][nyu]);
		}
		y4+=pow(p,myuu)*sigma;
	}
	y5=0.0;
	for(myuu=6;myuu<=8;myuu++)
	{
		bunshi=0.0;
		for(nyu=1;nyu<=n[myuu];nyu++)
		{
			sigma_bunshi=sigma_bunbo=0.0;
			for(lambda=1;lambda<=l[myuu];lambda++)
			{
				sigma_bunshi+=x[myuu][lambda]*b[myuu][lambda]*pow(X,x[myuu][lambda]);
			}
			for(lambda=1;lambda<=l[myuu];lambda++)
			{
				sigma_bunbo+=b[myuu][lambda]*pow(X,x[myuu][lambda]);
			}
			sigma=1.0+z[myuu][nyu]*bb*t-bb*t*sigma_bunshi/(pow(p,2.0-myuu)+sigma_bunbo);
			bunshi+=B[myuu][nyu]*pow(X,z[myuu][nyu])*sigma;
		}
		sigma_bunbo=0.0;
		for(lambda=1;lambda<=l[myuu];lambda++)
		{
			sigma_bunbo+=b[myuu][lambda]*pow(X,x[myuu][lambda]);
		}
		y5+=bunshi/(pow(p,2.0-myuu)+sigma_bunbo);
	}
	y6=0.0;
	for(nyu=0;nyu<=6;nyu++)
	{
		y6+=(1.0+t*(10.0*Ld_func(t)/L_func(t)+nyu*bb))*B[9][nyu]*pow(X,nyu);
	}
	y6*=p*pow(p/L_func(t),10.0);

	h=y1+y2-y3-y4-y5+y6;

	return h;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━ 定圧比熱 cp の式 ━━━━━━━━━━━━━━━━━━━━━━━━*/

/****************************
p120 φ1
入力:換算圧力[-],換算温度[-]
出力:換算定圧比熱[-]
****************************/
double phi1(double p,double t)
{
	USING_CONST_VALUE_A
	int nyu;
	double Y,Y1,Y2,Z0,Z,Z1,Z2,x1,x2,x3,x4,x5,x6,x7,phi;

	Y =1.0-a[1]*pow(t,2.0)-a[2]*pow(t,-6.0);
	Y1=-2.0*a[1]*t+6.0*a[2]*pow(t,-7.0);
	Y2=-2.0*a[1]-42.0*a[2]*pow(t,-8.0);
	Z0=a[3]*Y*Y-2.0*a[4]*t+2.0*a[5]*p;
	Z =Y+sqrt(Z0);
	Z1=Y1+(a[3]*Y*Y1-a[4])*pow(Z0,-0.5);
	Z2=Y2+a[3]*pow(Z0,-0.5)*(Y1*Y1+Y*Y2)-pow(a[3]*Y*Y1-a[4],2.0)*pow(Z0,-1.5);

	x1=A[0]/t;
	x2=0.0;
	for(nyu=0;nyu<=7;nyu++)
	{
		x2+=(nyu+1.0)*(nyu+2.0)*A[nyu+3]*pow(t,nyu);
	}
	x3=A[11]*(((12.0/29.0)*Z-Y)*(pow(Z,-5.0/17.0)*Z2-5.0/17.0*pow(Z,-22.0/17.0)*pow(Z1,2.0))+(24.0/29.0*Z1-2.0*Y1)*pow(Z,-5.0/17.0)*Z1+(17.0/29.0*Z2-17.0/12.0*Y2)*pow(Z,12.0/17.0));
	x4=p*(2.0*A[14]+90.0*A[15]*pow(a[6]-t,8.0)+722.0*A[16]*pow(t,36.0)*pow(a[7]+pow(t,19.0),-3.0)-342.0*A[16]*pow(t,17.0)*pow(a[7]+pow(t,19.0),-2.0));
	x5=(242.0*pow(t,20.0)*pow(a[8]+pow(t,11.0),-3.0)-110.0*pow(t,9.0)*pow(a[8]+pow(t,11.0),-2.0))*(A[17]*p+A[18]*p*p+A[19]*p*p*p);
	x6=A[20]*pow(t,16.0)*(306.0*a[9]+380.0*t*t)*(pow(a[10]+p,-3.0)+a[11]*p);
	x7=420.0*A[22]*pow(t,-22.0)*pow(p,4.0);

	phi=-t*(-x1+x2+x3+x4-x5-x6+x7);

	return phi;
}
/****************************
p120 φ2
入力:換算圧力[-],換算温度[-]
出力:換算定圧比熱[-]
****************************/
double phi2(double p,double t)
{
	USING_CONST_VALUE_B
	int nyu;
	double phi,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,sigma1,sigma2,sigma3;
	double X=exp(bb*(1.0-t));
	double F[9][3],G[9][3];

	#define _X0(i) (pow(X,i))
	#define _X1(i) (-i*bb*pow(X,i))
	#define _X2(i) (i*i*bb*bb*pow(X,i))
	F[6][0]=(B[6][1]*_X0(12.0)+B[6][2]*_X0(11.0))*pow(p,4.0);
	F[7][0]=(B[7][1]*_X0(24.0)+B[7][2]*_X0(18.0))*pow(p,5.0);
	F[8][0]=(B[8][1]*_X0(24.0)+B[8][2]*_X0(14.0))*pow(p,6.0);
	F[6][1]=(B[6][1]*_X1(12.0)+B[6][2]*_X1(11.0))*pow(p,4.0);
	F[7][1]=(B[7][1]*_X1(24.0)+B[7][2]*_X1(18.0))*pow(p,5.0);
	F[8][1]=(B[8][1]*_X1(24.0)+B[8][2]*_X1(14.0))*pow(p,6.0);
	F[6][2]=(B[6][1]*_X2(12.0)+B[6][2]*_X2(11.0))*pow(p,4.0);
	F[7][2]=(B[7][1]*_X2(24.0)+B[7][2]*_X2(18.0))*pow(p,5.0);
	F[8][2]=(B[8][1]*_X2(24.0)+B[8][2]*_X2(14.0))*pow(p,6.0);
	G[6][0]=1.0+b[6][1]*_X0(14.0)*pow(p,4.0);
	G[7][0]=1.0+b[7][1]*_X0(19.0)*pow(p,5.0);
	G[8][0]=1.0+(b[8][1]*_X0(54.0)+b[8][2]*_X0(27.0))*pow(p,6.0);
	G[6][1]=b[6][1]*_X1(14.0)*pow(p,4.0);
	G[7][1]=b[7][1]*_X1(19.0)*pow(p,5.0);
	G[8][1]=(b[8][1]*_X1(54.0)+b[8][2]*_X1(27.0))*pow(p,6.0);
	G[6][2]=b[6][1]*_X2(14.0)*pow(p,4.0);
	G[7][2]=b[7][1]*_X2(19.0)*pow(p,5.0);
	G[8][2]=(b[8][1]*_X2(54.0)+b[8][2]*_X2(27.0))*pow(p,6.0);
	#undef _X0
	#undef _X1
	#undef _X2

	#define B_X(h,i,j,k) ((h)*B[i][j]*pow(X,(k)))
	y1=B0/t;
	y2=2.0*B[0][3];
	y3=6.0*B[0][4]*t;
	y4=12.0*B[0][5]*t*t;
	y5=(B_X(169.0,1,1,13.0)+B_X(9.0,1,2,3.0))*bb*bb*p;
	y6=(B_X(324.0,2,1,18.0)+B_X(4.0,2,2,2.0)+B_X(1.0,2,3,1.0))*bb*bb*p*p;
	y7=(B_X(324.0,3,1,18.0)+B_X(100.0,3,2,10.0))*bb*bb*p*p*p;
	y8=(B_X(625.0,4,1,25.0)+B_X(196.0,4,2,14.0))*bb*bb*p*p*p*p;
	y9=(B_X(1024.0,5,1,32.0)+B_X(784.0,5,2,28.0)+B_X(576.0,5,3,24.0))*bb*bb*p*p*p*p*p;
	y10=0.0;
	for(nyu=6;nyu<=8;nyu++)
	{
		y10+=F[nyu][2]/G[nyu][0]-(2.0*F[nyu][1]*G[nyu][1]+F[nyu][0]*G[nyu][2])/pow(G[nyu][0],2.0)+2.0*F[nyu][0]*pow(G[nyu][1],2.0)/pow(G[nyu][0],3.0);
	}
	sigma1=sigma2=sigma3=0.0;
	for(nyu=0;nyu<=6;nyu++)
	{
		sigma1+=B_X(1.0,9,nyu,nyu);
		sigma2+=B_X(nyu,9,nyu,nyu);
		sigma3+=B_X(nyu*nyu,9,nyu,nyu);
	}
	y11=pow(p/L_func(t),11.0)*((110.0/L_func(t)*pow(Ld_func(t),2.0)-10.0*2.0*L2)*sigma1+20.0*bb*Ld_func(t)*sigma2+bb*bb*L_func(t)*sigma3);
	#undef B_X

	phi=-t*(-y1+y2+y3+y4-y5-y6-y7-y8-y9-y10+y11);

	return phi;
}

/*━━━━━━━━━━━━━━━━━━━━━━ 粘性係数の国際補間式(1975) ━━━━━━━━━━━━━━━━━━━━━*/

/*******************************
p121 μ
入力:温度[K],密度[kg/m3]
出力:粘性係数[μPas]
*******************************/
double troumyuu(double t,double rou)
{
	int i,j,k;
	double Tast,rouast,a[4],b[6][5],sigma,myuu0,myuu;
	Tast   = 647.27;  a[0] = 0.0181583; a[2] = 0.0105287;
	rouast = 317.763; a[1] = 0.0177624; a[3] =-0.0036744;
	b[0][0]= 0.501938; b[1][0]= 0.162888; b[2][0]=-0.130356; b[3][0]= 0.907919; b[4][0]=-0.551119; b[5][0]= 0.146543;
	b[0][1]= 0.235622; b[1][1]= 0.789393; b[2][1]= 0.673665; b[3][1]= 1.207552; b[4][1]= 0.0670665;b[5][1]=-0.0843370;
	b[0][2]=-0.274637; b[1][2]=-0.743539; b[2][2]=-0.959456; b[3][2]=-0.687343; b[4][2]=-0.497089; b[5][2]= 0.195286;
	b[0][3]= 0.145831; b[1][3]= 0.263129; b[2][3]= 0.347247; b[3][3]= 0.213486; b[4][3]= 0.100754; b[5][3]=-0.032932;
	b[0][4]=-0.0270448;b[1][4]=-0.0253093;b[2][4]=-0.0267758;b[3][4]=-0.0822904;b[4][4]= 0.0602253;b[5][4]=-0.0202595;

	sigma=0.0;
	for(k=0;k<=3;k++)
	{
		sigma+=a[k]*pow(Tast/t,k);
	}
	myuu0=sqrt(t/Tast)/sigma;
	sigma=0.0;
	for(i=0;i<=5;i++)
	{
		for(j=0;j<=4;j++)
		{
			sigma+=b[i][j]*pow(Tast/t-1.0,i)*pow(rou/rouast-1.0,j);
		}
	}

	myuu=myuu0*exp(rou/rouast*sigma);

	return(myuu);
}

/*━━━━━━━━━━━━━━━━━━━━━ 熱伝導率の実用国際補間式(1977) ━━━━━━━━━━━━━━━━━━━━*/

/*******************************
p122 λ
入力:温度[K],密度[kg/m3]
出力:熱伝達率[W/mK]
*******************************/
double troulambda(double t,double rou)
{
	double lambda0,lambda_,lambdad,lambda;
	double Tast,rouast,B1,B2,C1,C2,C3,C4,C5,C6,a0,a1,a2,a3,b0,b1,b2,d1,d2,d3,d4;
	Tast   = 647.30;       C3 = -6.17937;     a2 =  1.56146E-02; d2 =  1.18520E-02;
	rouast = 317.7;        C4 =  3.08976E-03; a3 = -4.22464E-03; d3 =  1.69937E-03;
	B1     = -1.71587E-01; C5 =  8.22994E-02; b0 = -3.97070E-01; d4 = -1.02000;
	B2     =  2.39219;     C6 =  1.00932E+01; b1 =  4.00302E-01;
	C1     =  6.42857E-01; a0 =  1.02811E-02; b2 =  1.06000;
	C2     = -4.11717;     a1 =  2.99621E-02; d1 =  7.01309E-02;

	t/=Tast;
	rou/=rouast;

	double Q,R,S,dT;
	dT=fabs(t-1.0)+C4;
	Q=2.0+C5*pow(dT,-0.6);
	R=Q+1.0;
	S=((t>=1.0)?(1.0/dT):(C6*pow(dT,-0.6)));

	lambda0=sqrt(t)*(a0+a1*t+a2*t*t+a3*t*t*t);
	lambda_=b0+b1*rou+b2*exp(B1*pow(rou+B2,2.0));
	lambdad=(d1*pow(t,-10.0)+d2)*pow(rou,1.8)*exp(C1*(1.0-pow(rou,2.8)))+d3*S*pow(rou,Q)*exp(Q/R*(1.0-pow(rou,R)))+d4*exp(C2*pow(t,1.5)+C3*pow(rou,-5.0));

	lambda=lambda0+lambda_+lambdad;

	return lambda;
}

/*━━━━━━━━━━━━━━━━━━ 上記関数の引数及び戻り値をSIに変換した式 ━━━━━━━━━━━━━━━━━━*/

/*********************************
飽和温度[℃]と飽和圧力[kPa]の関係式
*********************************/
//飽和圧力[kPa]
inline double t_p(double t)
{
	double p;
	t=(t+T0)/Tc;
	p=k_func(t)*Pc;
	return p/1000.0;
}
//飽和温度[℃]
double p_t(double p)
{
	int i;double x[2],y[2],dec=1.0+1e-7,eps=1e-7;
	double t=pow(p,0.2)*14.13-44.0;
	do{
		for(i=0;i<=1;i++)
		{
			x[i] = t*((i==0)?(dec):(1.0));
			y[i] = t_p(x[i])-p;
		}
		t = x[0]-y[0]*(x[1]-x[0])/(y[1]-y[0]);
	}while(fabs(y[1])>=eps);
	return(t);
}

//旧関数p_rg,t_rgの定義(各々t_pとp_tと同様)
inline double p_rg(double t){ return( t_p(t) ); }
inline double t_rg(double p){ return( p_t(p) ); }

/*********************************
比容積[m^3/kg]
*********************************/
//圧縮液
inline double sc_vl(double p,double t)
{
	p*=1000;
	double v;
	p=p/Pc;
	t=(t+T0)/Tc;
	v=chi1(p,t)*vc;
	return v;
}
//過熱蒸気
inline double sh_vv(double p,double t)
{
	p*=1000;
	double v;
	p=p/Pc;
	t=(t+T0)/Tc;
	v=chi2(p,t)*vc;
	return v;
}
//飽和液
inline double sat_vl(double t)
{
	double p,v;
	p=t_p(t);
	v=sc_vl(p,t);
	return v;
}
//飽和蒸気
inline double sat_vv(double t)
{
	double p,v;
	p=t_p(t);
	v=sh_vv(p,t);
	return v;
}
/*********************************
密度[kg/m^3]
*********************************/
//圧縮液
inline double sc_roul(double p,double t)
{
	return 1.0/sc_vl(p,t);
}
//過熱蒸気
inline double sh_rouv(double p,double t)
{
	return 1.0/sh_vv(p,t);
}
//飽和液
inline double sat_roul(double t)
{
	return 1.0/sat_vl(t);
}
//飽和蒸気
inline double sat_rouv(double t)
{
	return 1.0/sat_vv(t);
}
/*********************************
エントロピ[kJ/kgK]
*********************************/
//圧縮液
inline double sc_sl(double p,double t)
{
	p*=1000;
	double s;
	p=p/Pc;
	t=(t+T0)/Tc;
	s=sigma1(p,t)*(Pc*vc/Tc);
	return s/1000;
}
//過熱蒸気
inline double sh_sv(double p,double t)
{
	p*=1000;
	double s;
	p=p/Pc;
	t=(t+T0)/Tc;
	s=sigma2(p,t)*(Pc*vc/Tc);
	return s/1000;
}
//飽和液
inline double sat_sl(double t)
{
	double p,s;
	p=t_p(t);
	s=sc_sl(p,t);
	return s;
}
//飽和蒸気
inline double sat_sv(double t)
{
	double p,s;
	p=t_p(t);
	s=sh_sv(p,t);
	return s;
}
/*********************************
エンタルピ[kJ/kg]
*********************************/
//圧縮液
inline double sc_hl(double p,double t)
{
	p*=1000;
	double h;
	p=p/Pc;
	t=(t+T0)/Tc;
	h=epsilon1(p,t)*(Pc*vc);
	return h/1000;
}
//過熱蒸気
inline double sh_hv(double p,double t)
{
	p*=1000;
	double h;
	p=p/Pc;
	t=(t+T0)/Tc;
	h=epsilon2(p,t)*(Pc*vc);
	return h/1000;
}
//飽和液
inline double sat_hl(double t)
{
	double p,h;
	p=t_p(t);
	h=sc_hl(p,t);
	return h;
}
//飽和蒸気
inline double sat_hv(double t)
{
	double p,h;
	p=t_p(t);
	h=sh_hv(p,t);
	return h;
}
/*********************************
定圧比熱[kJ/kgK]
*********************************/
//圧縮液
inline double sc_cpl(double p,double t)
{
	p*=1000;
	double cp;
	p=p/Pc;
	t=(t+T0)/Tc;
	cp=phi1(p,t)*(Pc*vc/Tc);
	return cp/1000;
}
//過熱蒸気
inline double sh_cpv(double p,double t)
{
	p*=1000;
	double cp;
	p=p/Pc;
	t=(t+T0)/Tc;
	cp=phi2(p,t)*(Pc*vc/Tc);
	return cp/1000;
}
//飽和液
inline double sat_cpl(double t)
{
	double p,cp;
	p=t_p(t);
	cp=sc_cpl(p,t);
	return cp;
}
//飽和蒸気
inline double sat_cpv(double t)
{
	double p,cp;
	p=t_p(t);
	cp=sh_cpv(p,t);
	return cp;
}
/*********************************
粘性係数[Pas]
*********************************/
//圧縮液
inline double sc_myul(double p,double t)
{
	double myuu,rou;
	rou=sc_roul(p,t);
	t=t+T0;
	myuu=troumyuu(t,rou)*1e-6;
	return myuu;
}
//過熱蒸気
inline double sh_myuv(double p,double t)
{
	double myuu,rou;
	rou=sh_rouv(p,t);
	t=t+T0;
	myuu=troumyuu(t,rou)*1e-6;
	return myuu;
}
//飽和液
inline double sat_myul(double t)
{
	double p,myuu;
	p=t_p(t);
	myuu=sc_myul(p,t);
	return myuu;
}
//飽和蒸気
inline double sat_myuv(double t)
{
	double p,myuu;
	p=t_p(t);
	myuu=sh_myuv(p,t);
	return myuu;
}
/*********************************
熱伝導率[W/mK]
*********************************/
//圧縮液
inline double sc_laml(double p,double t)
{
	double lambda,rou;
	rou=sc_roul(p,t);
	t=t+T0;
	lambda=troulambda(t,rou);
	return lambda;
}
//過熱蒸気
inline double sh_lamv(double p,double t)
{
	double lambda,rou;
	rou=sh_rouv(p,t);
	t=t+T0;
	lambda=troulambda(t,rou);
	return lambda;
}
//飽和液
inline double sat_laml(double t)
{
	double p,lambda;
	p=t_p(t);
	lambda=sc_laml(p,t);
	return lambda;
}
//飽和蒸気
inline double sat_lamv(double t)
{
	double p,lambda;
	p=t_p(t);
	lambda=sh_lamv(p,t);
	return lambda;
}
/*********************************
動粘性係数[m^2/s]
*********************************/
//圧縮液
inline double sc_nyul(double p,double t)
{
	double nyu;
	nyu=sc_myul(p,t)*sc_vl(p,t);
	return nyu;
}
//過熱蒸気
inline double sh_nyuv(double p,double t)
{
	double nyu;
	nyu=sh_myuv(p,t)*sh_vv(p,t);
	return nyu;
}
//飽和液
inline double sat_nyul(double t)
{
	double nyu;
	nyu=sat_myul(t)*sat_vl(t);
	return nyu;
}
//飽和蒸気
inline double sat_nyuv(double t)
{
	double nyu;
	nyu=sat_myuv(t)*sat_vv(t);
	return nyu;
}
/*********************************
プラントル数[-]
*********************************/
//圧縮液
inline double sc_Prl(double p,double t)
{
	double pr;
	pr=sc_nyul(p,t)*sc_cpl(p,t)*1000/sc_vl(p,t)/sc_laml(p,t);
	return pr;
}
//過熱蒸気
inline double sh_Prv(double p,double t)
{
	double pr;
	pr=sh_nyuv(p,t)*sh_cpv(p,t)*1000/sh_vv(p,t)/sh_lamv(p,t);
	return pr;
}
//飽和液
inline double sat_Prl(double t)
{
	double pr;
	pr=sat_nyul(t)*sat_cpl(t)*1000/sat_vl(t)/sat_laml(t);
	return pr;
}
//飽和蒸気
inline double sat_Prv(double t)
{
	double pr;
	pr=sat_nyuv(t)*sat_cpv(t)*1000/sat_vv(t)/sat_lamv(t);
	return pr;
}
/*********************************
比エクセルギ[kJ/kg]
*********************************/
//!!基準は常に圧縮液であることに注意!!
//圧縮液
inline double sc_el(double p,double t,double p0=101.325,double t0=20.0)
{
	p*=1000;
	double e;
	e=(sc_hl(p,t)-sc_hl(p0,t0))-(t0+T0)*(sc_sl(p,t)-sc_sl(p0,t0));
	return e/1000;
}
//過熱蒸気
inline double sh_ev(double p,double t,double p0=101.325,double t0=20.0)
{
	p*=1000;
	double e;
	e=(sh_hv(p,t)-sc_hl(p0,t0))-(t0+T0)*(sh_sv(p,t)-sc_sl(p0,t0));
	return e/1000;
}
//飽和液
inline double sat_el(double t,double p0=101.325,double t0=20.0)
{
	double e;
	e=(sat_hl(t)-sc_hl(p0,t0))-(t0+T0)*(sat_hl(t)-sc_sl(p0,t0));
	return e/1000;
}
//飽和蒸気
inline double sat_ev(double t,double p0=101.325,double t0=20.0)
{
	double e;
	e=(sat_hv(t)-sc_hl(p0,t0))-(t0+T0)*(sat_hv(t)-sc_sl(p0,t0));
	return e/1000;
}
#endif//__PROPERTY_WATER_H

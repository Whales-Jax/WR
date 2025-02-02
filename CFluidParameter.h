#ifndef __CFluidParameter_h
#define __CFluidParameter_h

#include <iostream>
#include <fstream>
#include "DXProperty_ver06.h"

/*!
流体状態変数
*/
class CFluidParameter{
public:

	CFluidParameter();
	~CFluidParameter();

	int CycleNUM;

	double G;//!<質量流量[kg/s]
	double P;//!<圧力[kPa]
	double h;//!<比エンタルピ[kJ/kg]
	double x;//!<乾き度
	double X;//!<濃度
	double s;//!<比エントロピ[]
	double D;//!<密度[kg/m3]
	double u;//!<比内部エネルギー[kJ/kg]
	double M;//!<質量[kg]
	double Q;//!<熱量[kW]
	double q;//!<熱流束[kW/m2]
	double T;//!<温度[degC]
	double Tw;//!<湿球温度[degC]
	double dT;//!<温度差[degC]
	double v;//!<比体積[m3/kg]
	double V;//!<体積流量[m3/s]
	double V_lpm;//!<体積流量[l/min]
	double velocity;//!<流速[m/s]
	double cp;//!<定圧比熱[kJ/(kg*K)]
	double thc;//!<熱伝導率[W/(m*K)]
	double visc;//!<粘性係数[uPa*s]
	double st;//!<表面張力[N/m]
	double cv;
	double satT;//!飽和温度
	double phi;//!相対湿度
	double nu;//!動粘度
	
	double checkT;//!<温度[degC]

	double cpL;//!<定圧比熱[kJ/(kg*K)]
	double thcL;//!<熱伝導率[W/(m*K)]
	double viscL;//!<粘性係数[uPa*s]
	double cpV;//!<定圧比熱[kJ/(kg*K)]
	double thcV;//!<熱伝導率[W/(m*K)]
	double viscV;//!<粘性係数[uPa*s]


	double alpha;//!<熱伝達率
	double alpharate;//!<熱伝達率
	double alphaL;//!<熱伝達率
	double alphaV;//!<熱伝達率
	double alphaLO;//!<熱伝達率
	double alphaAS;//!<熱伝達率
	double alphaA;//!<熱伝達率
	double alphaDS;//!<熱伝達率
	double alphaD;//!<熱伝達率
	double pd;//!<圧力損失[kPa]
	double pd_H;//!<ヘッド分の圧力損失[kPa]
	double dp;//!<圧力損失[kPa/m]

	double Gi;//!<入口流量[kg/s]
	double Pi;//!<入口圧力[kPs]
	double hi;//!<入口比エンタルピ[kJ/kg]
	double Xi;//!<入口濃度
	double Di;//!<入口密度
	double si;//!<入口エントロピ
	double Ti;//!<入口温度
	double Twi;//!<入口湿球温度
	double xi;//!<入口乾き度
	double Vi;//!<入口体積流量
	double Qi;//!<熱量[kW]
	double SHi;//!<過熱度[degC]
	double SCi;//!<サブクール度[degC]

	double Gj;//!<入口流量[kg/s]
	double Pj;//!<入口圧力[kPs]
	double hj;//!<入口比エンタルピ[kJ/kg]
	double Xj;//!<入口濃度
	double Dj;//!<入口密度
	double sj;//!<入口エントロピ
	double Tj;//!<入口温度
	double Twj;//!<入口湿球温度
	double xj;//!<入口乾き度
	double Vj;//!<入口体積流量
	double Qj;//!<熱量[kW]
	double SHj;//!<過熱度[degC]
	double SCj;//!<サブクール度[degC]

	double Go;//!<出口流量[kg/s]
	double Po;//!<出口圧力[kPa]
	double ho;//!<出口比エンタルピ[kJ/kg]
	double Xo;//!<出口濃度
	double Do;//!<出口密度
	double so;//!<出口エントロピ
	double To;//!<出口温度
	double Two;//!<出口湿球温度
	double xo;//!<出口乾き度
	double Vo;//!<出口体積流量
	double Qo;//!<熱量[kW]
	double SHo;//!<過熱度[degC]
	double SCo;//!<サブクール度[degC]

	double Gp;//!<出口流量[kg/s]
	double Pp;//!<出口圧力[kPa]
	double hp;//!<出口比エンタルピ[kJ/kg]
	double Xp;//!<出口濃度
	double Dp;//!<出口密度
	double sp;//!<出口エントロピ
	double Tp;//!<出口温度
	double Twp;//!<出口湿球温度
	double xp;//!<出口乾き度
	double Vp;//!<出口体積流量
	double Qp;//!<熱量[kW]
	double SHp;//!<過熱度[degC]
	double SCp;//!<サブクール度[degC]

	double Vi_m3ph;//!<入口流量[m3/h]
	double Vj_m3ph;//!<入口流量[m3/h]
	double Vo_m3ph;//!<出口流量[m3/h]
	double Vp_m3ph;//!<出口流量[m3/h]

	double Vi_m3ps;//!<入口流量[m3/s]
	double Vj_m3ps;//!<入口流量[m3/s]
	double Vo_m3ps;//!<出口流量[m3/s]
	double Vp_m3ps;//!<出口流量[m3/s]


	double Vi_lpm;//!<入口流量[l/min]
	double Vj_lpm;//!<入口流量[l/min]
	double Vo_lpm;//!<出口流量[l/min]
	double Vp_lpm;//!<出口流量[l/min]

	double PV;//!<気体の圧力[kPs]
	double hV;//!<気体の比エンタルピ[kJ/kg]
	double XV;//!<気体の濃度
	double DV;//!<気体の密度
	double sV;//!<気体のエントロピ
	double uV;//!<気体の内部エネルギー
	double TV;//!<気体の温度
	double MV;//!<気体の質量

	double PL;//!<気体の圧力[kPs]
	double hL;//!<気体の比エンタルピ[kJ/kg]
	double XL;//!<気体の濃度
	double DL;//!<気体の密度
	double sL;//!<気体のエントロピ
	double uL;//!<気体の内部エネルギー
	double TL;//!<気体の温度
	double ML;//!<気体の質量

	double Xtt;
	double gxi;
	double Hgxi;
	double phiv; 

	double ep;
	double W;

	double Re;//!<レイノルズ数

	double Re1;//!<レイノルズ数
	double Re2;//!<レイノルズ数
	double Re3;//!<レイノルズ数
	double ReL;//!<レイノルズ数
	double ReV;//!<レイノルズ数
	double ReLO;


	double Pr1;//プラントル数
	double Pr2;//プラントル数
	double Pr3;//プラントル数
	double PrL;//プラントル数
	double PrV;//プラントル数

	double Bo1;//!<ボイリング数
	double Bo2;//!<
	double Bo3;//!<
	double BoL;//!<
	double BoV;//!<


	double Nu1;
	double Nu2;
	double Nu3;
	double NuF;
	double NuB;

	double Ga1;
	double Ga2;
	double Ga3;
	double Ga4;


	double C_cavallini;
	double C1_cavallini;
	double A_cavallini;
	double Fr_cavallini;
	double Rx_cavallini;

	double Jg_cavallini;
	double Jgstar_cavallini;

	char IOCode;//!<0:出口 1:入口
	int  IONum;//!<出入口の番号
	char *FluidCode;//!<動作冷媒の名前"AMMONIA.FLD"や"R134a"など
	int  RouteCode;//!<経路の番号

	int ConnectModule;
	int ConnectPort;
	CFluidParameter *CP;


	double MaxStep;

	void setOutletConnection( int _Module , int _Port );

	void copy( CFluidParameter );
	bool valueCheck( CFluidParameter );

	void subset( DXFluid );
	void subset_i( DXFluid );
	void subset_j( DXFluid );
	void subset_o( DXFluid );
	void subset_p( DXFluid );

};


#endif
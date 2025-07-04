/*

M`BฆฦณอนธvZvO
hantei == 0 ฌฬ๐มMตฝ๊
hantei == 1 ฌฬ๐โpตฝ๊


gpตฝฎ

Pๆ
ณอนธ@F
M`B@@FDittus_Boelter

๑ๆ
ณอนธ@F
ฆซM`BFGungor_Winterton
รkM`BFFujii


*/

#ifndef __CPressureDropAndHeatTransfer_h
#define __CPressureDropAndHeatTransfer_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
using namespace std;

#include "DXProperty_ver06.h"
#include "CFluidParameter.h"
#include "CNewtonRaphsonMethod.h"
#include "PropertyWater_ver05.h"
#include "PropertyLithiumBromide_ver05.h"
#include "CStructureParameter.h"



#include "CPDHT0001DittusBoelter.h"
#include "CPDHT0002Blasius.h"
#include "CPDHT0003carnavos.h"
#include "CPDHT0004ito.h"
#include "CPDHT0005yoshida.h"
#include "CPDHT0006chisholm.h"
#include "CPDHT0007fujii.h"
#include "CPDHT0008nozu.h"
#include "CPDHT0009GungorWinterton.h"
#include "CPDHT0010SeshimoFujii.h"
#include "CPDHT0011InsideSmoothTube.h"
#include "CPDHT0012InsideGroovedTube.h"
#include "CPDHT0013Sawai.h"
#include "CPDHT0014Kubota.h"
#include "CPDHT0015Yonemoto_ht.h"
#include "CPDHT0016Mori.h"
#include "CPDHT0017Yonemoto_pd.h"
#include "CPDHT0018Nusselt.h"
#include "CPDHT0019Meisenburg.h"
#include "CPDHT0020Fujita.h"
#include "CPDHT0021Islam.h"
#include "CPDHT0022JeongGarimella.h"
#include "CPDHT0023Kamoshida.h"
#include "CPDHT0024StephanAbdelsalam.h"
#include "CPDHT0025Goto.h"
#include "CPDHT0027Chen.h"
#include "CPDHT0028Kandlikar.h"
#include "CPDHT0029Shah.h"
#include "CPDHT0030JongTaekOh.h"
#include "CPDHT0031LiuWinterton.h"
#include "CPDHT0032gnielinski.h"
#include "CPDHT0033Thrope.h"
#include "CPDHT0034Nishikawa.h"

/*
class chisholm{
public :

	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double D;//ผa[m]

	double S;//fสฯ[m3]
	double u;//ฌฌ[m/s]

	double G2;
	double f;
	double rho;
	double mu;
	double x;
	double Re;

	double rho_l;
	double v_l;
	double mu_l;
	double rho_v;
	double v_v;
	double mu_v;

	double G2_l;
	double G2_v;

	double f_l;
	double f_v;

	double Re_l;
	double Re_v;


	double phai2_l;
	double phai2_v;

	double dpdz_l;
	double dpdz_v;

	double Xtt;


	double dpdz;
	int kanetu;

	DXProperty ref;
	void calc();

};
*/

/*
class itou_pd{
public:
	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double D;//ผa[m]

	double S;//fสฯ[m3]
	double u;//ฌฌ[m/s]

	double G2;
	double f;
	double rho;
	double mu;
	double mu_v;
	double mu_l;
	double x;
	double Re;
	double R;
	double alpha;
	double ramda;
	double theta;

	double dpdz;
	int kanetu;

	DXProperty ref;
	void calc();

};*/

/*
class kinsitsuryu{
public:
	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double D;//ผa[m]

	double S;//fสฯ[m3]
	double u;//ฌฌ[m/s]

	double G2;
	double f;
	double rho;
	double mu;
	double mu_v;
	double mu_l;
	double x;
	double Re;

	double dpdz;
	int kanetu;

	DXProperty ref;
	void calc();

};
*/
/*
class kinsitsuryu2{
public:
	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double Dn;//๖ฤเa[m]
	double Dh;//อผa[m]
	double Da;//ฝฯผa[m]

	double Sfa;//ฌHฬภfสฯ[m2]
	double u;//ฌฌ[m/s]

	double beta;//๙j[rad]
	double G2;
	double f;
	double rho;
	double mu;
	double mu_v;
	double mu_l;
	double x;
	double Re;

	double Aef;

	double dpdz;
	int kanetu;

	DXProperty ref;
	void calc();

};*/

/*
class yoshida{
public:

	yoshida(){
		mnm.setup( 3 , 1.0+1e-3 , 1e-8 );

	}

	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double D;//ผa[m]

	double S;//fสฯ[m3]
	double u;//ฌฌ[m/s]
	double mu;//Sx[Pa s]
	double lambda;//M`ฑฆ[W/(m*K)]
	double Re;//CmY
	double Pr;//vg
	double Cp;//ไM[J/(kg*K)]
	double rho;//งx[kg/m3]
	double n;//
	double alpha;//

	double T_wi;
	double lambda_l;
	double T;
	double G2;
	double rho_l;
	double rho_v;
	double mu_l;
	double mu_v;
	double v_l;
	double v_v;
	double x;
	double hv;
	double Cp_l;
	double Pr_l;
	double Re_l;
	double alpha_l;
	double Bo;
	double _Bo;
	double E_new;
	int kanetu;
	double Xtt;

	CNewtonRaphsonMethod mnm;


	DXProperty ref;
	void calc();
};
*/
/*
class dittus_boelter{
public:
	CFluidParameter		*fp;
	CStructureParameter *sp;

	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double T;//ทx[degC]
	double D;//ผa[m]

	double S;//fสฯ[m3]
	double u;//ฌฌ[m/s]
	double mu;//Sx[Pa s]
	double lambda;//M`ฑฆ[W/(m*K)]
	double Re;//CmY
	double Pr;//vg
	double Cp;//ไM[J/(kg*K)]
	double rho;//งx[kg/m3]
	double n;//
	double alpha;//M`Bฆ[kJ/(m2*K)]
	int kanetu;

	DXProperty ref;
	void calc();
	void calc2();
};
*/
/*
class gungor_winterton{
public:
	CFluidParameter		*fp;
	CStructureParameter *sp;

	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double D;//ผa[m]

	double S;//fสฯ[m3]
	double u;//ฌฌ[m/s]
	double mu;//Sx[Pa s]
	double lambda;//M`ฑฆ[W/(m*K)]
	double Re;//CmY
	double Pr;//vg
	double Cp;//ไM[J/(kg*K)]
	double rho;//งx[kg/m3]
	double n;//
	double alpha;//

	double T_wi;
	double lambda_l;
	double T;
	double G2;
	double rho_l;
	double rho_v;
	double mu_l;
	double mu_v;
	double x;
	double hv;
	double hfg;
	double Cp_l;
	double Pr_l;
	double Re_l;
	double alpha_l;
	double Bo;
	double _Bo;
	double E_new;
	int kanetu;


	DXProperty ref;
	void calc();
	void calc2();
};*/

/*
class fujii{
public:
	CFluidParameter		*fp;
	CStructureParameter *sp;

	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double D;//ผa[m]

	double g;//dอมฌx
	double G2;//โ}ฟสฌx[kg/(m2*s)]


	double S;//fสฯ[m3]
	double u;//ฌฌ[m/s]
	double mu;//Sx[Pa s]
	double lambda;//M`ฑฆ[W/(m*K)]
	double Re;//CmY
	double Pr;//vg
	double Nu;//
	double Cp;//ไM[J/(kg*K)]
	double rho;//งx[kg/m3]
	double n;//
	double alpha;//

	double Nu_f;
	double Nu_b;
	double Re_l;
	double v_g;
	double v_l;
	double x;
	double Pr_l;
	double C1;
	double C2;
	double C3;
	double C4;
	double C5;
	double H;
	double Ga;
	double W;
	double mu_l;
	double mu_v;
	double T_sat;
	double lambda_l;
	
	double rho_l;
	double rho_v;
	double taw_w;
	double taw_wv;
	double fai_v;
	double Cp_l;
	double T_i;
	double T_wi;//ูเวสทx
	double q;
	double X_tt;
	double gzai;
	double A;
	double H_l;
	double hv;
	double hfg;

	double Re_lstar;
	double T_iplus;

	double m;

	int kanetu;

	DXProperty ref;
	void calc();
	void calc2();
};

*/
/*

class fujii2{
public:

	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double D;//ผa[m]

	double g;//dอมฌx
	double G2;//โ}ฟสฌx[kg/(m2*s)]


	double S;//fสฯ[m3]
	double u;//ฌฌ[m/s]
	double mu;//Sx[Pa s]
	double lambda;//M`ฑฆ[W/(m*K)]
	double Re;//CmY
	double Pr;//vg
	double Ja;
	double Ja_v;
	double Nu;//
	double Cp;//ไM[J/(kg*K)]
	double rho;//งx[kg/m3]
	double n;//
	double alpha;//
	double superheat;//

	double Nu_f;
	double Nu_b;
	double Re_l;
	double v_g;
	double v_l;
	double x;
	double Pr_l;
	double Pr_v;
	double C1;
	double C2;
	double C3;
	double C4;
	double C5;
	double H;
	double Ga;
	double W;
	double mu_l;
	double mu_v;
	double T_sat;
	double lambda_l;
	double lambda_v;
	
	double rho_l;
	double rho_v;
	double taw_w;
	double taw_wv;
	double fai_v;
	double Cp_l;
	double Cp_v;
	double T_i;
	double T_wi;//ูเวสทx
	double q;
	double X_tt;
	double gzai;
	double A;
	double H_l;
	double hv;

	double Re_lstar;
	double T_iplus;

	double m;

	int kanetu;

	DXProperty ref;
	void calc();
};
*/

/*class fujita{
public:
	CFluidParameter		*fp;
	CStructureParameter *sp;

	double gamma;	//ะคฌส[kg/ms]
	double Re_f;	//CmY
	double Nu_m;	//kbZg
	double g;		//dอมฌx

	void calc();
	void calc2();
};*/
/*class nusselt{
public:
	CFluidParameter		*fp;
	CStructureParameter *sp;

	double h_fg;	//รk๖M
	double g;		//dอมฌx

	void calc();
	void calc2();
};*/

/*class islam{
public:
	CBasicParameter B;
	CFluidParameter		*fp;
	CStructureParameter *sp;

	double gamma;	//ะคฌส[kg/ms]
	double Ja;
	double Ka;
	double Sc;
	double Sh;
	double Po;		//๎ณอ(=1.0)[kPa]
	double i_ab;
	double g;
	double b;

	void calc();
	void calc2();
};*/

/*class garimella{
	PropertyLithiumBromide libr;
public:
	CFluidParameter		*fp;
	CStructureParameter *sp;

	double X_b;
	double X_ii;

	double g;
	double pi;
	double c1;
	double c2;
	double c3;
	double c4;
	double da;
	double dd;
	double m_d;
	double n;
	double t_form;
	double WR;
	double X_i;
	double delta_c;
	double gamma;	//ะคฌส[kg/ms]
	double lambda;
	double xi;

	void calc();
	void calc2();
	void calc_form();
};*/


/*
class seshimo{
public:

	//๓CM`Bฆ

	CBasicParameter B;

	//	<input>
//	double AO;
	double Di;
	double Do;
	double L_x;
	double L_y;
	double L_y_pitch;
	double L;
	double pNum_y;

	double Fp;//tBsb`
	double Ft;//tB๚ณ

	double L_f;
//	double ramda_f;
	double ramda_fin;
//	double eps;

//	double T_air;
	double rho_air;
	double rho_airi;
	double visc_air;
	double v_air;
	double thc_air;
	double Pr_air;
//	double cp_air;

	//ผa
	double Df;
	//tBJ[a
	double Dc;
	//tB`Mสฯ
	double Ae;
	//ว`Mสฯ
	double Ap;
	//S`Mสฯ
	double Ao;
	//Sเ`Mสฯ
	double Ai;
	//Oสสฯ
	double Af;
	//ฝฯฉRส฿fสฯ
	double Ac;
	//ใ\ฌ
	double Vac;
	//ใ\ก@
	double Dec;
	//CmY
	double Re_air;
	//kbZgi๊๑j
	double Nu;
	//M`Bฆ
	double alpha;
	//fin๘ฆ
	double e_fin;


	double De_min;
	double Va_max;
	double Re_air_star;

	//ณน
	double dpdz;
	double f;

	DXProperty air;
	void calc();
};

*/

/*
class Grooved_tube{
public:

	CFluidParameter *fp;
	CStructureParameter *sp;

	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double D;//ผa[m]
	double T;//ทx[degC]
	double T_wi;//วสทx[degC]

	double x;

	double alpha;
	double dpdz;

	DXProperty ref;

	CPDHT0001DittusBoelter dbo;
	CPDHT0002Blasius blas;
	CPDHT0005yoshida ysh;
	CPDHT0006chisholm chi;
	CPDHT0008nozu noz;


	int kanetu;

	void init();
	void setup();
	void calc_pd();
	void calc_ht();



};
*/


class CPressureDropAndHeatTransfer{
public:

	CFluidParameter fp;
	CFluidParameter _fp;
	CStructureParameter sp;

	double G;//ฟสฌส[kg/s]
	double P;//ณอ[kPa]
	double h;//ไG^s[kJ/kg]
	double D;//ผa[m]
	double T;//ทx[degC]
	double T_wi;//วสทx[degC]

	double x;

	double alpha;
	double dpdz;

	DXProperty ref;
	PropertyWater wat2;
	PropertyLithiumBromide libr;

	CPDHT0001DittusBoelter dbo;
	CPDHT0002Blasius blas;
	CPDHT0003carnavos car;
	CPDHT0004ito ito;
	CPDHT0005yoshida ysh;
	CPDHT0006chisholm chi;
	CPDHT0007fujii fji;
	CPDHT0008nozu noz;
	CPDHT0009GungorWinterton gwi;
	CPDHT0010SeshimoFujii ssm;
	CPDHT0011InsideSmoothTube IST;
	CPDHT0012InsideGroovedTube Grv;
	CPDHT0013Sawai Saw;
	CPDHT0014Kubota Kub;
	CPDHT0015Yonemoto_ht Yoh;
	CPDHT0016Mori Mor;
	CPDHT0017Yonemoto_pd Yon;
	CPDHT0018Nusselt nsl;
	CPDHT0019Meisenburg msb;
	CPDHT0020Fujita fjt;
	CPDHT0021Islam isl;
	CPDHT0022JeongGarimella grm;
	CPDHT0023Kamoshida kms;
	CPDHT0024StephanAbdelsalam sab;
	CPDHT0025Goto got;
	CPDHT0027Chen chn;
	CPDHT0028Kandlikar knd;
	CPDHT0029Shah sha;
	CPDHT0030JongTaekOh jon;
	CPDHT0031LiuWinterton liu;
	CPDHT0032gnielinski gni;
	CPDHT0033Thrope Thr;
	CPDHT0034Nishikawa Nis;
	//fujita fjt;
	//nusselt nsl;
	//islam isl;
	//garimella grm;

	int kanetu;

	void init();
	void setup();
	void setup_water(double P, double h);
	void setup_libr(double P, double h, double X);
//	void calc_pd();
//	void calc_ht();
	void calc_ht2();

};



#endif
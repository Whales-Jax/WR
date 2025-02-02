#ifndef _PropertyWater_ver05_h_
#define _PropertyWater_ver05_h_


#include <cmath>



class PropertyWater{
public:

	PropertyWater();
	~PropertyWater();

	double Tc;
	double Pc;
	double vc;
	double I1;
	double T0;
	double L0;
	double L1;
	double L2;
	double B0;
	double alpha0;
	double alpha1;

	int n[9],z[9][4],l[9],x[9][3];
	double A[23],a[13];
	double bb,B[10][7],b[9][9];
	double k[10];


	double k_func  ( double _T );
	double L_func  ( double _T );
	double Ld_func ( double _T );
	double chi1    ( double _P , double _T );
	double sigma1  ( double _P , double _T );

	double epsilon1  ( double _P , double _T );
	double chi2      ( double _P , double _T );
	double sigma2    ( double _P , double _T );
	double epsilon2  ( double _P , double _T );
	double phi1      ( double _P , double _T );
	double phi2      ( double _P , double _T );
	double troumyuu  ( double _T , double _rho );
	double troulambda( double _T , double _rho );
	double t_p       ( double _T );
	double p_t       ( double _P );

	double sc_vl( double _P , double _T );
	double sh_vv( double _P , double _T );
	double sat_vl(double t);
	double sat_vv(double t);

	double sc_roul(double p,double t);
	double sh_rouv(double p,double t);
	double sat_roul(double t);
	double sat_rouv(double t);

	double sc_sl(double p,double t);
	double sh_sv(double p,double t);
	double sat_sl(double t);
	double sat_sv(double t);

	double sc_hl(double p,double t);
	double sh_hv(double p,double t);
	double sat_hl(double t);
	double sat_hv(double t);

	double sc_cpl(double p,double t);
	double sh_cpv(double p,double t);
	double sat_cpl(double t);
	double sat_cpv(double t);

	double sc_myul(double p,double t);
	double sh_myuv(double p,double t);
	double sat_myul(double t);
	double sat_myuv(double t);

	double sc_laml(double p,double t);
	double sh_lamv(double p,double t);
	double sat_laml(double t);
	double sat_lamv(double t);

	double sc_nyul(double p,double t);
	double sh_nyuv(double p,double t);
	double sat_nyul(double t);
	double sat_nyuv(double t);

	double sc_Prl(double p,double t);
	double sh_Prv(double p,double t);
	double sat_Prl(double t);
	double sat_Prv(double t);

	double sc_el(double p,double t,double p0=101.325,double t0=20.0);
	double sh_ev(double p,double t,double p0=101.325,double t0=20.0);
	double sat_el(double t,double p0=101.325,double t0=20.0);
	double sat_ev(double t,double p0=101.325,double t0=20.0);

	double sc_tl(double p,double h);
	double sh_tv(double p,double h);


	double rou_Ph( double P , double h );
	double T_Ph( double P , double h );







};


#endif
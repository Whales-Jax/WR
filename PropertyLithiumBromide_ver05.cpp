#include "PropertyLithiumBromide_ver05.h"




PropertyLithiumBromide::PropertyLithiumBromide(){

	A[0] = -2024.33;
	A[1] = 163.309;
	A[2] = -4.88161;
	A[3] = 6.302948E-2;
	A[4] = -2.913705E-4;

	B[0] = 18.2829;
	B[1] = -1.1691757;
	B[2] = 3.248041E-2;
	B[3] = -4.034184E-4;
	B[4] = 1.8520569E-6;

	C[0] = -3.7008214E-2;
	C[1] = 2.8877666E-3;
	C[2] = -8.1313015E-5;
	C[3] = 9.9116628E-7;
	C[4] = -4.4441207E-9;

	D[0] = -2.00755;
	D[1] = 0.16976;
	D[2] = -3.133362e-3;
	D[3] = 1.97668e-5;

	E[0] = 124.937;
	E[1] = -7.71649;
	E[2] = 0.152286;
	E[3] = -7.95090e-4;

	dec=1.0+1e-7;
	eps=1.0e-10;

}

PropertyLithiumBromide::~PropertyLithiumBromide(){
}


double PropertyLithiumBromide::sc_h_XT(double _X , double _T ){
	_X *= 100.0;
	An = 0;
	for( int i = 0 ; i < 5 ; i++ ){
		An += A[i] * pow( _X , i );
	}
	Bn = 0;
	for( int j = 0 ; j < 5 ; j++ ){
		Bn += B[j] * pow( _X , j );
	}
	Cn = 0;
	for( int k = 0 ; k < 5 ; k++ ){
		Cn += C[k] * pow( _X , k );
	}
	h  = An + _T * Bn + pow( _T , 2 ) * Cn;
	return( h );
}

double PropertyLithiumBromide::sc_T_Xh(double _X , double _h ){
	_X *= 100;
	An = 0;
	for( int i = 0 ; i < 5 ; i++ ){
		An += A[i] * pow( _X , i );
	}
	Bn = 0;
	for( int j = 0 ; j < 5 ; j++ ){
		Bn += B[j] * pow( _X , j );
	}
	Cn = 0;
	for( int k = 0 ; k < 5 ; k++ ){
		Cn += C[k] * pow( _X , k );
	}
	T = ( - Bn + pow( (pow(Bn,2)-4*Cn*(An-_h)),0.5) ) / 2.0 / Cn;
	return (T);
}

double PropertyLithiumBromide::sc_T_XTsat(double _X , double _Tsat ){
	_X *= 100;
	Dn = 0;
	for( int i = 0 ; i < 4 ; i++ ){
		Dn += D[i] * pow( _X , i );
	}
	En = 0;
	for( int j = 0 ; j < 4 ; j++ ){
		En += E[j] * pow( _X , j );
	}
	T = _Tsat * Dn + En;
	return( T );
}

double PropertyLithiumBromide::sc_X_TTsat(double _T , double _Tsat){
	X = 0.50;
	do{
		for(int i=0;i<=1;i++)
		{
			F[i] = X*((i==0)?(dec):(1.0));
			G[i] = sc_T_XTsat( F[i] , _Tsat ) - _T;
		}
		X = F[0]-G[0] * (F[1]-F[0]) / (G[1]-G[0]);
	}while( fabs( G[1] ) >= eps );
	return(X);
}

double PropertyLithiumBromide::sc_Tsat_XT(double _X , double _T){
	int i=0;
	Tsat = 50.0;
	do{
		i++;
		if(i>1000){return(999);}

		for(int j=0;j<=1;j++)
		{
			F[j] = Tsat*((j==0)?(dec):(1.0));
			G[j] = sc_T_XTsat( _X , F[j] ) - _T;
		}
		Tsat = F[0]-G[0] * (F[1]-F[0]) / (G[1]-G[0]);
	}while( fabs( G[1] ) >= eps );
	return(Tsat);
}

double PropertyLithiumBromide::sc_rho_XT( double _X , double _T ){
	_X *= 100;
	a = 1169.2 + 4.4304 * _X + 0.12096 * pow ( _X, 2 );
	b = -0.48269 - 6.7967e-4 * _X;
	rho = a + b * ( _T + 273.15 );
	return ( rho );
}

double PropertyLithiumBromide::sc_X_Trho(double _T , double _rho){
	int i=0;
	X = 0.50;
	do{
		i++;
		if(i>1000){return(999);}

		for(int j=0;j<=1;j++)
		{
			F[j] = X*((j==0)?(dec):(1.0));
			G[j] = sc_rho_XT( F[j] , _T ) - _rho;
		}
		X = F[0]-G[0] * (F[1]-F[0]) / (G[1]-G[0]);
	}while( fabs( G[1] ) >= eps );
	return(X);
}

double PropertyLithiumBromide::sc_visc_XT( double _X , double _T ){
   T = 273.15 + _T;
   _X = _X * 100.0;
   double A1, A2, A3;
   A1 = - 494.122 + 16.3967 * _X - 0.14511 * pow( _X , 2 );
   A2 =   28606.4 - 934.568 * _X + 8.52755 * pow( _X , 2 );
   A3 =   70.3848 - 2.35014 * _X + 0.0207809 * pow( _X , 2 );
   visc = exp( A1 + A2 / T + A3 * log( T ) )/ 1000.0;

   return( visc );
}

double PropertyLithiumBromide::sc_visc_XT2( double _X , double _T ){
   t10 = 323.0 / ( 273.15 + _T );
   _X = _X * 100.0;
   _X = _X / 30;
   a10 =  5.0783 - 34.995 * _X + 78.000 * pow( _X , 2 ) - 91.150 * pow( _X , 3 ) + 50.531 * pow( _X , 4 ) - 10.530 * pow( _X , 5 );
   b10 = -16.242 + 66.917 * _X - 144.34 * pow( _X , 2 ) + 168.34 * pow( _X , 3 ) - 93.511 * pow( _X , 4 ) + 19.579 * pow( _X , 5 );
   c10 =  10.610 - 32.023 * _X + 67.893 * pow( _X , 2 ) - 78.856 * pow( _X , 3 ) + 43.789 * pow( _X , 4 ) - 9.1576 * pow( _X , 5 );
   visc  = pow( 2.7183 , ( a10 + b10 * t10 + c10 * t10 * t10 ) ) / 1000.0;

   return( visc );
}

double PropertyLithiumBromide:: sc_thc_XT( double _X, double _T )
{
	_T = _T + 273.15;
	_X *= 100.0;
	double A, B, C;
	A = - 1407.53  + 11.0513    * _T - 1.46741e-2 * pow( _T , 2 );
	B =   38.9855  - 0.240475   * _T + 3.48073e-4 * pow( _T , 2 );
	C = - 0.265025 + 1.51915e-3 * _T - 2.32262e-6 * pow( _T , 2 );
	thc = ( A + B * _X + C * pow( _X, 2 ) ) / 1000.0; 
	return thc;
}

double PropertyLithiumBromide:: sc_thc_XT2( double _X, double _T )
{
	_T = _T + 273.15;
	_X *= 100.0;
	thc = -8.9012e-1 + 9.0301e-3 * _X
		 + _T * (  8.3279e-3 - 7.2638e-5 * _X )
		 + pow( _T, 2 ) * ( -1.0937e-5 + 1.0213e-7 * _X );
	return thc;
}

double PropertyLithiumBromide:: sc_cp_XT ( double _X, double _T )
{
	double A, B;

	A = 3.067819 - 2.15232 * _X;
	B = 0.006018 - 0.00731 * _X;
	
	cp = A + B * _T;

	return cp;
}

double PropertyLithiumBromide:: sc_cp_XT2 ( double _X, double _T )
{
	_X *= 100.0;
	cp = 1.098 - 1.529e-2 * _X + 6.220e-5 * pow( _X , 2.0 )
		+ _T * ( -3.651e-3 + 4.204e-5 * _X )
		+ pow ( _T, 2.0 ) * ( 3.576e-5 - 4.238e-7 * _X );
	cp *= 4.18605;
	return cp;
}

double PropertyLithiumBromide:: sc_st_XT (double _X, double _T)
{
	_X *= 100.0;

	st =                       7.626234e+1 + _T * ( -1.507474e-1 - _T * 1.107075e-5  )
		+ pow ( _X, 1.0 ) * (  4.583900e-1 + _T * ( -9.057263e-3 + _T * 7.238986e-5  ) )
		+ pow ( _X, 2.0 ) * ( -1.463071e-2 + _T * (  4.459087e-4 - _T * 3.822731e-6  ) )
		+ pow ( _X, 3.0 ) * (  3.834735e-4 + _T * ( -9.542318e-6 + _T * 8.077592e-8  ) )
		+ pow ( _X, 4.0 ) * ( -2.733854e-6 + _T * (  6.610416e-8 - _T * 5.681625e-10 ) );

	st = st / 1000.0;//’PˆÊ‚ÌŠ·ŽZ

	return st;
}

double PropertyLithiumBromide:: sc_st_XT2 (double _X, double _T)
{
	_X = _X * 100.0;

	E0 = 1.50949e-5*_T*_T - 1.11758e-3*_T + 0.340993;
	E1 = -3.08467e-7*_T*_T - 5.48023e-7*_T + 8.40993e-3;
	E2 = 7.76341e-5*_T*_T - 0.158370*_T + 76.1867;
	
	st = (exp(E0 + E1*_X) + E2 -30.0) * 1e-3;
	
	return st;
}

double PropertyLithiumBromide:: sc_d_X (double _X)
{
	_X = _X * 100.0;
	d = 2.46e-9 - 2.7e-11*_X;
	return d;
}

double PropertyLithiumBromide:: sc_d_XT (double _X, double _T)
{
	double mu_sl;
	double mu_sl_25;
	double A_ll;
	double D_sl;

	mu_sl = sc_visc_XT( _X , _T );
	mu_sl_25 = sc_visc_XT( _X , 25.0 );

	_X = _X * 100.0;

	A_ll =  1.106 * 10000.0 + 4.684 * 100.0 * _X - 1.062 * 10.0 * _X*_X + 3.468/100.0 * _X*_X*_X;
	D_sl = 3.6*1e-10 * A_ll * ((_T+273.15)/298.15) * (mu_sl_25/mu_sl) / 3600.0;
	return D_sl;

	
}




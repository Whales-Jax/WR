/***********************************************************
 * Air_Property.cpp
 *
 * 2008/07 ������O
 ***********************************************************/
/*
	��C�����̌v�Z�́wM.CONDE ENGINEERING, 2007�x���Q�Ƃ��邱��
	������C�Ɛ����C�̔M�`�����ƔS���W���̎��͏o�����킩�炸���x���s���Ȃ̂ŁC
	�ʂ̂��̂��g�p�D
	������C��
		Journal of Physical and chemical reference data, 1985, 14,(4), P947-970
		�wViscosity and Thermal Conductivity of Dry Air in the Gaseous Phase�x
		K. Kadoya, N. Matsunaga, A. Nagashima
	���Q��
	�����C��
		�wproperty_water_ver2.h�x	�i1980 ���{�@�B�w����C�\�j
	���Q�Ƃ��邱��
	���̕����Ɋւ��Ă͐V�������ۊ�ɂȂ��Ă��Ȃ��̂ŁC�O�ꂵ�Ă�肽���l�͍ŐV
	�̓��{�@�B�w����C�\������Ύ����킩��܂��D�܂��C���̕�����0.0�x�ȉ��ł͎g
	���Ȃ��̂ŁC�v�Z�͈͂����l��refprop���g�����Ƃ��������ł��D�ł��C�M�`������
	�S���W�����v�Z���邽�߂����Ƀt�@�C����������̂������Ȃ̂ŁC���̂Ƃ���͎���
	����C�͂Ȃ��ł��D

	���jproperty_water_ver2.h���g�p���Ă���̂�include����̂�Y�ꂸ�ɁI�I�I
*/

// 07/08 �ꕔ�ύX by ����

#include <iostream>
#include <cmath>
#include "Moist_Air_Property.h"
#include "property_water_ver2.h"

const double Error = 1e-9; // ��������
const double Delta = 1.0 + 1e-3; // ���ݕ�
const int Loop = 1000; // �J��Ԃ��v�Z��
int LoopCount;

const double Ma = 28.9645; // ��C�̕��q�� g/mol
const double Mw = 18.01528; // ���̕��q�� g/mol
const double R = 8.31441; // �C�̒萔 J/(mol/K)

/********************************************************************************
 * Air::Air_Ptx -- �C�ӏ�Ԋ֐�(���́E���x�E��Ύ��x����)
 *                 �C�ӂ̈��́E���x�E��Ύ��x�ɂ������Ԃ̎Z�o
 *
 * ����(���͒l)
 *     pressure     -- ���� kPa
 *     temperature  -- ���x deg.C
 *     humidity     -- ��Ύ��x kg/kg(DA)
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
 * Air::Air_Phx -- �C�ӏ�Ԋ֐�(���́E��G���g���s�[�E��Ύ��x����)
 *                 �C�ӂ̈��́E��G���g���s�[�E��Ύ��x�ɂ������Ԃ̎Z�o
 *
 * ����(���͒l)
 *     pressure    -- ���� kPa
 *     enthalpy    -- ��G���^���s�[ kJ/kg(DA)
 *     humidity    -- ��Ύ��x kg/kg(DA)
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
 * Air::Air_tx -- �C�ӏ�Ԋ֐�(���x�E��Ύ��x����)
 *                 �C�ӂ̉��x�E��Ύ��x�ɂ������Ԃ̎Z�o
 *
 * ����(���͒l)
 *     temperature  -- ���x deg.C
 *     humidity     -- ��Ύ��x kg/kg(DA)
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
 * Air::Air_hx -- �C�ӏ�Ԋ֐�(��G���g���s�[�E��Ύ��x����)
 *                 �C�ӂ̔�G���g���s�[�E��Ύ��x�ɂ������Ԃ̎Z�o
 *
 * ����(���͒l)
 *     enthalpy    -- ��G���^���s�[ J/kg(DA)
 *     humidity     -- ��Ύ��x kg/kg(DA)
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
* �����C��e�ϊ֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F��e��[m^3/kg(DA)]
************************************************/
double Air_Ptxv ( double P, double T, double x )
{
	double v, v_;
	double na, nw;
	double Xa, Xw;
	double Bm, Cm;

	// �P�ʕϊ�
	P *= 1e3; // kPa -> Pa
	T += 273.15; // deg.C -> K

	// ����
	na = 1.0 / Ma;
	nw = x / Mw;

	// ��������
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// Bm, Cm�v�Z
	Bm_Cm ( T, Xa, Xw, Bm, Cm );
	
	// ��e�όv�Z
	// ����
	v = R * T / P;

	LoopCount = 0;
	while ( LoopCount < Loop )
	{
		v_ = R * T / P * ( 1.0 + Bm * 1e-6 / v + Cm * 1e-12 / ( v * v ) );
		
		// ��������
		if ( fabs ( 1.0 - v_ / v ) < Error ) break;

		// �l�̍X�V
		v = v_;

		LoopCount++;
	}

	// �P�ʊ��Z
	v = v_ / ( Xa * Ma * 1e-3 ); // m^3/mol -> m^3/kg(DA)

	return v;
}

/***********************************************
* �����C��e�ϊ֐�	-60 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ���Ύ��x[-]
* �o�́F��e��[m^3/kg(DA)]
************************************************/
double Air_Ptphiv ( double P, double T, double phi )
{
	double x;
	double v;

	// ��Ύ��x
	x = Air_Ptphix ( P, T, phi );

	// ��e��
	v = Air_Ptxv ( P, T, x );

	return ( v );
}

/***********************************************
* �����C���x�֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F���x[kg(DA)/m^3]
************************************************/
double Air_Ptxrho ( double P, double T, double x )
{
	double rho;

	// ���x
	rho = 1.0 / Air_Ptxv ( P, T, x );

	return ( rho );
}

/***********************************************
* �����C���x�֐�	-60 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ���Ύ��x[-]
* �o�́F���x[kg(DA)/m^3]
************************************************/
double Air_Ptphirho ( double P, double T, double phi )
{
	double rho;

	// ���x
	rho = 1.0 / Air_Ptphiv ( P, T, phi );

	return ( rho );
}

/***********************************************
* �����C��G���^���s�[�֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F��G���^���s�[[kJ/kg(DA)]
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

	// ������
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

	// ����
	na = 1.0 / Ma;
	nw = x / Mw;

	// ��������
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// ��e��v
	v = Air_Ptxv ( P, T, x );

	// �P�ʊ��Z
	v *= Xa * Ma * 1e-3 * 1e6; // m^3/kg(DA) -> cm^3/mol

	// �P�ʕϊ�
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

	// ��G���^���s�[
	hm = Xa * ( GT + ha ) + Xw * ( HT + hw ) + R * T * ( ( Bm - T * ( Bm_delta - Bm ) / ( T_delta - T ) ) / v 
		+ ( Cm - 0.5 * T * ( Cm_delta - Cm ) / ( T_delta - T ) ) / ( v * v ) );

	// �P�ʊ��Z
	hm /= ( Xa * Ma ); // J/mol -> kJ/kg(DA)

	return ( hm );
}

/***********************************************
* �����C��G���^���s�[�֐�	-60 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ���Ύ��x[-]
* �o�́F��G���^���s�[[kJ/kg(DA)]
************************************************/
double Air_Ptphih ( double P, double T, double phi )
{
	double x;
	double h;

	// ��Ύ��x
	x = Air_Ptphix ( P, T, phi );

	// ��G���^���s�[
	h = Air_Ptxh ( P, T, x );

	return ( h );
}

/***********************************************
* �����C��G���g���s�[�֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F��G���g���s�[[kJ/(kg(DA)*K)]
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

	// ������
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

	sa = -196.125465; // J/(mol�K)
	sw = -63.31449; // J/(mol�K)

	// ����
	na = 1.0 / Ma;
	nw = x / Mw;

	// ��������
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// ��e��v
	v = Air_Ptxv ( P, T, x );

	// �P�ʊ��Z
	v *= Xa * Ma * 1e-3 * 1e6; // m^3/kg(DA) -> cm^3/mol

	// �P�ʕϊ�
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

	// ��G���g���s�[
	sm = Xa * ( MT + M[5] * log ( T ) + sa ) + Xw * ( NT + N[6] * log ( T ) + sw ) - R * log ( P / 101325.0 ) + Xa * R * log ( ( P * v ) / ( Xa * R0 * T ) )
		+ Xw * R * log ( ( P * v ) / ( Xw * R0 * T ) ) - R * ( ( Bm + T * ( Bm_delta - Bm ) / ( T_delta - T ) ) / v
		+ 0.5 * ( Cm + T * ( Cm_delta - Cm ) / ( T_delta - T ) ) / ( v * v ) );

	// �P�ʊ��Z
	sm /= ( Xa * Ma ); // J/(mol�K) -> kJ/(kg(DA)�K)

	return ( sm );
}

/***********************************************
* �����C��G���g���s�[�֐�	-60 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ���Ύ��x[-]
* �o�́F��G���g���s�[[kJ/(kg(DA)*K)]
************************************************/
double Air_Ptphis ( double P, double T, double phi )
{
	double x;
	double s;

	// ��Ύ��x
	x = Air_Ptphix ( P, T, phi );

	// ��G���g���s�[
	s = Air_Ptxs ( P, T, x );

	return ( s );
}

/***********************************************
* �����C���Ύ��x�֐�	-60 <= T <= 365
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F���Ύ��x[-]
************************************************/
double Air_Ptxphi ( double P, double T, double x )
{
	double phi;
	double Pvs;
	double f;
	double na, nw;
	double Xw;

	// ����
	na = 1.0 / Ma;
	nw = x / Mw;

	// ��������
	Xw = nw / ( na + nw );

	// �O�a�����C��
	Pvs = Air_tPs ( T );

	// f
	f = f_ ( P, T );

	// ���Ύ��x
	phi = Xw * P / ( f * Pvs );

	return ( phi );
}

/***********************************************
* �����C���Ύ��x�֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       �������x[deg.C]
* �o�́F���Ύ��x[-]
************************************************/
double Air_Pttwbphi ( double P, double T, double Twb )
{
	double phi;
	double Twb_;
	double dT;

	int i, j;
	double X[2], Y[2];

	// ���x��
	dT = T - Twb;

	if ( dT < 0.0 )
	{
		std::cout << std::endl;
		std::cout << "Twb exceeds T" << std::endl;
		exit(1);
	}

	// ����l
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

			// �������x
			Twb_ = Air_Ptphitwb ( P, T, phi );

			X[j] = phi;
			Y[j] = 1.0 - Twb_ / Twb;
		} // end of j loop
		
		// ��������
		if ( fabs ( Y[1] ) < Error ) break;
		else phi = X[1] - Y[1] * ( X[0] - X[1] ) / ( Y[0] - Y[1] );

		// �J��Ԃ��v�Z��
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
* �����C��Ύ��x�֐�	-60 <= T <= 365
* ���́F����[kPa]
*       ���x[deg.C]
*       ���Ύ��x[-]
* �o�́F��Ύ��x[kg/kg(DA)]
************************************************/
double Air_Ptphix ( double P, double T, double phi )
{
	double x;
	double Pvs;
	double f;
	double na;
	double Xw;

	// ����
	na = 1.0 / Ma;//��C1kg������̃����� ��C1mol�Ŗ�28kg���炢�i�\�[�X�̏���ɁA�����l�Ƃ��ē��́j

	// �O�a�����C��  M.CONDE ENGINEERING, 2007	�̎�[10]���
	Pvs = Air_tPs ( T );//

	// f
	f = f_ ( P, T );

	// �������� M.CONDE ENGINEERING, 2007	�̎�[9]���
	Xw = phi * f * Pvs / P;

	if ( Xw > 1.0 )
	{
		std::cout << std::endl;
		std::cout << "Pvs exceeds P -> Air_Ptphix" << std::endl;
		exit(1);
	}

	// ��Ύ��x
	x = na * Xw * Mw / ( 1.0 - Xw );

	return ( x );
}

/***********************************************
* �����C��Ύ��x�֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ��G���^���s�[[kJ/kg(DA)]
*       ���x[deg.C]
* �o�́F��Ύ��x[kg/kg(DA)]
************************************************/
double Air_Phtx ( double P, double h, double T )
{
	double x;
	double T_;
	int i, j;
	double X[2], Y[2];

	// ����
	x = ( h - 1.006 * T ) / ( 1.805 * T + 2501.0 );

	for ( i = 1 ; i <= Loop ; i++ )
	{
		for ( j = 0, x *= Delta ; j < 2 ; j++ )
		{
			if ( j == 1 )	x /= Delta;

			// ���x
			T_ = Air_Phxt ( P, h, x );

			X[j] = x;
			Y[j] = T - T_;
		} // end of j loop

		// ��������
		if ( fabs ( Y[1] ) < Error ) break;
		else x = X[1] - Y[1] * ( X[0] - X[1] ) / ( Y[0] - Y[1] );

		// �J��Ԃ��v�Z��
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
* �����C��Ύ��x�֐�	-60 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       �������x[deg.C]
* �o�́F��Ύ��x[kg/kg(DA)]
************************************************/
double Air_Pttwbx ( double P, double T, double Twb )
{
	double x;
	double phi;

	// ���Ύ��x
	phi = Air_Pttwbphi ( P, T, Twb );

	// ��Ύ��x
	x = Air_Ptphix ( P, T, phi );

	return ( x );
}

/***********************************************
* �����C���x�֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ��G���^���s�[[kJ/kg(DA)]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F���x[deg.C]
************************************************/
double Air_Phxt ( double P, double h, double x )
{
	double T;
	double h_;

	int i, j;
	double X[2], Y[2];

	// ����l
	T = ( h - 2501.0 * x ) / ( 1.006 + 1.805 * x );

	if( - 0.001 < T && T < 0.001 ){
		T = 0.001;
	}//����l��0�ɂȂ����Ƃ��Ɍv�Z�ł��Ȃ��o�O�C���@���

	for ( i = 1 ; i <= Loop ; i++ )
	{
		for ( j = 0, T *= Delta ; j < 2 ; j++ )
		{
			if ( j == 1 )	T /= Delta;

			// ��G���^���s�[
			h_ = Air_Ptxh ( P, T, x );

			X[j] = T;
			Y[j] = h - h_;
		} // end of j loop
	
		// ��������
		if ( fabs ( Y[1] ) < Error ) break;
		else T = X[1] - Y[1] * ( X[0] - X[1] ) / ( Y[0] - Y[1] );

		// �J��Ԃ��v�Z��
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
* �����C���x�֐�	-60 <= T <= 200
* ���́F����[kPa]
*       ��G���^���s�[[kJ/kg(DA)]
*       ���Ύ��x[-]
* �o�́F���x[deg.C]
************************************************/
double Air_Phphit ( double P, double h, double phi )
{
	double T;
	double h_;

	int i, j;
	double X[2], Y[2];

	// ����l
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

			// ��G���^���s�[
			h_ = Air_Ptphih ( P, T, phi );

			X[j] = T;
			Y[j] = h - h_;
		} // end of j loop
		
		// ��������
		if ( fabs ( Y[1] ) < Error ) break;
		else T = X[1] - Y[1] * ( X[0] - X[1] ) / ( Y[0] - Y[1] ) * 0.8;

		// �J��Ԃ��v�Z��
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
* �O�a�����C���֐�	-60 <= T <= 365
* ���́F���x[deg.C]����[kPa]
* �o�́F����[kPa]
************************************************/
double Air_tPs ( double T )
{
	double Pvs;//M.CONDE ENGINEERING, 2007��6�y�[�W�ڂ�Mole Fraction and Humidity Ratio��� Psv:�O�a�����C��
	double t, s;
	double Pcr, Tcr;
	double A, A1, A2, A3, A4, A5, A6;

	Pcr = 22064.0; // kPa M.CONDE ENGINEERING, 2007��6�y�[�W�ڂ�Mole Fraction and Humidity Ratio�̎���� Pcr,H2O
	Tcr = 647.14; // K  M.CONDE ENGINEERING, 2007��6�y�[�W�ڂ�Mole Fraction and Humidity Ratio�̎���� Tcr,H2O

	// �P�ʊ��Z
	T += 273.15;

	// �ƁC��
	s = T / Tcr;//��=T/Tcr,H2O
	t = 1.0 - s;//��=1-��

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

	// �O�a�����C��
	A = A1 * t + A2 * pow ( t, 1.5 ) + A3 * pow ( t, 3.0 ) + A4 * pow ( t, 3.5 ) + A5 * pow ( t, 4.0 ) + A6 * pow ( t, 7.5 );//M.CONDE ENGINEERING, 2007��6�y�[�W�ڂ̎�[10]���ln(Psv/Pcr,H2O)=1/��*(A1*��+A2*��^1.5+A3*��^3+A4*��^3.5+A5*��^4+A6*��^7.5)  
	Pvs = Pcr * exp ( A / s );

	return ( Pvs );
}

/***********************************************
* �����C�舳��M�֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F�舳��M[kJ/(kg(DA)*K)]
************************************************/
double Air_Ptxcp ( double P, double T, double x )
{
	double cp;
	double h1, h2;
	const double dT = 0.05;

	// h1, h2
	h1 = Air_Ptxh ( P, T - dT, x );
	h2 = Air_Ptxh ( P, T + dT, x );

	// �舳��M
	cp = ( h2 - h1 ) / ( 2 * dT );

	return ( cp );
}

/***********************************************
* �����C�M�`�����֐�	0.01 <= T <= 350
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F�M�`����[kW/(m*K)]
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

	// ����
	na = 1.0 / Ma;
	nw = x / Mw;

	// ��������
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// f
	f = f_ ( P, T );

	// ����
	Pw = P * Xw / f;
	Pa = P - Pw;

	// ������C�S���W���C�M�`�B��
	visc_a = Dair_Ptvisc ( Pa, T );
	thc_a = Dair_Ptthc ( Pa, T );

	if ( x == 0.0 )	thc = thc_a;
	else
	{
		// �����C�S���W���C�M�`�B��
		visc_w = sh_myuv ( Pw, T );
		thc_w = sh_lamv ( Pw, T );

		// Gaw, Gwa
		Gaw = sqrt ( 2.0 ) / 4.0 * sqrt ( Mw / ( Ma + Mw ) ) * pow ( 1.0 + sqrt ( visc_a / visc_w ) * pow ( Mw / Ma, 0.25 ), 2.0 );
		Gwa = sqrt ( 2.0 ) / 4.0 * sqrt ( Ma / ( Ma + Mw ) ) * pow ( 1.0 + sqrt ( visc_w / visc_a ) * pow ( Ma / Mw, 0.25 ), 2.0 );

		// �M�`����
		thc = thc_a / ( 1.0 + Gaw * Xw / Xa ) + thc_w / ( 1.0 + Gwa * Xa / Xw );
	}

	// �P�ʊ��Z
	thc *= 1e-3; // W/(m*K) -> kW/(m*K)

	return ( thc );
}

/***********************************************
* �����C�S���W���֐�	0.01 <= T <= 350
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F�S���W��[Pa*s]
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

	// ����
	na = 1.0 / Ma;
	nw = x / Mw;

	// ��������
	Xa = na / ( na + nw );
	Xw = nw / ( na + nw );

	// f
	f = f_ ( P, T );

	// ����
	Pw = P * Xw / f;
	Pa = P - Pw;

	// ������C�S���W��
	visc_a = Dair_Ptvisc ( Pa, T );

	if ( x == 0.0 )	visc = visc_a;
	else
	{
		// �����C�S���W���C�M�`�B��
		visc_w = sh_myuv ( Pw, T );

		// Gaw, Gwa
		Gaw = sqrt ( 2.0 ) / 4.0 * sqrt ( Mw / ( Ma + Mw ) ) * pow ( 1.0 + sqrt ( visc_a / visc_w ) * pow ( Mw / Ma, 0.25 ), 2.0 );
		Gwa = sqrt ( 2.0 ) / 4.0 * sqrt ( Ma / ( Ma + Mw ) ) * pow ( 1.0 + sqrt ( visc_w / visc_a ) * pow ( Ma / Mw, 0.25 ), 2.0 );

		// �S���W��
		visc = visc_a / ( 1.0 + Gaw * Xw / Xa ) + visc_w / ( 1.0 + Gwa * Xa / Xw );
	}

	return ( visc );
}

/***********************************************
* �����C���S���W���֐�	0.01 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F���S���W��[m^2/s]
************************************************/
double Air_Ptxkine ( double P, double T, double x )
{
	double kine;
	double visc, rho;

	// ���x
	rho = Air_Ptxrho ( P, T, x );

	// �S���W��
	visc = Air_Ptxvisc ( P, T, x );

	// ���S���W��
	kine = visc / rho;

	return ( kine );
}

/***********************************************
* �����C�M�g�U�W���֐�	0.01 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F�M�g�U�W��[m^2/s]
************************************************/
double Air_Ptxthd ( double P, double T, double x )
{
	double thd;
	double rho, cp, thc;

	// ���x
	rho = Air_Ptxrho ( P, T, x );

	// ��M
	cp = Air_Ptxcp ( P, T, x );

	// �M�`����
	thc = Air_Ptxthc ( P, T, x );

	thd = thc / ( rho * cp );

	return ( thd );
}

/***********************************************
* �����C�����C�̊g�U�W���֐�	-20 <= T <= 300
* ���́F����[kPa]
*       ���x[deg.C]
* �o�́F�����C�̊g�U�W��[m^2/s]
************************************************/
double Air_PtDh2o ( double P, double T )
{
	double D;

	// �P�ʊ��Z
	T += 273.15; // deg.C -> K
	P *= 1e3; // kPa -> Pa

	// �ꍇ����
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
* �����C�v�����g�����֐�	0.01 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F�v�����g����[-]
************************************************/
double Air_PtxPr ( double P, double T, double x )
{
	double Pr;
	double thd, kine;

	// �M�g�U�W��
	thd = Air_Ptxthd ( P, T, x );

	// ���S���W��
	kine = Air_Ptxkine ( P, T, x );

	// �v�����g����
	Pr = kine / thd;

	return ( Pr );
}

/***********************************************
* �����C�V���~�b�g���֐�	0.01 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F�V���~�b�g��[-]
************************************************/
double Air_PtxSc ( double P, double T, double x )
{
	double Sc;
	double D, kine;

	// �g�U�W��
	D = Air_PtDh2o ( P, T );

	// ���S���W��
	kine = Air_Ptxkine ( P, T, x );

	// �V���~�b�g��
	Sc = kine / D;

	return ( Sc );
}

/***********************************************
* �����C�I�_���x�֐�	x >= 0.0
* ���́F����[kPa]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F���x[deg.C]
************************************************/
double Air_Pxtdp ( double P, double x )
{
	double x_;
	double Tdp;
	double TL, TR, dT;

	// �������ݕ�
	dT = 1.0;

	// �I�_���x����l
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

	// �I�_���x
	if(x==0.0)	Tdp = -100.0;
	else
	{
		// ����I�_��Ύ��x
		x_ = Air_Ptphix ( P, Tdp, 1.0 );

		LoopCount = 0;
		while ( x_ >= x || LoopCount < Loop )
		{
			// �I�_���x����
			Tdp -= dT;

			// ����I�_��Ύ��x
			x_ = Air_Ptphix ( P, Tdp, 1.0 );

			LoopCount++;
		}

		TL = Tdp;
		TR = Tdp + dT;

		LoopCount = 0;
		while ( LoopCount < Loop )
		{
			// �I�_���x����
			Tdp = ( TL + TR ) / 2.0;

			// ����I�_��Ύ��x
			x_ = Air_Ptphix ( P, Tdp, 1.0 );

			// ��������
			if ( fabs ( 1.0 - x_ / x ) < Error ) break;

			// �l�̍X�V
			if ( x_ > x ) TR = Tdp;
			else TL = Tdp;

			LoopCount++;
		}
	}

	return ( Tdp );
}

/***********************************************
* �����C�Ìœ_���x�֐�	P <= 200.0[MPa]
* ���́F����[kPa]
* �o�́F���x[deg.C]
************************************************/
double Air_Ptmel ( double P )
{
	double Tm;
	double Ptp, Ttp;

	Ptp = 0.006112; // bar
	Ttp = 273.16; // K

	// �P�ʊ��Z
	P /= 1e8; // kPa -> bar

	// �Ìœ_���x
	Tm = Ttp * pow ( - ( P - Ptp ) / 3950 + 1.0, 1.0 / 9.0 );

	// �P�ʊ��Z
	Tm -= 273.15; // K -> deg.C

	return ( Tm );
}



// ���L�̊֐��͏�L�̊֐����Ŏg�p����֐�

/***********************************************
* �����C�������x�֐�	-60 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ��Ύ��x[kg/kg(DA)]
* �o�́F���x[deg.C]
************************************************/
double Air_Ptxtwb ( double P, double T, double x )
{
	double Twb;
	double phi;

	// ���Ύ��x
	phi = Air_Ptxphi ( P, T, x );

	// �������x
	Twb = Air_Ptphitwb ( P, T, phi );

	return ( Twb );
}

/***********************************************
* �����C�������x�֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ���Ύ��x[-]
* �o�́F���x[deg.C]
************************************************/
double Air_Ptphitwb ( double P, double T, double phi )
{
	double Twb;
	double x, xwb;
	double h, hwb, hc, h_;

	double TL, TR, dT;

	// �������ݕ�
	dT = 1.0;

	// �������x
	if ( phi == 1.0 )	Twb = T;
	else
	{
		// �������x����
		Twb = T;

		// ��Ύ��x
		x = Air_Ptphix ( P, T, phi );

		// ��G���^���s�[
		h = Air_Ptxh ( P, T, x );

		// �O�a��Ύ��x
		xwb = Air_Ptphix ( P, Twb, 1.0 );

		// �O�a��G���^���s�[
		hwb = Air_Ptxh ( P, Twb, xwb );

		// �Ïk����G���^���s�[
		hc = Air_Pthc ( P, Twb );

		LoopCount = 0;
		while ( h < hwb - ( xwb - x ) * hc  )
		{
			// �������x����
			Twb -= dT;

			// �O�a��Ύ��x
			xwb = Air_Ptphix ( P, Twb, 1.0 );

			// �O�a��G���^���s�[
			hwb = Air_Ptxh ( P, Twb, xwb );

			// �Ïk����G���^���s�[
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
			// �������x����
			Twb = ( TL + TR ) / 2.0;

			// �O�a��Ύ��x
			xwb = Air_Ptphix ( P, Twb, 1.0 );

			// �O�a��G���^���s�[
			hwb = Air_Ptxh ( P, Twb, xwb );

			// �Ïk����G���^���s�[
			hc = Air_Pthc ( P, Twb );

			// ��G���^���s�[
			h_ = hwb - ( xwb - x ) * hc;

			// ��������
			if ( fabs ( 1.0 - h_ / h ) < Error && fabs ( 1.0 - TL / Twb ) < Error ) break;

			// �l�̍X�V
			if ( h_ > h ) TR = Twb;
			else TL = Twb;

			LoopCount++;
		}

	}

	return ( Twb );
}

/***********************************************
* �����C�������x�֐�	-100 <= T <= 200
* ���́F����[kPa]
*       ���x[deg.C]
*       ���Ύ��x[-]
* �o�́F���x[deg.C]
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

	// ������
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

	// �O�a�����C��
	Psv = Air_tPs ( T );

	// �Ìœ_���x
	Tm = Air_Ptmel ( P );

	// �Ïk���i�ő́j
	if ( T < Tm )
	{
		// �P�ʊ��Z
		T += 273.15;

		for ( i = 0 ; i <= 3 ; i++ )	zT += z[i] * pow ( T, i );

		// ��G���^���s�[
		h = zT + z[4] * Psv;
	}
	// �Ïk���i�t�́j
	else
	{
		// �O�a�����C��
		T_delta = T * Delta;
		Psv_delta = Air_tPs ( T_delta );

		// �P�ʊ��Z
		T += 273.15; // deg.C -> K

		// �ꍇ����
		if ( T <= 373.15 )
		{
			for ( i = 0 ; i <= 4 ; i++ ) lT += l[i] * pow ( T, i );

			// ��
			a = lT + l[5] * pow ( 10.0, l[6] * ( T - 273.15 ) );
		}
		else if ( T > 373.15 && T <= 403.128 )
		{
			for ( i = 0 ; i <= 4 ; i++ ) oT += o[i] * pow ( T, i );

			// ��
			a = oT;
		}
		else if ( T > 403.128 && T <= 473.15 )
		{
			for ( i = 0 ; i <= 4 ; i++ ) oT += o[i] * pow ( T, i );

			// ��
			a = oT - o[5] * pow ( T - 403.128, 3.1 );
		}

		// ��e��
		for ( i = 0 ; i <= 5 ; i++ )	gT1 += g[i] * pow ( T, i );
		for ( i = 6 ; i <= 7 ; i++ )	gT2 += g[i] * pow ( T, i - 6 );
		v = 18015.28 * gT2 / gT1;

		// �P�ʊ��Z
		v *= 1e-6; // cm^3/kg -> m^3/kg

		// ��
		b = T * v * ( Psv_delta - Psv ) / ( T_delta - ( T - 273.15 ) );

		// ��G���^���s�[
		h = a - 0.01214 + b;
	}
	
	return ( h );
}

/***********************************************
* Bm, Cm�v�Z�p�֐�
* ���́F���x[K]
*       ������[-]
* �o�́FBm[cm^3/mol]
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

	// ������
	Baa = Baw = Bww = 0.0;
	Caaa = Caaw = Caww = Cwww = 0.0;

	//EOS for Air ���݋C�̂ɕ␳���邽�߁C�r���A���������@�@Z=PVa/(RT)=1+Baa/Va+Caaa/(Va^2)+�c�c�ɍ����悤�ɂ���DVw:1mol������̑̐�(mol�̐�)�D�r���A���W���͑�O�r���A���W���ȉ��͏ȗ����Ă���D
	//B1[i]:M.CONDE ENGINEERING, 2007��2�y�[�W�̎�[1]�ɓ����W��  Baa:[cm^3/mol]=��(i=0�`3)B1[i]*T^(-i)
	B1[0] = 0.349568 * 1e2;
	B1[1] = -0.668772 * 1e4;
	B1[2] = -0.210141 * 1e7;
	B1[3] = 0.924746 * 1e8;
	//C1[i]:M.CONDE ENGINEERING, 2007��2�y�[�W�̎�[1]�ɓ����W�� Caaa:[cm^6/(mol^2)]=��(i=0�`3)C1[i]*T^(-i)  C1[3]��M.CONDE ENGINEERING, 2007��2�y�[�W�̕\�ɂ͏�����Ă��Ȃ��c�c
	C1[0] = 0.125975 * 1e4;
	C1[1] = -0.190905 * 1e6;
	C1[2] = 0.632467 * 1e8;

	//EOS for Water Vapour ���݋C�̂ɕ␳���邽�߁C�r���A���������@�@Z=PVw/(RT)=1+Bww/Vw+Cwww/(Vw^2)+�c�c�ɍ����悤�ɂ���DVw:1mol������̑̐�(mol�̐�)�D�r���A���W���͑�O�r���A���W���ȉ��͏ȗ����Ă���D
	//B2[i]:M.CONDE ENGINEERING, 2007��2�y�[�W�̎�[2]�ɓ����W��  Bww:[cm^3/mol]=RT*(B2[0]+B2[1]*e^(B2[2]/T))
	B2[0] = 0.70 * 1e-8;
	B2[1] = -0.147184 * 1e-8;
	B2[2] = 1734.29;
	//C2[i]:M.CONDE ENGINEERING, 2007��2�y�[�W�̎�[2]�ɓ����W��  Cwww:[cm^6/(mol^2)]=((RT)^2)*((C2[0]+C2[1]*e^(C2[2]/T))^2)*Bww^2
	C2[0] = 0.104 * 1e-14;
	C2[1] = -0.335297 * 1e-17;
	C2[2] = 3645.09;

	//EOS for the Mixture(humid air) ���݋C�̂ɕ␳���邽�߁C�r���A���������@�@Z=P*����Vm/(RT)=1+Bm/����Vm+Cm/(����Vm^2)+�c�c�ɍ����悤�ɂ���D����Vm:1mol������̑̐�(mol�̐�)�D�r���A���W���͑�O�r���A���W���ȉ��͏ȗ����Ă���D
	//Bm[cm^3/mol]=(Xa^2)*Baa+2*Xa*Xw*Baw+(Xw^2)*Bww  M.CONDE ENGINEERING, 2007��3�y�[�W�̎�[3]
	//Baw[cm^3/mol]=��(i=0�`4)D[i]*(T^(-i))  M.CONDE ENGINEERING, 2007��3�y�[�W�̎�[5]
	//D[i]:M.CONDE ENGINEERING, 2007��2�y�[�W�̎�[2]�ɓ����W��  Bww:[cm^3/mol]=RT*(B2[0]+B2[1]*e^(B2[2]/T))
	D[0] = 0.32366097 * 1e2;
	D[1] = -0.141138 * 1e5;
	D[2] = -0.1244535 * 1e7;
	D[3] = 0.0;
	D[4] = -0.2348789 * 1e10;

	//Cm[cm^6/(mol^2)]=Xa^3*Caaa+3*(Xa^2)*Xw*Caaw+3*Xa*(Xw^2)*Caww+(Xw^3)*Cwww  M.CONDE ENGINEERING, 2007��3�y�[�W�̎�[3]
	//Caaw[(cm^6)/(mol^2)]=��(i=0�`4)E[i]*(T^-i)
	E[0] = 0.482737 * 1e3;
	E[1] = 0.105678 * 1e6;
	E[2] = -0.656394 * 1e8;
	E[3] = 0.299444 * 1e11;
	E[4] = -0.319317 * 1e13;
	//Caww[(cm^6)/(mol^2)]=F4*��(i=0�`4)E[i]*(T^-i)
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
* �C���W��f�Z�o�֐�	 M.CONDE ENGINEERING, 2007��7�y�[�W��[12]
* ���́F����[kPa]
*       ���x[deg.C]
* �o�́F�C���W��[-]
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

	// ������
	gT1 = gT2 = 0.0;
	sT = 0.0;
	xSo = 0.0;
	xSn = 0.0;
	Baa = Baw = Bww = 0.0;
	Caaa = Caaw = Caww = Cwww = 0.0;
	//g[i]:liquid water �Ïk�ߒ��̉t���̔�e��[cm^3/mol]�̌v�Z�ɗp����W�� M.CONDE ENGINEERING, 2007��7�y�[�W�ڂ̉��̕\���
	g[0] = -0.2403360201 * 1e4;
	g[1] = -0.140758895 * 1e1;
	g[2] = 0.1068287657;
	g[3] = -0.2914492351 * 1e-3;
	g[4] = 0.373497936 * 1e-6;
	g[5] = -0.21203787 * 1e-9;
	g[6] = -0.3424442728 * 1e1;
	g[7] = 0.1619785 * 1e-1;
	//s[i]:solid water �Ïk�ߒ��̌ő��̔�e��[cm^3/mol]�̌v�Z�ɗp����W�� M.CONDE ENGINEERING, 2007��7�y�[�W�ڂ̉��̕\���
	s[0] = 0.1070003 * 1e-2;
	s[1] = -0.249936 * 1e-7;
	s[2] = 0.371611 * 1e-9;
	//Ka�̒l�ɂ͋�C�̎�v�ȍ\���v�f�ł���_�f�ƒ��f��K��p����D�_�f��Ko�ƒ��f��Kn�̌v�Z�ɗp�����̕����l�̓��� Ka=Ko*Kn/(0.22*Ko+0.78*Kn)*101325.16  M.CONDE ENGINEERING, 2007��9�y�[�W�ڂ̎�[18]
	Kmax_o = 7.08 * 1e4;//M.CONDE ENGINEERING, 2007��9�y�[�W�ڂ̕\���Kmax[atm/mol fraction](oxygen)=7.08*10^4
	Kmax_n = 12.39 * 1e4;//M.CONDE ENGINEERING, 2007��9�y�[�W�ڂ̕\���Kmax[atm/mol fraction](nitrogen)=12.39*10^4
	//x[i]: Himmelblau����͂��\�z�������ɗp����W�� K=Kmax*10^(��(i=0�`4)x[i]*��^(-i) )�@M.CONDE ENGINEERING, 2007��9�y�[�W�ڂ̎�[19]
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

	// �O�a���C��
	Pvs = Air_tPs ( T );

	// �Ìœ_���x
	Tm = Air_Ptmel ( P );

	// �P�ʊ��Z
	P *= 1e3; // kPa -> Pa
	Pvs *= 1e3; // kPa -> Pa
	T += 273.15; // deg.C -> K
	Tm += 273.15; // deg.C -> K
	R0 = R * 1e6;

	// �Ïk���i�t�́jM.CONDE ENGINEERING, 2007��7�y�[�W��[13]
	if ( T > Tm )
	{
		for ( i = 0 ; i <= 5 ; i++ )	gT1 += g[i] * pow ( T, i );
		for ( i = 6 ; i <= 7 ; i++ )	gT2 += g[i] * pow ( T, i - 6 );

		// ��e��[cm^3/mol]
		v = 18015.28 * gT2 / gT1;
	}
	// �Ïk���i�ő́jM.CONDE ENGINEERING, 2007��7�y�[�W��[14]
	else
	{
		for ( i = 0 ; i <= 2 ; i++ )	sT += s[i] * pow ( T, i );

		// ��e��[cm^3/mol]
		v = 18015.28 * sT;
	}


	// 0.0 deg.C�ȏ�
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
	// 0.0 deg.C�ȉ�
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
	// ����
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

		// ��������
		if ( fabs ( 1 - f_ / f ) < Error ) break;

		// �l�̍X�V
		f = f_;

		LoopCount++;
	}
	return ( f_ );
}

/***********************************************
* ������C�M�`�����֐�		193 <= T <= 1727
* ���́F����[kPa]
*       ���x[deg.C]
* �o�́F�M�`����[kW/(m*K)]
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

	// ������
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

	// ���x
	rho = Air_Ptxrho ( P, T, 0.0 );

	// Tr, rhor
	Tr = ( T + 273.15 ) / Tc;
	rhor = rho / rhoc;

	// thc0
	thc0 = C1 * Tr + C05 * pow ( Tr, 0.5 );
	for ( i = 0 ; i <= 4 ; i++ )	thc0 += C[i] * pow ( Tr, -i );

	// thc1
	for ( i = 0 ; i <= 4 ; i++ )	thc1 += D[i] * pow ( rhor, i + 1 );

	// �M�`����
	thc = l * ( thc0 + thc1 );

	return ( thc );
}

/***********************************************
* ������C�S���W���֐�		193 <= T <= 1727
* ���́F����[kPa]
*       ���x[deg.C]
* �o�́F�S���W��[Pa*s]
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

	// ������
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

	// ���x
	rho = Air_Ptxrho ( P, T, 0.0 );

	// Tr, rhor
	Tr = ( T + 273.15 ) / Tc;
	rhor = rho / rhoc;

	// visc0
	visc0 = A1 * Tr + A05 * pow ( Tr, 0.5 );
	for ( i = 0 ; i <= 4 ; i++ )	visc0 += A[i] * pow ( Tr, -i );

	// visc1
	for ( i = 0 ; i <= 3 ; i++ )	visc1 += B[i] * pow ( rhor, i + 1 );

	// �M�`����
	visc = H * ( visc0 + visc1 );

	return ( visc );
}
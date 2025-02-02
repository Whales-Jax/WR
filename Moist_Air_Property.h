/***********************************************************
 * Air_Property.h
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
#include <cmath>

#ifndef __MOIST_AIR_PROPERTY_H__
#define __MOIST_AIR_PROPERTY_H__

/*********************************************************************************
 * Air �N���X -- �C�ӂ̏�Ԃ̕������i�[����N���X
 *
 * �����o�ϐ�
 *     T       -- ���x deg.C
 *     P       -- ���� kPa
 *     rho     -- ���x kg/m^3
 *     h       -- ��G���^���s kJ/kg(DA)
 *     s       -- ��G���g���s kJ/(kg(DA)*K)
 *     x       -- ��Ύ��x kg/kg(DA)
 *	   phi	   -- ���Ύ��x -
 *     cp      -- �舳��M kJ/(kg(DA)*K)
 *     thc     -- �M�`���� kW/(m*K)
 *     visc    -- �S�� Pa*s
 *     Psat    -- �O�a�����C���� kPa
 *     G       -- ���� kg(DA)/s
 *     dP      -- ���� kPa
 * �����o�֐�
 *     Air     -- �R���X�g���N�^
 *     Air_Ptx -- ��Ԋ֐��i���́E���x�E��Ύ��x���́j
 *     Air_Phx -- ��Ԋ֐��i���́E��G���^���s�[�E��Ύ��x���́j
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

		// �R���X�g���N�^
		Air()
		{
			T	 = 0.0;
			P	 = 101.325; // ��C��
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

		// ��Ԋ֐��i���́E���x�E��Ύ��x���́j
		void Air_Ptx ( double pressure, double temperature, double humidity );
		
		// ��Ԋ֐��i���́E��G���^���s�[�E��Ύ��x���́j
		void Air_Phx ( double pressure, double enthalpy, double humidity );

		// ���L�ȈՌv�Z�p
		// ��Ԋ֐��i���x�E��Ύ��x���́j
		void Air_tx ( double temperature, double humidity );

		// ��Ԋ֐��i��G���^���s�[�E��Ύ��x���́j
		void Air_hx ( double enthalpy, double humidity );
};

// �ȉ��C�����C�����֐�

// �����C��e�ϊ֐�	-100 <= T <= 200
double Air_Ptxv ( double P, double T, double x );
// �����C��e�ϊ֐�	-60 <= T <= 200
double Air_Ptphiv ( double P, double T, double phi );

// �����C���x�֐�	-100 <= T <= 200
double Air_Ptxrho ( double P, double T, double x );
// �����C���x�֐�	-60 <= T <= 200
double Air_Ptphirho ( double P, double T, double phi );

// �����C��G���^���s�[�֐�	-100 <= T <= 200
double Air_Ptxh ( double P, double T, double x );
// �����C��G���^���s�[�֐�	-60 <= T <= 200
double Air_Ptphih ( double P, double T, double phi );

// �����C��G���g���s�[�֐�	-100 <= T <= 200
double Air_Ptxs ( double P, double T, double x );
// �����C��G���g���s�[�֐�	-60 <= T <= 200
double Air_Ptphis ( double P, double T, double phi );

// �����C���Ύ��x�֐�	-60 <= T <= 365
double Air_Ptxphi ( double P, double T, double x );
// �����C���Ύ��x�֐�	-100 <= T <= 200
double Air_Pttwbphi ( double P, double T, double Twb );

// �����C��Ύ��x�֐�	-60 <= T <= 365
double Air_Ptphix ( double P, double T, double phi );
// �����C��Ύ��x�֐�	-100 <= T <= 200
double Air_Phtx ( double P, double h, double T );
// �����C��Ύ��x�֐�	-60 <= T <= 200
double Air_Pttwbx ( double P, double T, double Twb );

// �����C���x�֐�	-100 <= T <= 200
double Air_Phxt ( double P, double h, double x );
// �����C���x�֐�	-60 <= T <= 200
double Air_Phphit ( double P, double h, double phi );

// �O�a�����C���֐�	-60 <= T <= 365
double Air_tPs ( double T );

// �����C�舳��M�֐�	-100 <= T <= 200
double Air_Ptxcp ( double P, double T, double x );

// �����C�M�`�����֐�	0.01 <= T <= 350
double Air_Ptxthc ( double P, double T, double x );

// �����C�S���W���֐�	0.01 <= T <= 350
double Air_Ptxvisc ( double P, double T, double x );

// �����C���S���W���֐�	0.01 <= T <= 200
double Air_Ptxkine ( double P, double T, double x );

// �����C�M�g�U���֐�	0.01 <= T <= 200
double Air_Ptxthd ( double P, double T, double x );

// �����C�����C�̊g�U�W���֐�	-20 <= T <= 300
double Air_PtDh2o ( double P, double T );

// �v�����g�����֐�	0.01 <= T <= 200
double Air_PtxPr ( double P, double T, double x );

// �V���~�b�g���֐�	0.01 <= T <= 200
double Air_PtxSc ( double P, double T, double x );

// �I�_���x�֐�	x >= 0.0
double Air_Pxtdp ( double P, double x );

// �Ìœ_���x�֐�	P <= 200.0[MPa]
double Air_Ptmel ( double P );



// ���L�̊֐��͏�L�̊֐����Ŏg�p����֐�

// �������x�֐�	-60 <= T <= 200
double Air_Ptxtwb ( double P, double T, double x );
// �������x�֐�	-100 <= T <= 200
double Air_Ptphitwb ( double P, double T, double phi );

// �O�a�Ïk����G���^���s�[�֐�	-100 <= T <= 200
double Air_Pthc ( double P, double T );

// Bm, Cm�v�Z�p�֐�
void Bm_Cm ( double T, double Xa, double Xw, double &Bm, double &Cm );

// �C���W��f�Z�o�֐�
double f_ ( double P, double T );

// ������C�M�`�����֐�	-193 <= T <= 1727
double Dair_Ptthc ( double P, double T );

// ������C�S���W���֐�	-193 <= T <= 1727
double Dair_Ptvisc ( double P, double T );



// �ȉ��C��C���i101.325kPa�j�ł̊ȈՋ�C�֐�

//��̐ϊ֐�
inline double Air_txv ( double T, double x )
{
	return 4.554 * 1e-3 * ( 273.15 + T ) * ( 0.622 + x );
}

// ���x�֐�
inline double Air_txrho ( double T, double x )
{
	return 1.0 / Air_txv ( T, x );
}

// ��G���^���s�[�֐�
inline double Air_txh ( double T, double x )
{
	return ( 1005.22 + 0.02615 * T ) * T + x * ( 2500800 + 1868 * T );
}

// ��C���x�֐�
inline double Air_hxt ( double h, double x )
{
	double a,b,c,d;
	double T;

	a = 1005.22;
	b = 0.02615;
	c = 2500800.0;
	d = 1868.0;

	// ���x
	T = ( -( d * x + a ) + pow( pow( x * d + a, 2.0 ) - 4.0 * b * ( x * c - h ), 0.5 ) ) / ( 2.0 * b );

	return ( T );
}

// �M�`����
inline double Air_tthc ( double T )
{
	return 0.02624 + 7.58e-5 * ( T - 27.0 );
}

// �S���W��
inline double Air_tvisc ( double T )
{
	return 1.859e-5 + 4.32e-8 * ( T - 27.0 );
}


// �f�V�J���g�v�Z�p
// �����C�̏������M
inline double Wat_tlh(double t)//latentheat
{
	return 2500800 + 1868 * t;//2358500.0 - 2460.0 * ( t - 60.0 );
}

#endif
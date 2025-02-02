/********************************************************************************
 * DXProperty_ver04.h                                                           *
 *                                                                              *
 * 2008/04 �J���`�[��                                                           *
 * Property�ɂQ�����p�̃Z�b�g�A�b�v�֐��ǉ��Dby����                             *
 * Property�ɕ������p�̃Z�b�g�A�b�v�֐��ǉ��Dver.3.0 by����						*
 ********************************************************************************/

#ifndef __DXPROPERTY_VER04_H__
#define __DXPROPERTY_VER04_H__

#include <windows.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include "CCSVFile.h"
using namespace std;

/*********************************************************************************
 * Fluid �N���X -- �C�ӂ̏�Ԃ̕������i�[����N���X
 *
 * �����o�ϐ�
 *     T       -- ���x deg.C
 *     P       -- ���� kPa
 *     rho     -- ���x kg/m^3
 *     h       -- ��G���^���s kJ/kg
 *     s       -- ��G���g���s kJ/(kg*K)
 *     x       -- �����x -
 *     cp      -- �舳��M kJ/(kg*K)
 *     thc     -- �M�`���� W/(m*K)
 *     visc    -- �S�� uPa*s
 *     st      -- �\�ʒ��� N/m
 * �����o�ϐ�
 *     Fluid   -- �R���X�g���N�^
 *********************************************************************************/
class DXFluid
{
	public:
		double T;
		double P;
		double rho;
		double h;
		double s;
		double x;
		double cp;
		double cv;
		double thc;
		double visc;
		double st;
		double u;

		// �R���X�g���N�^
		DXFluid()
		{
			T    = 0.0;
			P    = 0.0;
			rho  = 0.0;
			h    = 0.0;
			s    = 0.0;
			x    = 0.0;
			cp   = 0.0;
			cv   = 0.0;
			thc  = 0.0;
			visc = 0.0;
			st   = 0.0;
		}

		~DXFluid(){}
};
/*********************************************************************************
 * Property �N���X -- refprop�𗘗p���ė�}�̕����l���v�Z����N���X
 *
 * �����o�֐�
 *     DllSetup         -- �R���X�g���N�^
 *     ~DllSetup        -- �f�X�g���N�^
 *     freeDll          -- DLL����֐�
 *     setup            -- �Z�b�g�A�b�v�֐�
 *     setref           -- ���Ԑݒ�֐�
 *     xmole            -- ���������E���q�ʊ֐�
 *     xmass            -- ���ʕ����E���q�ʊ֐�
 *     criticalPoint    -- �ՊE��Ԋ֐�
 *     state_tp         -- �P����Ԋ֐�(���x�E���͓���)
 *     state_ph         -- �C�ӏ�Ԋ֐�(���́E��G���^���s����)
 *     sat_t            -- �O�a��Ԋ֐�(���x����)
 *     sat_p            -- �O�a��Ԋ֐�(���͓���)
 *********************************************************************************/
// DLL�Ăяo���֐��̌^�錾
typedef void (__stdcall *SetupDllType) (long&,char*,char*,char*,long&,char*,long,long,long,long);
typedef void (__stdcall *SetrefDllType) (char*,long&,double*,double&,double&,double&,double&,long&,char*,long,long);
typedef void (__stdcall *XmassDllType) (double*,double*,double&);
typedef void (__stdcall *XmoleDllType) (double*,double*,double&);
typedef void (__stdcall *EntDllType) (double&,double&,double*,double&);
typedef void (__stdcall *CritpDllType) (double*,double&,double&,double&,long&,char*,long);
typedef void (__stdcall *SatDllType) (double&,double*,long&,double&,double&,double&,double*,double*,long&,char*,long);
typedef void (__stdcall *PhflshDllType) (double&,double&,double*,double&,double&,double&,double&,double*,double*,double&,double&,double&,double&,double&,double&,long&,char*,long);
typedef void (__stdcall *PtflshDllType) (double&,double&,double*,double&,double&,double&,double*,double*,double&,double&,double&,double&,double&,double&,double&,long&,char*,long);
typedef void (__stdcall *SurtenDllType) (double&,double&,double&,double*,double*,double&,long&,char*,long);
typedef void (__stdcall *SurftDllType)(double &,double &,double *,double &,long &,char*,long);
typedef void (__stdcall *ThermdllTYPE)  (double&,double&,double*,double&,double&,double&,double&,double&,double&,double&,double&);
typedef void (__stdcall *TpflshDllType) (double&,double&,double*,double&,double&,double&,double*,double*,double&, double&,double&,double&,double&,double&,double&,long&,char*,long);
typedef void (__stdcall *TprhoDllType)(double&,double&,double*,long&,long&,double&,long&,char*,long);
typedef void (__stdcall *TrnprpDllType)(double&,double&,double*,double&,double&,long&,char*,long);

class DXProperty
{
	private:
		// DLL�̃n���h��
		HINSTANCE hRefpropDll;
		// DLL���̊֐��̃G���g��
		SetupDllType  setupDll;
		SetrefDllType setrefDll;
		XmassDllType  xmassDll;
		XmoleDllType  xmoleDll;

		CritpDllType  critpDll;
		SatDllType    sattDll;
		SatDllType    satpDll;
		PhflshDllType phflshDll;
		PtflshDllType ptflshDll;
		SurtenDllType surtenDll;
		SurftDllType  surftDll;
		ThermdllTYPE  thermDll;
		TpflshDllType tpflshDll;
		TprhoDllType  tprhoDll;
		TrnprpDllType trnprpDll;

		// ������
		long numberOfComponents_;
		// ��������
		double *moleFraction;

	public:
		// ���ʕ����E���q��
		double *massFraction;
		double molecularWeight;
		ofstream refdata;

		DXFluid **data_table;	//!<state_ph�p�f�[�^�e�[�u��
		DXFluid **data_table_p;	//!<sat_p�p�f�[�^�e�[�u��
		DXFluid **data_table_tp;//!<state_tp�p�f�[�^�e�[�u��
		DXFluid **data_table_t;	//!<sat_t�p�f�[�^�e�[�u��
		DXFluid **data_table_du;//!<state_vu�p�f�[�^�e�[�u��
		int rep_h;
		int rep_p;
		int rep_t;
		int rep_d;
		int rep_u;
		bool _DLLsetup;

		int calccount;

		double min_P;
		double min_h;
		double min_T;
		double max_P;
		double max_h;
		double max_T;

		string PATH;
		ostringstream name;

		string reibai;
		string joutai;
		DXFluid Rc,Rl,Rv;


		// �R���X�g���N�^
		explicit DXProperty();
		// �f�X�g���N�^
		virtual ~DXProperty();
		// DLL����֐�
		void freeDll();
		void LoadDLL(std::string DllFileName);
		// �Z�b�g�A�b�v�֐�
		void setup();
		// ���������E���q�ʊ֐�
		void xmole(double *massFraction, double *moleFraction, double &molecularWeight);
		// ���ʕ����E���q�ʊ֐�
		void xmass(double *moleFraction, double *massFraction, double &molecularWeight);
		// ���q�ʊ֐�
		double molecularWeight_f();
		// �ՊE��Ԋ֐�
		void criticalPoint(DXFluid &crit);
		// �\�ʒ��͊֐�
		double surfaceTension_t(double temperature);

		// �P����Ԋ֐�
		void state_tp(double temperature, double pressure);
		// �C�ӏ�Ԋ֐�
		void state_ph(double pressure, double enthalpy);
		// �O�a��Ԋ֐�(���x����)
		void sat_t(double temperature);
		// �O�a��Ԋ֐�(���͓���)
		void sat_p(double pressure);




		void CopyTable( DXProperty &ref );
		void set_step( int h , int p , int t);
		void init_table(void);
		void load_table(std::string fname);
		void load_table2(std::string fname);
		void make_table(double start_h , double end_h , double start_p , double end_p ,double start_t , double end_t );
		bool state_ph2(double pressure, double enthalpy);
		bool sat_p2(double pressure);
		bool state_tp2(double temperature , double pressure );
		bool sat_t2(double temperature);

		bool state_ph3(double pressure, double enthalpy);
		bool sat_p3(double pressure);
		bool state_tp3(double temperature , double pressure );
		bool sat_t3(double temperature);


		//state_vu
		void init_table_du(void);
		void load_table_du(std::string fname);
		void make_table_du(double start_d , double end_d , double start_u , double end_u );
		bool state_du2(double density , double energy );

};

#endif
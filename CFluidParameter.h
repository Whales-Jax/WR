#ifndef __CFluidParameter_h
#define __CFluidParameter_h

#include <iostream>
#include <fstream>
#include "DXProperty_ver06.h"

/*!
���̏�ԕϐ�
*/
class CFluidParameter{
public:

	CFluidParameter();
	~CFluidParameter();

	int CycleNUM;

	double G;//!<���ʗ���[kg/s]
	double P;//!<����[kPa]
	double h;//!<��G���^���s[kJ/kg]
	double x;//!<�����x
	double X;//!<�Z�x
	double s;//!<��G���g���s[]
	double D;//!<���x[kg/m3]
	double u;//!<������G�l���M�[[kJ/kg]
	double M;//!<����[kg]
	double Q;//!<�M��[kW]
	double q;//!<�M����[kW/m2]
	double T;//!<���x[degC]
	double Tw;//!<�������x[degC]
	double dT;//!<���x��[degC]
	double v;//!<��̐�[m3/kg]
	double V;//!<�̐ϗ���[m3/s]
	double V_lpm;//!<�̐ϗ���[l/min]
	double velocity;//!<����[m/s]
	double cp;//!<�舳��M[kJ/(kg*K)]
	double thc;//!<�M�`����[W/(m*K)]
	double visc;//!<�S���W��[uPa*s]
	double st;//!<�\�ʒ���[N/m]
	double cv;
	double satT;//!�O�a���x
	double phi;//!���Ύ��x
	double nu;//!���S�x
	
	double checkT;//!<���x[degC]

	double cpL;//!<�舳��M[kJ/(kg*K)]
	double thcL;//!<�M�`����[W/(m*K)]
	double viscL;//!<�S���W��[uPa*s]
	double cpV;//!<�舳��M[kJ/(kg*K)]
	double thcV;//!<�M�`����[W/(m*K)]
	double viscV;//!<�S���W��[uPa*s]


	double alpha;//!<�M�`�B��
	double alpharate;//!<�M�`�B��
	double alphaL;//!<�M�`�B��
	double alphaV;//!<�M�`�B��
	double alphaLO;//!<�M�`�B��
	double alphaAS;//!<�M�`�B��
	double alphaA;//!<�M�`�B��
	double alphaDS;//!<�M�`�B��
	double alphaD;//!<�M�`�B��
	double pd;//!<���͑���[kPa]
	double pd_H;//!<�w�b�h���̈��͑���[kPa]
	double dp;//!<���͑���[kPa/m]

	double Gi;//!<��������[kg/s]
	double Pi;//!<��������[kPs]
	double hi;//!<������G���^���s[kJ/kg]
	double Xi;//!<�����Z�x
	double Di;//!<�������x
	double si;//!<�����G���g���s
	double Ti;//!<�������x
	double Twi;//!<�����������x
	double xi;//!<���������x
	double Vi;//!<�����̐ϗ���
	double Qi;//!<�M��[kW]
	double SHi;//!<�ߔM�x[degC]
	double SCi;//!<�T�u�N�[���x[degC]

	double Gj;//!<��������[kg/s]
	double Pj;//!<��������[kPs]
	double hj;//!<������G���^���s[kJ/kg]
	double Xj;//!<�����Z�x
	double Dj;//!<�������x
	double sj;//!<�����G���g���s
	double Tj;//!<�������x
	double Twj;//!<�����������x
	double xj;//!<���������x
	double Vj;//!<�����̐ϗ���
	double Qj;//!<�M��[kW]
	double SHj;//!<�ߔM�x[degC]
	double SCj;//!<�T�u�N�[���x[degC]

	double Go;//!<�o������[kg/s]
	double Po;//!<�o������[kPa]
	double ho;//!<�o����G���^���s[kJ/kg]
	double Xo;//!<�o���Z�x
	double Do;//!<�o�����x
	double so;//!<�o���G���g���s
	double To;//!<�o�����x
	double Two;//!<�o���������x
	double xo;//!<�o�������x
	double Vo;//!<�o���̐ϗ���
	double Qo;//!<�M��[kW]
	double SHo;//!<�ߔM�x[degC]
	double SCo;//!<�T�u�N�[���x[degC]

	double Gp;//!<�o������[kg/s]
	double Pp;//!<�o������[kPa]
	double hp;//!<�o����G���^���s[kJ/kg]
	double Xp;//!<�o���Z�x
	double Dp;//!<�o�����x
	double sp;//!<�o���G���g���s
	double Tp;//!<�o�����x
	double Twp;//!<�o���������x
	double xp;//!<�o�������x
	double Vp;//!<�o���̐ϗ���
	double Qp;//!<�M��[kW]
	double SHp;//!<�ߔM�x[degC]
	double SCp;//!<�T�u�N�[���x[degC]

	double Vi_m3ph;//!<��������[m3/h]
	double Vj_m3ph;//!<��������[m3/h]
	double Vo_m3ph;//!<�o������[m3/h]
	double Vp_m3ph;//!<�o������[m3/h]

	double Vi_m3ps;//!<��������[m3/s]
	double Vj_m3ps;//!<��������[m3/s]
	double Vo_m3ps;//!<�o������[m3/s]
	double Vp_m3ps;//!<�o������[m3/s]


	double Vi_lpm;//!<��������[l/min]
	double Vj_lpm;//!<��������[l/min]
	double Vo_lpm;//!<�o������[l/min]
	double Vp_lpm;//!<�o������[l/min]

	double PV;//!<�C�̂̈���[kPs]
	double hV;//!<�C�̂̔�G���^���s[kJ/kg]
	double XV;//!<�C�̂̔Z�x
	double DV;//!<�C�̖̂��x
	double sV;//!<�C�̂̃G���g���s
	double uV;//!<�C�̂̓����G�l���M�[
	double TV;//!<�C�̂̉��x
	double MV;//!<�C�̂̎���

	double PL;//!<�C�̂̈���[kPs]
	double hL;//!<�C�̂̔�G���^���s[kJ/kg]
	double XL;//!<�C�̂̔Z�x
	double DL;//!<�C�̖̂��x
	double sL;//!<�C�̂̃G���g���s
	double uL;//!<�C�̂̓����G�l���M�[
	double TL;//!<�C�̂̉��x
	double ML;//!<�C�̂̎���

	double Xtt;
	double gxi;
	double Hgxi;
	double phiv; 

	double ep;
	double W;

	double Re;//!<���C�m���Y��

	double Re1;//!<���C�m���Y��
	double Re2;//!<���C�m���Y��
	double Re3;//!<���C�m���Y��
	double ReL;//!<���C�m���Y��
	double ReV;//!<���C�m���Y��
	double ReLO;


	double Pr1;//�v�����g����
	double Pr2;//�v�����g����
	double Pr3;//�v�����g����
	double PrL;//�v�����g����
	double PrV;//�v�����g����

	double Bo1;//!<�{�C�����O��
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

	char IOCode;//!<0:�o�� 1:����
	int  IONum;//!<�o�����̔ԍ�
	char *FluidCode;//!<�����}�̖��O"AMMONIA.FLD"��"R134a"�Ȃ�
	int  RouteCode;//!<�o�H�̔ԍ�

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
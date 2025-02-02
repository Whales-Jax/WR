/********************************************************************************
 * DXProperty_ver05.cpp                                                         *
 *                                                                              *
 * 2008/04 �J���`�[��                                                           *
 * Property�ɂQ�����p�̃Z�b�g�A�b�v�֐��ǉ��Dby����                             *
 * Property�ɕ������p�̃Z�b�g�A�b�v�֐��ǉ��Dver.3.0 by����						*
 ********************************************************************************/
#include <iostream>
#include <string>
#include <windows.h>
const long MaxComponent = 20;
#include "DXProperty_ver06.h"

#pragma warning(disable:4996)

// �����񒷂�
const size_t LFluidFile = 10000;// ���̃t�@�C���̃p�X������
const long   LMixMoldel = 255;  // �������f���t�@�C���̃p�X������
const long   LReference = 3;    // ���Ԗ�����
const long   LError     = 255;  // �G���[���b�Z�[�W����

/*!
*******************************************************************************
* Property::Property -- �R���X�g���N�^  DLL�̃��[�h�Ɗ֐��A�h���X�̎擾
*
* ����
*     DllFileName  -- DLL�t�@�C����
*******************************************************************************
 */
DXProperty::DXProperty()
{
/*	// DLL�̃��[�h
	hRefpropDll = LoadLibrary(DllFileName.c_str());
	// ���[�h�G���[
	if(hRefpropDll == 0)
	{
		std::cout << DllFileName << "�̃��[�h�Ɏ��s" << std::endl;
		exit(1);
	}
	else
	{
		// DLL�̊֐��A�h���X�̎擾
		setupDll  = (SetupDllType)  GetProcAddress(hRefpropDll, "SETUPdll");
		setrefDll = (SetrefDllType) GetProcAddress(hRefpropDll, "SETREFdll");
		xmassDll  = (XmassDllType)  GetProcAddress(hRefpropDll, "XMASSdll");
		xmoleDll  = (XmoleDllType)  GetProcAddress(hRefpropDll, "XMOLEdll");
		critpDll  = (CritpDllType)  GetProcAddress(hRefpropDll, "CRITPdll");
		tpflshDll = (TpflshDllType) GetProcAddress(hRefpropDll, "TPFLSHdll");
		phflshDll = (PhflshDllType) GetProcAddress(hRefpropDll, "PHFLSHdll");
		sattDll   = (SatDllType)    GetProcAddress(hRefpropDll, "SATTdll");
		satpDll   = (SatDllType)    GetProcAddress(hRefpropDll, "SATPdll");
		surftDll  = (SurftDllType)  GetProcAddress(hRefpropDll, "SURFTdll");
		thermDll  = (ThermdllTYPE)  GetProcAddress(hRefpropDll, "THERMdll");
		trnprpDll = (TrnprpDllType) GetProcAddress(hRefpropDll, "TRNPRPdll");
	}*/

	// ������
	numberOfComponents_ = 0;
	molecularWeight = 0.0;
	massFraction = new double [20];
	massFraction[0] = 1.0;
	moleFraction = new double [20];
	moleFraction[0] = 1.0;

	_DLLsetup = 0;

	rep_h = 300;
	rep_p = 300;
	rep_t = 300;
	rep_d = 300;
	rep_u = 300;

	min_P = 0.0;
	min_h = 0.0;
	min_T = 0.0;
	max_P = 0.0;
	max_h = 0.0;
	max_T = 0.0;

	calccount = 0;

//	PATH = "C:\\lib\\energyflow\\";
	PATH = "";
//	PATH = "../";
//	PATH = "../../../../../lib/energyflow/";

}



/*!
*******************************************************************************
* Property::~Property -- �f�X�g���N�^  DLL�̉��
*******************************************************************************
*/
DXProperty::~DXProperty()
{
	// DLL�̊J��
	FreeLibrary(hRefpropDll);
	delete [] massFraction;
}

/*!
*******************************************************************************
* Property::freeDll -- DLL����֐�  �����I��DLL�̉��
*******************************************************************************
*/
void DXProperty::freeDll()
{
	// DLL�̊J��
	FreeLibrary(hRefpropDll);
}


/********************************************************************************
 * Property::setup -- �Z�b�g�A�b�v�֐� ���̖�����ԂȂǂ̐ݒ�
 *
 * ����(���͒l)
 *     FluidName -- ���̖�
 *     Reference -- ����
 *
 * ����
 *     1)fluids�f�B���N�g�����J�����g�f�B���N�g���ɒu���Ă���
 *     2)�����̖��O��fluids�f�B���N�g�����̃t�@�C��������I������
 *     3)���Ԃ͊��Ԑݒ�֐��̕��Ő����ɐݒ肷��
 *       �Z�b�g�A�b�v�ł̓_�~�[�Ƃ���DEF�����Ă���
 *     4)���Ԃ�NBP,ASH,IIR,DEF����I������(OTH�͑I���ł��Ȃ�)
 *       ���ꂼ��̈Ӗ��̂��Ă̓}�j���A�����Q�Ƃ��ꂽ�� 
 ********************************************************************************/
/*void DXProperty::setup(std::string FluidName, std::string Reference)
{
	////////////////////
	// ���̂̐ݒ�
	////////////////////

	// �������̐ݒ�
	numberOfComponents_ = 1;

	// ���̃t�@�C���̃p�X��
	std::string FluidFile = "fluids\\" + FluidName;

	// C������
	char cFluidFile[LFluidFile+1];// ���̃t�@�C���̃p�X��
	char cMixMoldel[LMixMoldel+1];// �������f���t�@�C���̃p�X��
	char cDummy[LReference+1];    // ����(�_�~�[)
	char cError[LError+1];        // �G���[���b�Z�[�W

	// ������̑��(string�����񂩂�z�񕶎���ւ̕ϊ�)
	std::strcpy(cFluidFile,  FluidFile.c_str());
	std::strcpy(cMixMoldel, "fluids\\HMX.BNC\0");
	std::strcpy(cDummy,     "DEF\0");//�_�~�[
	std::strcpy(cError,     "Ok\0");

	// �G���[�ԍ�
	long isetup = 0;

	// �Z�b�g�A�b�v�֐�
	setupDll(numberOfComponents_, cFluidFile, cMixMoldel, cDummy, isetup,
	         cError, LFluidFile, LMixMoldel, LReference, LError);

	////////////////////
	// ���Ԃ̐ݒ�
	////////////////////

	// ���ʕ����̐ݒ�
	massFraction[0] = 1.0;

	// ���������E���q�ʂ̐ݒ�
	xmole(massFraction, moleFraction, molecularWeight);

	// �G���[�ԍ�
	long isetref = 0;

	// �G���[���b�Z�[�W
	char cError2[LError+1];
	std::strcpy(cError2, "Ok");

	// ���ԕ�����
	char cReference[LReference+1];
	std::strcpy(cReference, Reference.c_str());

	// �_�~�[
	double dummy[4] = {0.0};
	long iDummy = 0;

	// �������t���b�O
	long ixflag = 1;

	// ���Ԃ̐ݒ�֐�
	setrefDll(cReference, ixflag, moleFraction, dummy[0], dummy[1], dummy[2], dummy[3], isetref, cError2, LError, iDummy);
}*/

/*!
*******************************************************************************
* DXProperty::LoadDLL -- DLL�̃��[�h�Ɗ֐��A�h���X�̎擾
*
* ����
*     DllFileName  -- DLL�t�@�C����
*******************************************************************************
*/
void DXProperty::LoadDLL(std::string DllFileName)
{
	// DLL�̃��[�h
	hRefpropDll = LoadLibrary( string( PATH + DllFileName ).c_str() );
	// ���[�h�G���[
	if(hRefpropDll == 0)
	{
		std::cout << DllFileName << "�̃��[�h�Ɏ��s" << std::endl;
		exit(1);
	}
	else
	{
		// DLL�̊֐��A�h���X�̎擾
		setupDll  = (SetupDllType)  GetProcAddress(hRefpropDll, "SETUPdll");
		setrefDll = (SetrefDllType) GetProcAddress(hRefpropDll, "SETREFdll");
		xmassDll  = (XmassDllType)  GetProcAddress(hRefpropDll, "XMASSdll");
		xmoleDll  = (XmoleDllType)  GetProcAddress(hRefpropDll, "XMOLEdll");
		critpDll  = (CritpDllType)  GetProcAddress(hRefpropDll, "CRITPdll");
		tpflshDll = (TpflshDllType) GetProcAddress(hRefpropDll, "TPFLSHdll");
		phflshDll = (PhflshDllType) GetProcAddress(hRefpropDll, "PHFLSHdll");
		sattDll   = (SatDllType)    GetProcAddress(hRefpropDll, "SATTdll");
		satpDll   = (SatDllType)    GetProcAddress(hRefpropDll, "SATPdll");
		surftDll  = (SurftDllType)  GetProcAddress(hRefpropDll, "SURFTdll");
		thermDll  = (ThermdllTYPE)  GetProcAddress(hRefpropDll, "THERMdll");
		trnprpDll = (TrnprpDllType) GetProcAddress(hRefpropDll, "TRNPRPdll");
		surtenDll = (SurtenDllType) GetProcAddress(hRefpropDll, "SURTENdll");

	}
}

/*!
*******************************************************************************
* DXProperty::setup -- �Z�b�g�A�b�v�֐� ���̖�����ԂȂǂ̐ݒ�
*
* ����(���͒l)
*     FluidName -- ���̖�
*     Reference -- ����
*
* ����
*     1)fluids�f�B���N�g�����J�����g�f�B���N�g���ɒu���Ă���
*     2)�����̖��O��fluids�f�B���N�g�����̃t�@�C��������I������
*     3)���Ԃ͊��Ԑݒ�֐��̕��Ő����ɐݒ肷��
*       �Z�b�g�A�b�v�ł̓_�~�[�Ƃ���DEF�����Ă���
*     4)���Ԃ�NBP,ASH,IIR,DEF����I������(OTH�͑I���ł��Ȃ�)
*       ���ꂼ��̈Ӗ��̂��Ă̓}�j���A�����Q�Ƃ��ꂽ�� 
*******************************************************************************
*/
void DXProperty::setup()
{
	////////////////////
	// ���̂̐ݒ�
	////////////////////

	std::string FluidNames = reibai;
	std::string Reference  = joutai;

	_DLLsetup = 1;

	int i, j, k;
	std::string FluidFile;

	//������
	i = 0;
	j = 0;
	k = 0;
	FluidFile = "\0";

	// ���̃t�@�C���̃p�X��
	while(1)
	{
		i++;
		k = FluidNames.find( ",", j );
		if ( k == -1 )
		{
			FluidFile += PATH + "fluids\\" + FluidNames.substr( j ) ;//+ ".fld";
			break;
		}
		FluidFile +=  PATH + "fluids\\" + FluidNames.substr( j, k - j ) /*+ ".fld"*/ + "|";
		j = k + 1;
	}

	// �������̐ݒ�
	numberOfComponents_ = i;
//	moleFraction = new double [i];

	// C������
	char cFluidFile[LFluidFile+1];// ���̃t�@�C���̃p�X��
	char cMixMoldel[LMixMoldel+1];// �������f���t�@�C���̃p�X��
	char cDummy[LReference+1];    // ����(�_�~�[)
	char cError[LError+1];        // �G���[���b�Z�[�W

	// ������̑��(string�����񂩂�z�񕶎���ւ̕ϊ�)
	strcpy_s(cFluidFile,  FluidFile.c_str());
	strcpy_s(cMixMoldel, "fluids\\HMX.BNC\0");
	strcpy_s(cDummy,     "DEF\0");//�_�~�[
	strcpy_s(cError,     "Ok\0");

	// �G���[�ԍ�
	long isetup = 0;

	// �Z�b�g�A�b�v�֐�
	setupDll(numberOfComponents_, cFluidFile, cMixMoldel, cDummy, isetup,
	         cError, LFluidFile, LMixMoldel, LReference, LError);

	////////////////////
	// ���Ԃ̐ݒ�
	////////////////////

	// ���������E���q�ʂ̐ݒ�
	xmole(massFraction, moleFraction, molecularWeight);

	// �G���[�ԍ�
	long isetref = 0;

	// �G���[���b�Z�[�W
	char cError2[LError+1];
	strcpy_s(cError2, "Ok");

	// ���ԕ�����
	char cReference[LReference+1];
	strcpy_s(cReference, Reference.c_str());

	// �_�~�[
	double dummy[4] = {0.0};
	long iDummy = 0;

	// �������t���b�O
	long ixflag = i;

	// ���Ԃ̐ݒ�֐�
	setrefDll(cReference, ixflag, moleFraction, dummy[0], dummy[1], dummy[2], dummy[3], isetref, cError2, LError, iDummy);
}

/*!
*******************************************************************************
* DXProperty::xmole -- ���������E���q�ʊ֐�  ���ʕ������烂�������ƕ��q�ʂ̎Z�o
*
* ����(���͒l)
*     massFraction    -- ���ʕ����z��
* ����(�o�͒l)
*     moleFraction    -- ���������z��
*     molecularWeight -- ���q�� g/mol
*
* ����
*     1)�P�����̏ꍇxmoleDll���烂���������Z�o�ł��Ȃ�(�����s��)
*       ���̖���������邽�߁C�P�����̏ꍇ���������ɒ���1.0�������Ă���
*******************************************************************************
*/
void DXProperty::xmole(double *massFraction, double *moleFraction, double &molecularWeight)
{
	// ���������E���q��
	xmoleDll(massFraction, moleFraction, molecularWeight);
	// �P�����̏ꍇ
	if(numberOfComponents_ == 1)
	{
		moleFraction[0] = 1.0;
	}
}

/*!
*******************************************************************************
* Property::xmass -- ���ʕ����E���q�ʊ֐�  �����������玿�ʕ����ƕ��q�ʂ̎Z�o
*
* ����(���͒l)
*     moleFraction    -- ���������z��
* ����(�o�͒l)
*     massFraction    -- ���ʕ����z��
*     molecularWeight -- ���q�� g/mol
*
* ����
*     1)�P�����̏ꍇxmassDll���玿�ʕ������Z�o�ł��Ȃ�(�����s��)
*       ���̖���������邽�߁C�P�����̏ꍇ���ʕ����ɒ���1.0�������Ă���
*******************************************************************************
*/
void DXProperty::xmass(double *moleFraction, double *massFraction, double &molecularWeight)
{
	// ���ʕ����E���q��
	xmassDll(moleFraction, massFraction, molecularWeight);
	// �P�����̏ꍇ
	if(numberOfComponents_ == 1)
	{
		massFraction[0] = 1.0;
	}
}

/********************************************************************************
 * Property::molecularWeight -- ���q�ʂ̎Z�o�Dmol����kg�ւ̕ϊ��ŕK�v           *
 *                                                                              *
 * �߂�l                                                                       *
 *     molecularWeight  -- ���q��                                               *
 *                                                                              *
 * ����                                                                         *
 *     1)�������쓮���̂ł̎g�p�͂ł��Ȃ�                                       *
 *                                                                              *
 ********************************************************************************/
double DXProperty::molecularWeight_f()
{
	// ��������
	double moleFraction[MaxComponent] = {0.0};
	// �e�����̃��������i�����ł͒P�����̂݁j
	moleFraction[0] = 1.0;
	// ���ʕ���
	double massFraction[MaxComponent] = {0.0};
	// ���q�� g/mol
	double molecularWeight = 0.0;
	// ���q�ʊ֐�
	xmassDll(moleFraction, massFraction, molecularWeight);
	return (molecularWeight);
}



/*!
*******************************************************************************
* DXProperty::critilcalPoint -- �ՊE��Ԋ֐�  �ՊE��Ԃ̎Z�o
*
* ����(�o�͒l)
*     crit -- ���̃N���X
*
* ����
*     1)��M�C�M�`���C�S���͒�`�ł��Ȃ����߁C0���o�͂���
********************************************************************************
*/
void DXProperty::criticalPoint(DXFluid &crit)
{
	// �G���[�ԍ�
	long icritical = 0;

	// �G���[���b�Z�[�W
	char cError[LError+1];
	strcpy_s(cError, "Ok");

	// �ՊE��Ԋ֐�
	critpDll(moleFraction, crit.T, crit.P, crit.rho, icritical, cError, LError);

	// �_�~�[
	double dummy[6] = {0.0};

	// ��G���^���s�E��G���g���s
	thermDll(crit.T, crit.rho, moleFraction, dummy[0], dummy[1], crit.h, crit.s, dummy[2], dummy[3], dummy[4], dummy[5]);

	// �P�ʕϊ�
	crit.T   -= 273.15;          // K -> deg.C
	crit.rho *= molecularWeight; // mol/L -> kg/m^3
	crit.h   /= molecularWeight; // J/mol -> kJ/kg
	crit.s   /= molecularWeight; // J/(mol*K) -> kJ/(kg*K)

	// �l�̑��
	crit.x   = 0.0;
	crit.cp  = 0.0;
	crit.thc   = 0.0;
	crit.visc  = 0.0;
	crit.st   = 0.0;
}
/*!
*******************************************************************************
* DXProperty::state_tp -- �P����Ԋ֐�(���x����͓���)
*                       �C�ӂ̈��ͥ���x�ɂ�����P����Ԃ̎Z�o
*
* ����(���͒l)
*     temperature -- ���x deg.C
*     pressure    -- ���� KPa
* ����(�o�͒l)
*     ref         -- ���̃N���X
*******************************************************************************
*/
void DXProperty::state_tp(double temperature, double pressure)
{
	DXFluid ref;

	// �l�̑��
	ref.T = temperature;
	ref.P = pressure;

	// �P�ʕϊ�
	ref.T += 273.15;// deg.C -> K

	// �_�~�[
	double dummy[5] = {0.0};
	double xldummy[MaxComponent] = {0.0};
	double xvdummy[MaxComponent] = {0.0};

	// �G���[�ԍ�
	long itpflsh = 0;

	// �G���[���b�Z�[�W
	char cError[LError+1];
	strcpy_s(cError, "Ok");

	// ��Ԋ֐�
//	tpflshDll(ref.T, ref.P, moleFraction, ref.rho, dummy[0], dummy[1], xldummy, xvdummy,
//	          ref.x, dummy[2], ref.h, ref.s, dummy[3], ref.cp, dummy[4], itpflsh, cError, LError);
	tpflshDll(ref.T, ref.P, moleFraction, ref.rho, dummy[0], dummy[1], xldummy, xvdummy,
			  ref.x, dummy[2], ref.h, ref.s, ref.cv, ref.cp, dummy[4], itpflsh, cError, LError);

	// �G���[�ԍ�
	long itrnprp = 0;

	// �ړ������֐�
	trnprpDll(ref.T, ref.rho, moleFraction, ref.visc, ref.thc, itrnprp, cError, LError);

	// �P�ʕϊ�
	ref.T   -= 273.15;         // K -> deg.C
	ref.rho *= molecularWeight;// mol/L -> kg/m^3
	ref.h   /= molecularWeight;// J/mol -> kJ/kg
	ref.s   /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	ref.cp  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	ref.cv  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	Rc = ref;
}
/*!
*******************************************************************************
* DXProperty::state_ph -- �C�ӏ�Ԋ֐�(���ͥ��G���^���s����)
*                       �C�ӂ̈��ͥ��G���^���s�ɂ������Ԃ̎Z�o
*
* ����(���͒l)
*     pressure    -- ���� KPa
*     enthalpy    -- �G���^���s kJ/kg
* ����(�o�͒l)
*     ref         -- ���̃N���X
*
* ����
*     1)�񑊈�O�ł̊����x�͖��Ӗ��Ȓl
*     2)�񑊈���ł̓o���N�̒l���o�͂���
*     3)�񑊈���ł͔�M�C�M�`���C�S���͒�`�ł��Ȃ����߁C0���o�͂���
*******************************************************************************
*/
void DXProperty::state_ph(double pressure, double enthalpy)
{

	DXFluid ref;


	// �l�̑��
	ref.P = pressure;
	ref.h = enthalpy;

	// �P�ʕϊ�
	ref.h *= molecularWeight;// kJ/kg -> J/mol

	// �_�~�[
	double dummy[5] = {0.0};
	double xldummy[MaxComponent] = {0.0};
	double xvdummy[MaxComponent] = {0.0};

	// �G���[�ԍ�
	long iphflsh = 0;

	// �G���[���b�Z�[�W
	char cError[LError+1];
	strcpy_s(cError, "Ok");

	// ��Ԋ֐�
//	phflshDll(ref.P, ref.h, moleFraction, ref.T, ref.rho, dummy[0], dummy[1], xldummy, xvdummy,
//	          ref.x, dummy[2], ref.s, dummy[3], ref.cp, dummy[4], iphflsh, cError, LError);
	phflshDll(ref.P, ref.h, moleFraction, ref.T, ref.rho, dummy[0], dummy[1], xldummy, xvdummy,
		ref.x, dummy[2], ref.s, ref.cv, ref.cp, dummy[4], iphflsh, cError, LError);

	// �G���[�ԍ�
	long itrnprp = 0;

	// �ړ������֐�
	trnprpDll(ref.T, ref.rho, moleFraction, ref.visc, ref.thc, itrnprp, cError, LError);

	// �l�̑��
	if(0 < ref.x && ref.x < 1)
	{
		ref.cp   = 0.0;
		ref.thc  = 0.0;
		ref.visc = 0.0;
		ref.st   = 0.0;
	}

	// �P�ʕϊ�
	ref.T   -= 273.15;         // K -> deg.C
	ref.rho *= molecularWeight;// mol/L -> kg/m^3
	ref.h   /= molecularWeight;// J/mol -> kJ/kg
	ref.s   /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	ref.cp  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	ref.cv  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)

	Rc = ref;

}
/*!
*******************************************************************************
* DXProperty::sat_t -- �O�a��Ԋ֐��i���x���́j
*                    �C�ӂ̉��x�ɂ�����O�a�t�ƖO�a���C�̏�Ԃ̎Z�o
*
* ����(���͒l) 
*     temperature  -- ���x deg.C
* ����(�o�͒l) 
*     refl         -- �t�����̃N���X
*     refv         -- ���C�����̃N���X
*******************************************************************************
*/
void DXProperty::sat_t(double temperature)
{

	DXFluid refl,refv;

	// �_�~�[
	double dummy[7] = {0.0};
	double xldummy[MaxComponent] = {0.0};
	double xvdummy[MaxComponent] = {0.0};

	////////////////////
	// �O�a�t
	////////////////////

	// �l�̑��
	refl.T = temperature;

	// �P�ʕϊ�
	refl.T += 273.15;// deg.C -> K

	// �t�F�[�Y�t���b�O
	long liquid = 1;

	// �G���[�ԍ�
	long isattl = 0;
	long itrnprpl = 0;
	long isurft   = 0;

	// �G���[���b�Z�[�W
	char cErrorl[LError+1];
	strcpy_s(cErrorl, "Ok");

	// �O�a�t�֐�
	sattDll(refl.T, moleFraction, liquid, refl.P, refl.rho, dummy[0], xldummy, xvdummy, isattl, cErrorl, LError);
	// ��G���^���s�E��G���g���s�E�舳��M
	thermDll(refl.T, refl.rho, moleFraction, dummy[1], dummy[2], refl.h, refl.s, refl.cv, refl.cp, dummy[4], dummy[5]);
	// �ړ������֐�
	trnprpDll(refl.T, refl.rho, moleFraction, refl.visc, refl.thc, itrnprpl, cErrorl, LError);
	// �\�ʒ���
	surftDll(refl.T, dummy[6] , moleFraction, refl.st, isurft, cErrorl, LError);

	// �����x
	refl.x = 0.0;

	// �P�ʕϊ�
	refl.T   -= 273.15;         // K -> deg.C
	refl.rho *= molecularWeight;// mol/L -> kg/m^3
	refl.h   /= molecularWeight;// J/mol -> kJ/kg
	refl.s   /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	refl.cp  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	refl.cv  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)


	////////////////////
	// �O�a���C
	////////////////////

	// �l�̑��
	refv.T = temperature;

	// �P�ʕϊ�
	refv.T += 273.15;// deg.C -> K

	// �t�F�[�Y�t���b�O
	long vapor = 2;

	// �G���[�ԍ�
	long isattv = 0;
	long itrnprpv = 0;

	// �G���[���b�Z�[�W
	char cErrorv[LError+1];
	strcpy_s(cErrorv, "Ok");

	// �O�a�t�֐�
	sattDll(refv.T, moleFraction, vapor, refv.P, dummy[0], refv.rho, xldummy, xvdummy, isattv, cErrorv, LError);
	// ��G���^���s�E��G���g���s�E�舳��M
	thermDll(refv.T, refv.rho, moleFraction, dummy[1], dummy[2], refv.h, refv.s, refv.cv, refv.cp, dummy[4], dummy[5]);
	// �ړ������֐�
	trnprpDll(refv.T, refv.rho, moleFraction, refv.visc, refv.thc, itrnprpv, cErrorv, LError);

	// �����x
	refv.x = 1.0;

	// �P�ʕϊ�
	refv.T   -= 273.15;         // K -> deg.C
	refv.rho *= molecularWeight;// mol/L -> kg/m^3
	refv.h   /= molecularWeight;// J/mol -> kJ/kg
	refv.s   /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	refv.cp  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	refv.cv  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)

	Rl = refl;
	Rv = refv;
}

/*!
*******************************************************************************
* DXProperty::sat_p -- �O�a��Ԋ֐��i���͓��́j
*                    �C�ӂ̈��͂ɂ�����O�a�t�ƖO�a���C�̏�Ԃ̎Z�o
*
* ����(���͒l) 
*     pressure     -- ���� kPa
* ����(�o�͒l) 
*     refl         -- �t�����̃N���X
*     refv         -- ���C�����̃N���X
*******************************************************************************
*/
void DXProperty::sat_p(double pressure)
{
	DXFluid refl,refv;
	// �_�~�[
	double dummy[7] = {0.0};
	double xldummy[MaxComponent] = {0.0};
	double xvdummy[MaxComponent] = {0.0};

	////////////////////
	// �O�a�t
	////////////////////

	// �l�̑��
	refl.P = pressure;

	// �t�F�[�Y�t���b�O
	long liquid = 1;

	// �G���[�ԍ�
	long isatpl = 0;
	long itrnprpl = 0;
	long isurft   = 0;

	// �G���[���b�Z�[�W
	char cErrorl[LError+1];
	strcpy_s(cErrorl, "Ok");

	// �O�a�t�֐�
	satpDll(refl.P, moleFraction, liquid, refl.T, refl.rho, dummy[0], xldummy, xvdummy, isatpl, cErrorl, LError);
	// ��G���^���s�E��G���g���s�E�舳��M
	thermDll(refl.T, refl.rho, moleFraction, dummy[1], dummy[2], refl.h, refl.s, refl.cv, refl.cp, dummy[4], dummy[5]);
	// �ړ������֐�
	trnprpDll(refl.T, refl.rho, moleFraction, refl.visc, refl.thc, itrnprpl, cErrorl, LError);
	// �\�ʒ���
	surftDll(refl.T, dummy[6] , moleFraction, refl.st, isurft, cErrorl, LError);

	// �����x
	refl.x = 0.0;

	// �P�ʕϊ�
	refl.T   -= 273.15;         // K -> deg.C
	refl.rho *= molecularWeight;// mol/L -> kg/m^3
	refl.h   /= molecularWeight;// J/mol -> kJ/kg
	refl.s   /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	refl.cp  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	refl.cv  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)

	////////////////////
	// �O�a���C
	////////////////////

	// �l�̑��
	refv.P = pressure;

	// �t�F�[�Y�t���b�O
	long vapor = 2;

	// �G���[�ԍ�
	long isatpv = 0;
	long itrnprpv = 0;

	// �G���[���b�Z�[�W
	char cErrorv[LError+1];
	strcpy_s(cErrorv, "Ok");

	// �O�a���C�֐�
	satpDll(refv.P, moleFraction, vapor, refv.T, dummy[0], refv.rho, xldummy, xvdummy, isatpv, cErrorv, LError);
	// ��G���^���s�E��G���g���s�E�舳��M
	thermDll(refv.T, refv.rho, moleFraction, dummy[1], dummy[2], refv.h, refv.s, refv.cv, refv.cp, dummy[4], dummy[5]);
	// �ړ������֐�
	trnprpDll(refv.T, refv.rho, moleFraction, refv.visc, refv.thc, itrnprpv, cErrorv, LError);

	// �����x
	refv.x = 1.0;

	// �P�ʕϊ�
	refv.T   -= 273.15;         // K -> deg.C
	refv.rho *= molecularWeight;// mol/L -> kg/m^3
	refv.h   /= molecularWeight;// J/mol -> kJ/kg
	refv.s   /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	refv.cp  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)
	refv.cv  /= molecularWeight;// J/(mol*K) -> kJ/(kg*K)

	Rl = refl;
	Rv = refv;
}

/*!
*******************************************************************************
* DXProperty::surfaceTension_t -- �\�ʒ��͊֐��D                               
*                   
* ����(���͒l)     
*     temperature -- ���x deg.C 
* �߂�l               
*     surfaceTension -- �\�ʒ��� N/m 
* 
*******************************************************************************
*/
double DXProperty::surfaceTension_t(double temperature)
{
	// �P�ʕϊ�
	temperature += 273.15;// deg.C -> K

	// ��������
	double moleFraction[MaxComponent] = {0.0};
	// �e�����̃��������i�����ł͒P�����̂݁j
	moleFraction[0] = 1.0;
	// �t�F�[�Y�t���b�O�i�P�����̏ꍇ��1or2�j
	long phase = 1;

	// �C���E�t�����������i�P�����̏ꍇ�݂͌��ɓ������j
	double liquidFraction[MaxComponent] = {0.0};
	double vaporFraction[MaxComponent] = {0.0};
	
	// �O�a����
	double saturationPressure;
	
	// �C���E�t�����x
	double liquidDensity, vaporDensity;

	// �G���[�ԍ�
	long isatt = 0;
	// �G���[���b�Z�[�W
	char cErrorMessage[LError+1];
	strncpy_s(cErrorMessage, "Ok", sizeof(cErrorMessage)-1);

	// �O�a��Ԋ֐�
	sattDll(temperature, moleFraction, phase, saturationPressure, liquidDensity, vaporDensity,
	        liquidFraction, vaporFraction, isatt, cErrorMessage, LError);

	// �G���[����
	if(isatt != 0)
	{
		std::cout << "sattDll�֐��G���[" << std::endl;
		std::cout << isatt << " : " << cErrorMessage << std::endl;
		// DLL�̊J��
		freeDll();
		exit(1);
	}

	// �G���[�ԍ�
	long isurten = 0;

	// �\�ʒ��͊֐�
	double surfaceTension;
	surtenDll(temperature, liquidDensity, vaporDensity, liquidFraction, vaporFraction,
	          surfaceTension, isurten, cErrorMessage, LError);

	// �G���[����
	if(isurten != 0)
	{
		std::cout << "surten�֐��G���[" << std::endl;
		std::cout << isurten << " : " << cErrorMessage << std::endl;
		// DLL�̊J��
		freeDll();
		exit(1);
	}

	// �P�ʕϊ�
	temperature -= 273.15;           // K -> deg.C

	return (surfaceTension);
}




/*************************************************************************************************************************/
/*---------------------------------------------------------�v�Z���x����̗̈�--------------------------------------------*/
/*************************************************************************************************************************/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	state_ph2,state_tp2,sat_p2,sat_t2�̎g����
//	�`�����`
//	300�~300�̊i�q�𗘗p���Đ��refprop�̊֐����v�Z��,���̓_�𗘗p���Ĕ���Ԃ��ďo�͂���֐�.
//	��Ɍv�Z�����Ă��邽�ߕ��ʂ�state_ph���������v�Z�������Ԃ������������.
//	���̂��߂��̊֐��������g�����߂ɂ�,�\�ߌv�Z�������Ă��̊i�q�_�����̃f�[�^�����Ȃ��Ă͂Ȃ�Ȃ�.
//	
//	�`�g�����`
//	1.setup��
//	�����ł͗�}�Ƃ���CO2���g���ꍇ������.
//------------------------------------------------------------------
//	DXProperty Ref;
//	Ref.reibai = "CO2.FLD";
//	Ref.joutai = "IIR";
//	Ref.LoadDLL("refprop.DLL");
//	Ref.setup();
//-----------------�������܂ł����ƈꏏ---------------------------
//
//	�܂��͊i�q�f�[�^�𐻍삷��.
//------------------------------------------------------------------
//	Ref.make_table( 200, 600, 2000, 20000, -20, 200); //�O��̓G���^���s�[�̕�,����͈��͂̕�,���͉��x�̕�.
//------------------------------------------------------------------
//
//	���̂悤�Ɋi�q�̕������Ă����Ă��,���̕��œ��Ԋu��300������,���̓_���v�Z��make_table�̏ꍇ��refdate.csv�Ƃ������̂�,�ۑ������.
//	���ɂ��̍��ꂽ�t�@�C����fluids�t�H���_�ɗl�X�ȗ�}�̊֐��������Ă���悤��,fluids2�Ƃ����t�H���_�������ō��,
//	���̒��ɂ���csv�t�@�C��������.
//	(����)�ق��̗�}�̏�������Ƃ�,make_table�����͏��refdata.csv�Ƃ������O�ŕۑ�����̂�,�������O��ύX��,�����ŊǗ����Ă������ق����悢
//	
//	2.������
//------------------------------------------------------------------
//	DXProperty Ref;
//	Ref.reibai = "CO2.FLD";
//	Ref.joutai = "IIR";
//	Ref.LoadDLL("refprop.DLL");
//	Ref.setup();
//
//	Ref.init_table();//table�̏�����
//	Ref.load_table("refdate.csv");//csv�t�@�C���̓ǂݍ���
//
//	Ref.state_ph2( 300, 5000 );//�g�p�\.
//------------------------------------------------------------------
//
//	3.���ӓ_
//	�ꌩ�v�Z���x���オ�邽�ߕ֗��Ȋ֐��Ɍ����邪���̓_�ɒ��ӂ��邱��.
//	* �ߎ��l����邽��,�l�̐��x���r���Ȃ�(�i�q�̕���傫������΂���قǍr���Ȃ�)
//	* �����v�Z�Ȃǂł��̊i�q�̕�����l���o���ꍇ�v�Z�s�\�ƂȂ�.
//
//		�ȏ� 2009/05/21 nakajima
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/*!
*******************************************************************************
 DXProperty::set_step( int h, int p, int t ) -- �i�q�̐���set����֐�.
*******************************************************************************
 */
void DXProperty::set_step( int h , int p , int t){
	rep_h = h;
	rep_p = p;
	rep_t = t;

}

/*!
*******************************************************************************
 DXProperty::init_table() -- table������������֐�.
*******************************************************************************
*/
void DXProperty::init_table(){
	int i;

	data_table = new DXFluid*[rep_h];
	for( i = 0 ; i < rep_h ; i++ ){
		data_table[i] = new DXFluid[rep_p];
	}

	data_table_tp = new DXFluid*[rep_t];
	for( i = 0 ; i < rep_t ; i++ ){
		data_table_tp[i] = new DXFluid[rep_p];
	}


	data_table_p = new DXFluid*[2];
	data_table_p[0] = new DXFluid[rep_p];
	data_table_p[1] = new DXFluid[rep_p];

	data_table_t = new DXFluid*[2];
	data_table_t[0] = new DXFluid[rep_t];
	data_table_t[1] = new DXFluid[rep_t];


}
void DXProperty::init_table_du(){
	int i;

	data_table_du = new DXFluid*[rep_u];
	for( i = 0 ; i < rep_u ; i++ ){
		data_table_du[i] = new DXFluid[rep_d];
	}


}
/*!
*******************************************************************************
 DXProperty::load_table( std::string fname) -- csv�t�@�C���ǂݍ��݊֐�.

 ����(���͒l)                                                                 
     fname    -- csv�̃t�@�C���̖��O
*******************************************************************************
*/
void DXProperty::load_table( std::string fname){

	CCSVFile cf;
	ostringstream os;
	os << PATH << "fluids2/" << fname;
	if(!cf.LoadBuffer( os.str() ) )
	{
		return;
	}



	int i, j;
	int k = 0;


	int m = 0 ;
	string temp_massage;
	while( 1 ){

		temp_massage = cf.GetCSVBuffer(k,0);

		if( temp_massage[0] == 47 && temp_massage[0] == 47 ){
//			std::cout << "comment" << std::endl;
			k++;
		}else{
//			std::cout << "not comment" << std::endl;
			break;
		}


		m++;
		if( 1000 < m ){
			std::cout << "long comment" << std::endl;
			exit(0);
		}
	}


	min_h = atof(cf.GetCSVBuffer(k++,0));
	max_h = atof(cf.GetCSVBuffer(k++,0));
	min_P = atof(cf.GetCSVBuffer(k++,0));
	max_P = atof(cf.GetCSVBuffer(k++,0));
	min_T = atof(cf.GetCSVBuffer(k++,0));
	max_T = atof(cf.GetCSVBuffer(k++,0));

	k++;//���x���̃X�L�b�v

	int max = rep_p * rep_h;
	for( i = 0 ; i < rep_h ; i++ ){
		for( j = 0 ; j < rep_p ; j++ , k++){
			data_table[i][j].T   =  atof(cf.GetCSVBuffer(k,0));
			data_table[i][j].P   =  atof(cf.GetCSVBuffer(k,1));
			data_table[i][j].rho =  atof(cf.GetCSVBuffer(k,2));
			data_table[i][j].h   =  atof(cf.GetCSVBuffer(k,3));
			data_table[i][j].s   =  atof(cf.GetCSVBuffer(k,4));
			data_table[i][j].x   =  atof(cf.GetCSVBuffer(k,5));
			data_table[i][j].cp  =  atof(cf.GetCSVBuffer(k,6));
			data_table[i][j].thc =  atof(cf.GetCSVBuffer(k,7));
			data_table[i][j].visc=  atof(cf.GetCSVBuffer(k,8));
			data_table[i][j].st  =  atof(cf.GetCSVBuffer(k,9));
			data_table[i][j].cv  =  atof(cf.GetCSVBuffer(k,10));

		}
		cout << "\rNOW LOADING " << double(k)/double(max);
	}

	std::cout << std::endl;

	for( i = 0 ; i < rep_t ; i++ ){
		for( j = 0 ; j < rep_p ; j++ , k++){
			data_table_tp[i][j].T   =  atof(cf.GetCSVBuffer(k,0));
			data_table_tp[i][j].P   =  atof(cf.GetCSVBuffer(k,1));
			data_table_tp[i][j].rho =  atof(cf.GetCSVBuffer(k,2));
			data_table_tp[i][j].h   =  atof(cf.GetCSVBuffer(k,3));
			data_table_tp[i][j].s   =  atof(cf.GetCSVBuffer(k,4));
			data_table_tp[i][j].x   =  atof(cf.GetCSVBuffer(k,5));
			data_table_tp[i][j].cp  =  atof(cf.GetCSVBuffer(k,6));
			data_table_tp[i][j].thc =  atof(cf.GetCSVBuffer(k,7));
			data_table_tp[i][j].visc=  atof(cf.GetCSVBuffer(k,8));
			data_table_tp[i][j].st  =  atof(cf.GetCSVBuffer(k,9));
			data_table_tp[i][j].cv  =  atof(cf.GetCSVBuffer(k,10));

		}
		cout << "\rNOW LOADING " << double(k)/double(max);
	}

	i = k;
	for( j = 0 ; j < rep_p ; j++ ){
		data_table_p[0][j].T   =  atof(cf.GetCSVBuffer(k,0));
		data_table_p[0][j].P   =  atof(cf.GetCSVBuffer(k,1));
		data_table_p[0][j].rho =  atof(cf.GetCSVBuffer(k,2));
		data_table_p[0][j].h   =  atof(cf.GetCSVBuffer(k,3));
		data_table_p[0][j].s   =  atof(cf.GetCSVBuffer(k,4));
		data_table_p[0][j].x   =  atof(cf.GetCSVBuffer(k,5));
		data_table_p[0][j].cp  =  atof(cf.GetCSVBuffer(k,6));
		data_table_p[0][j].thc =  atof(cf.GetCSVBuffer(k,7));
		data_table_p[0][j].visc=  atof(cf.GetCSVBuffer(k,8));
		data_table_p[0][j].st  =  atof(cf.GetCSVBuffer(k,9));
		data_table_p[0][j].cv  =  atof(cf.GetCSVBuffer(k,10));
		k++;
		data_table_p[1][j].T   =  atof(cf.GetCSVBuffer(k,0));
		data_table_p[1][j].P   =  atof(cf.GetCSVBuffer(k,1));
		data_table_p[1][j].rho =  atof(cf.GetCSVBuffer(k,2));
		data_table_p[1][j].h   =  atof(cf.GetCSVBuffer(k,3));
		data_table_p[1][j].s   =  atof(cf.GetCSVBuffer(k,4));
		data_table_p[1][j].x   =  atof(cf.GetCSVBuffer(k,5));
		data_table_p[1][j].cp  =  atof(cf.GetCSVBuffer(k,6));
		data_table_p[1][j].thc =  atof(cf.GetCSVBuffer(k,7));
		data_table_p[1][j].visc=  atof(cf.GetCSVBuffer(k,8));
		data_table_p[1][j].st  =  atof(cf.GetCSVBuffer(k,9));
		data_table_p[1][j].cv  =  atof(cf.GetCSVBuffer(k,10));
		k++;
	}

	i = k;
	for( j = 0 ; j < rep_t ; j++ ){
		data_table_t[0][j].T   =  atof(cf.GetCSVBuffer(k,0));
		data_table_t[0][j].P   =  atof(cf.GetCSVBuffer(k,1));
		data_table_t[0][j].rho =  atof(cf.GetCSVBuffer(k,2));
		data_table_t[0][j].h   =  atof(cf.GetCSVBuffer(k,3));
		data_table_t[0][j].s   =  atof(cf.GetCSVBuffer(k,4));
		data_table_t[0][j].x   =  atof(cf.GetCSVBuffer(k,5));
		data_table_t[0][j].cp  =  atof(cf.GetCSVBuffer(k,6));
		data_table_t[0][j].thc =  atof(cf.GetCSVBuffer(k,7));
		data_table_t[0][j].visc=  atof(cf.GetCSVBuffer(k,8));
		data_table_t[0][j].st  =  atof(cf.GetCSVBuffer(k,9));
		data_table_t[0][j].cv  =  atof(cf.GetCSVBuffer(k,10));
		k++;
		data_table_t[1][j].T   =  atof(cf.GetCSVBuffer(k,0));
		data_table_t[1][j].P   =  atof(cf.GetCSVBuffer(k,1));
		data_table_t[1][j].rho =  atof(cf.GetCSVBuffer(k,2));
		data_table_t[1][j].h   =  atof(cf.GetCSVBuffer(k,3));
		data_table_t[1][j].s   =  atof(cf.GetCSVBuffer(k,4));
		data_table_t[1][j].x   =  atof(cf.GetCSVBuffer(k,5));
		data_table_t[1][j].cp  =  atof(cf.GetCSVBuffer(k,6));
		data_table_t[1][j].thc =  atof(cf.GetCSVBuffer(k,7));
		data_table_t[1][j].visc=  atof(cf.GetCSVBuffer(k,8));
		data_table_t[1][j].st  =  atof(cf.GetCSVBuffer(k,9));
		data_table_t[1][j].cv  =  atof(cf.GetCSVBuffer(k,10));
		k++;
	}

	cout << "\nfinish" << endl;



}
void DXProperty::load_table2( std::string fname)
{
	FILE			*fp;
	ostringstream	os;
	double			dV0, dV1, dV2, dV3, dV4, dV5, dV6, dV7, dV8, dV9, dV10;
	char			szValue[11][32];
	char			temp_massage[300];
	int				i, j, ret;
	int				k = 0;
	int				m = 0;
	int				max = rep_p * rep_h;

	os << PATH << "fluids2/" << fname;

	std::string filename = os.str();
	fp = fopen(filename.c_str(), "r");
	if (fp == NULL){
		return;
	}



//	string temp_massage;
	while( 1 ){


//		temp_massage = cf.GetCSVBuffer(k,0);
//		ret = fscanf(fp,"%s",szValue[0]);
		fgets( temp_massage , 255 , fp );


		if( temp_massage[0] == 47 && temp_massage[1] == 47 ){
//			std::cout << "comment" << std::endl;
			k++;
		}else{
//			std::cout << "not comment" << std::endl;
			break;
		}


		m++;
		if( 1000 < m ){
			std::cout << "long comment" << std::endl;
			exit(0);
		}
	}




	min_h = atof(temp_massage);
//	ret = fscanf(fp,"%s",temp_massage);
	fgets( temp_massage , 255 , fp );
	max_h = atof(temp_massage);
	fgets( temp_massage , 255 , fp );
	min_P = atof(temp_massage);
	fgets( temp_massage , 255 , fp );
	max_P = atof(temp_massage);
	fgets( temp_massage , 255 , fp );
	min_T = atof(temp_massage);
	fgets( temp_massage , 255 , fp );
	max_T = atof(temp_massage);

	fgets( temp_massage , 255 , fp );



	for( i = 0 ; i < rep_h ; i++ ){
		for( j = 0 ; j < rep_p ; j++ , k++){
			ret = fscanf(fp, "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],", 
						szValue[0], szValue[1], szValue[2], szValue[3], szValue[4],
						szValue[5], szValue[6], szValue[7], szValue[8], szValue[9], szValue[10]);
			dV0  = atof(szValue[0]);
			dV1  = atof(szValue[1]);
			dV2  = atof(szValue[2]);
			dV3  = atof(szValue[3]);
			dV4  = atof(szValue[4]);
			dV5  = atof(szValue[5]);
			dV6  = atof(szValue[6]);
			dV7  = atof(szValue[7]);
			dV8  = atof(szValue[8]);
			dV9  = atof(szValue[9]);
			dV10 = atof(szValue[10]);
			if (ret != EOF){
				data_table[i][j].T   =  dV0;
				data_table[i][j].P   =  dV1;
				data_table[i][j].rho =  dV2;
				data_table[i][j].h   =  dV3;
				data_table[i][j].s   =  dV4;
				data_table[i][j].x   =  dV5;
				data_table[i][j].cp  =  dV6;
				data_table[i][j].thc =  dV7;
				data_table[i][j].visc=  dV8;
				data_table[i][j].st  =  dV9;
				data_table[i][j].cv  =  dV10;
			}
			else{
				cout << "read error!" << endl;
				return;
			}
		}

		if( i % 10 == 0 ){
			cout << "\rNOW LOADING " << double(k)/double(max);
		}
	}

	std::cout << std::endl;

	for( i = 0 ; i < rep_t ; i++ ){
		for( j = 0 ; j < rep_p ; j++ , k++){
			ret = fscanf(fp, "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],", 
						szValue[0], szValue[1], szValue[2], szValue[3], szValue[4],
						szValue[5], szValue[6], szValue[7], szValue[8], szValue[9], szValue[10]);
			dV0  = atof(szValue[0]);
			dV1  = atof(szValue[1]);
			dV2  = atof(szValue[2]);
			dV3  = atof(szValue[3]);
			dV4  = atof(szValue[4]);
			dV5  = atof(szValue[5]);
			dV6  = atof(szValue[6]);
			dV7  = atof(szValue[7]);
			dV8  = atof(szValue[8]);
			dV9  = atof(szValue[9]);
			dV10 = atof(szValue[10]);
			if (ret != EOF){
				data_table_tp[i][j].T   =  dV0;
				data_table_tp[i][j].P   =  dV1;
				data_table_tp[i][j].rho =  dV2;
				data_table_tp[i][j].h   =  dV3;
				data_table_tp[i][j].s   =  dV4;
				data_table_tp[i][j].x   =  dV5;
				data_table_tp[i][j].cp  =  dV6;
				data_table_tp[i][j].thc =  dV7;
				data_table_tp[i][j].visc=  dV8;
				data_table_tp[i][j].st  =  dV9;
				data_table_tp[i][j].cv  =  dV10;
			}
			else{
				cout << "read error!" << endl;
				return;
			}
		}
		if( i % 10 == 0 ){
			cout << "\rNOW LOADING " << double(k)/double(max);
		}
	}

	i = k;
	for( j = 0 ; j < rep_p ; j++ ){
		ret = fscanf(fp, "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],", 
					szValue[0], szValue[1], szValue[2], szValue[3], szValue[4],
					szValue[5], szValue[6], szValue[7], szValue[8], szValue[9], szValue[10]);
		dV0  = atof(szValue[0]);
		dV1  = atof(szValue[1]);
		dV2  = atof(szValue[2]);
		dV3  = atof(szValue[3]);
		dV4  = atof(szValue[4]);
		dV5  = atof(szValue[5]);
		dV6  = atof(szValue[6]);
		dV7  = atof(szValue[7]);
		dV8  = atof(szValue[8]);
		dV9  = atof(szValue[9]);
		dV10 = atof(szValue[10]);
		if (ret != EOF){
			data_table_p[0][j].T   =  dV0;
			data_table_p[0][j].P   =  dV1;
			data_table_p[0][j].rho =  dV2;
			data_table_p[0][j].h   =  dV3;
			data_table_p[0][j].s   =  dV4;
			data_table_p[0][j].x   =  dV5;
			data_table_p[0][j].cp  =  dV6;
			data_table_p[0][j].thc =  dV7;
			data_table_p[0][j].visc=  dV8;
			data_table_p[0][j].st  =  dV9;
			data_table_p[0][j].cv  =  dV10;
		}
		else{
			cout << "read error!" << endl;
			return;
		}
		k++;
		ret = fscanf(fp, "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],", 
					szValue[0], szValue[1], szValue[2], szValue[3], szValue[4],
					szValue[5], szValue[6], szValue[7], szValue[8], szValue[9], szValue[10]);
		dV0  = atof(szValue[0]);
		dV1  = atof(szValue[1]);
		dV2  = atof(szValue[2]);
		dV3  = atof(szValue[3]);
		dV4  = atof(szValue[4]);
		dV5  = atof(szValue[5]);
		dV6  = atof(szValue[6]);
		dV7  = atof(szValue[7]);
		dV8  = atof(szValue[8]);
		dV9  = atof(szValue[9]);
		dV10 = atof(szValue[10]);
		if (ret != EOF){
			data_table_p[1][j].T   =  dV0;
			data_table_p[1][j].P   =  dV1;
			data_table_p[1][j].rho =  dV2;
			data_table_p[1][j].h   =  dV3;
			data_table_p[1][j].s   =  dV4;
			data_table_p[1][j].x   =  dV5;
			data_table_p[1][j].cp  =  dV6;
			data_table_p[1][j].thc =  dV7;
			data_table_p[1][j].visc=  dV8;
			data_table_p[1][j].st  =  dV9;
			data_table_p[1][j].cv  =  dV10;
		}
		else{
			cout << "read error!" << endl;
			return;
		}
		k++;
	}

	i = k;
	for( j = 0 ; j < rep_t ; j++ ){
		ret = fscanf(fp, "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],", 
					szValue[0], szValue[1], szValue[2], szValue[3], szValue[4],
					szValue[5], szValue[6], szValue[7], szValue[8], szValue[9], szValue[10]);
		dV0  = atof(szValue[0]);
		dV1  = atof(szValue[1]);
		dV2  = atof(szValue[2]);
		dV3  = atof(szValue[3]);
		dV4  = atof(szValue[4]);
		dV5  = atof(szValue[5]);
		dV6  = atof(szValue[6]);
		dV7  = atof(szValue[7]);
		dV8  = atof(szValue[8]);
		dV9  = atof(szValue[9]);
		dV10 = atof(szValue[10]);
		if (ret != EOF){
			data_table_t[0][j].T   =  dV0;
			data_table_t[0][j].P   =  dV1;
			data_table_t[0][j].rho =  dV2;
			data_table_t[0][j].h   =  dV3;
			data_table_t[0][j].s   =  dV4;
			data_table_t[0][j].x   =  dV5;
			data_table_t[0][j].cp  =  dV6;
			data_table_t[0][j].thc =  dV7;
			data_table_t[0][j].visc=  dV8;
			data_table_t[0][j].st  =  dV9;
			data_table_t[0][j].cv  =  dV10;
		}
		else{
			cout << "read error!" << endl;
			return;
		}
		k++;
		ret = fscanf(fp, "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],", 
					szValue[0], szValue[1], szValue[2], szValue[3], szValue[4],
					szValue[5], szValue[6], szValue[7], szValue[8], szValue[9], szValue[10]);
		dV0  = atof(szValue[0]);
		dV1  = atof(szValue[1]);
		dV2  = atof(szValue[2]);
		dV3  = atof(szValue[3]);
		dV4  = atof(szValue[4]);
		dV5  = atof(szValue[5]);
		dV6  = atof(szValue[6]);
		dV7  = atof(szValue[7]);
		dV8  = atof(szValue[8]);
		dV9  = atof(szValue[9]);
		dV10 = atof(szValue[10]);
		if (ret != EOF){
			data_table_t[1][j].T   =  dV0;
			data_table_t[1][j].P   =  dV1;
			data_table_t[1][j].rho =  dV2;
			data_table_t[1][j].h   =  dV3;
			data_table_t[1][j].s   =  dV4;
			data_table_t[1][j].x   =  dV5;
			data_table_t[1][j].cp  =  dV6;
			data_table_t[1][j].thc =  dV7;
			data_table_t[1][j].visc=  dV8;
			data_table_t[1][j].st  =  dV9;
			data_table_t[1][j].cv  =  dV10;
		}
		else{
			cout << "read error!" << endl;
			return;
		}
		k++;
	}
	fclose(fp);

	cout << "\nfinish" << endl;
}


void DXProperty::load_table_du( std::string fname){

	CCSVFile cf;
	ostringstream os;
	os << "fluids2/" << fname;
	if(!cf.LoadBuffer( os.str() ) )
	{
		return;
	}

	int i, j;
	int k = 0;
	int max = rep_d * rep_u;
	for( i = 0 ; i < rep_u ; i++ ){
		for( j = 0 ; j < rep_d ; j++ , k++){
			data_table_du[i][j].T   =  atof(cf.GetCSVBuffer(k,0));
			data_table_du[i][j].P   =  atof(cf.GetCSVBuffer(k,1));
			data_table_du[i][j].rho =  atof(cf.GetCSVBuffer(k,2));
			data_table_du[i][j].h   =  atof(cf.GetCSVBuffer(k,3));
			data_table_du[i][j].s   =  atof(cf.GetCSVBuffer(k,4));
			data_table_du[i][j].x   =  atof(cf.GetCSVBuffer(k,5));
			data_table_du[i][j].cp  =  atof(cf.GetCSVBuffer(k,6));
			data_table_du[i][j].thc =  atof(cf.GetCSVBuffer(k,7));
			data_table_du[i][j].visc=  atof(cf.GetCSVBuffer(k,8));
			data_table_du[i][j].st  =  atof(cf.GetCSVBuffer(k,9));
			data_table_du[i][j].cv  =  atof(cf.GetCSVBuffer(k,10));
			data_table_du[i][j].u   =  atof(cf.GetCSVBuffer(k,11));

		}
		cout << "\rNOW LOADING " << double(k)/double(max);
	}

	cout << "\nfinish" << endl;



}
/*!
*******************************************************************************
 DXProperty::DXProperty::make_table -- csv����֐�.
 ����(���͒l)                                                                 
     start_h    -- �G���^���s�[�J�n�ʒu
     end_h	  -- �G���^���s�[�I���ʒu
     start_p	  -- ���͊J�n�ʒu
     end_p	  -- ���͏I���ʒu
     start_t	  -- ���x�J�n�ʒu
     end_t	  -- ���x�I���ʒu
*******************************************************************************
*/
void DXProperty::make_table(double _start_h , double _end_h , double _start_p , double _end_p , double _start_t , double _end_t ){

	double start_h;
	double start_p;
	double start_t;
	double end_h;
	double end_p;
	double end_t;

	start_h = _start_h - ( _end_h - _start_h ) * 0.005;
	end_h   = _end_h   + ( _end_h - _start_h ) * 0.005;
	start_p = _start_p - ( _end_p - _start_p ) * 0.005;
	end_p   = _end_p   + ( _end_p - _start_p ) * 0.005;
	start_t = _start_t - ( _end_t - _start_t ) * 0.005;
	end_t   = _end_t   + ( _end_t - _start_t ) * 0.005;

	if( start_p < _start_p / 2.0 ){
		start_p = _start_p / 2.0;
	}


	int h_num = 0;
	int p_num = 0;
	int t_num = 0;

	double step_h;
	double step_p;
	double step_t;

	step_h = ( end_h - start_h ) / double(rep_h);
//	step_p = ( end_p - start_p ) / double(rep_p);
	step_p = pow( ( end_p / start_p ) , 1.0 / double(rep_p-1) );

	step_t = ( end_t - start_t ) / double(rep_t);

	double h;
	double p;
	double t;

	int pp;
	int hh;
	int tt;



	double crt_pressure;
	criticalPoint(Rc);
	crt_pressure = Rc.P;

	refdata.open( "refdata.csv" , ios::out );
	if( refdata.fail() ){
		std::cout << "file open error" << endl;
		exit(1);    // �t�@�C�����J���܂���ł����B
	}


	refdata << "// DXProperty table data " << endl;
	refdata << "// enthalpy lange "    << _start_h << " - " << _end_h << " kJ/kg " << endl;
	refdata << "// pressure lange "    << _start_p << " - " << _end_p << " kPa " << endl;
	refdata << "// temperature lange " << _start_t << " - " << _end_t << " deg.C " << endl;


	refdata << _start_h << endl;
	refdata << _end_h   << endl;
	refdata << _start_p << endl;
	refdata << _end_p   << endl;
	refdata << _start_t << endl;
	refdata << _end_t   << endl;


	refdata << "//";
	refdata << "temperature[degC],";
	refdata << "pressure[kPa],";
	refdata << "density[kg/m3],";
	refdata << "enthalpy[kJ/kg],";
	refdata << "entorpy[kJ/kgK],";
	refdata << "vapor quality[-],";
	refdata << "CP[kJ/kgK],";
	refdata << "thc[W/(m*K)],";
	refdata << "visc[uPa*s],";
	refdata << "st[N/m],";
	refdata << "CV[kJ/kgK],";
	refdata << endl;




	refdata << std::setprecision(20);

	hh = 0;
	for( h = start_h ; hh < rep_h ; h += step_h , hh++ ){
		pp = 0;
//		for( p = start_p ; pp < rep_p ; p += step_p  ,pp++ ){
		for( p = start_p ; pp < rep_p ; p *= step_p  ,pp++ ){
			state_ph( p , h );

			if( p < crt_pressure ){
				sat_p( p );
				Rc.x = ( Rc.h - Rl.h ) / ( Rv.h - Rl.h );

				if( 0.0 <= Rc.x && Rc.x <= 1.0 ){
					Rc.cp = Rl.cp*(1.0-Rc.x)+Rv.cp*(Rc.x);
					Rc.cv = Rl.cv*(1.0-Rc.x)+Rv.cv*(Rc.x);
					Rc.visc = Rl.visc*(1.0-Rc.x)+Rv.visc*(Rc.x);
					Rc.thc = Rl.thc*(1.0-Rc.x)+Rv.thc*(Rc.x);
				}
			}



			data_table[hh][pp] = Rc;

			refdata << data_table[hh][pp].T << ",";
			refdata << data_table[hh][pp].P << ",";
			refdata << data_table[hh][pp].rho << ",";
			refdata << data_table[hh][pp].h << ",";
			refdata << data_table[hh][pp].s << ",";
			refdata << data_table[hh][pp].x << ",";
			refdata << data_table[hh][pp].cp << ",";
			refdata << data_table[hh][pp].thc << ",";
			refdata << data_table[hh][pp].visc << ",";
			refdata << data_table[hh][pp].st << ",";
			refdata << data_table[hh][pp].cv << ",";
			refdata << endl;
		}
		std::cout << "\r" << (double)hh/(double)rep_h ;
	}
	std::cout << std::endl;

	tt = 0;
	for( t = start_t ; tt < rep_t ; t += step_t , tt++ ){
		pp = 0;
//		for( p = start_p ; pp < rep_p ; p += step_p  ,pp++ ){
		for( p = start_p ; pp < rep_p ; p *= step_p  ,pp++ ){
			state_tp( t , p );

			if( p < crt_pressure ){
				sat_p( p );
				Rc.x = ( Rc.h - Rl.h ) / ( Rv.h - Rl.h );

				if( 0.0 <= Rc.x && Rc.x <= 1.0 ){
					Rc.cp = Rl.cp*(1.0-Rc.x)+Rv.cp*(Rc.x);
					Rc.cv = Rl.cv*(1.0-Rc.x)+Rv.cv*(Rc.x);
				}
			}

			data_table[tt][pp] = Rc;

			refdata << data_table[tt][pp].T << ",";
			refdata << data_table[tt][pp].P << ",";
			refdata << data_table[tt][pp].rho << ",";
			refdata << data_table[tt][pp].h << ",";
			refdata << data_table[tt][pp].s << ",";
			refdata << data_table[tt][pp].x << ",";
			refdata << data_table[tt][pp].cp << ",";
			refdata << data_table[tt][pp].thc << ",";
			refdata << data_table[tt][pp].visc << ",";
			refdata << data_table[tt][pp].st << ",";
			refdata << data_table[tt][pp].cv << ",";
			refdata << endl;
		}
		std::cout << "\r" << (double)tt/(double)rep_h ;
	}

	criticalPoint(Rc);
	step_p = (Rc.P-130.0) / double(rep_p);
	pp = 0;
	for( p = Rc.P-120.0 ; pp < rep_p ; p -= step_p , pp++ ){
		sat_p( p );

		refdata << Rl.T << ",";
		refdata << Rl.P << ",";
		refdata << Rl.rho << ",";
		refdata << Rl.h << ",";
		refdata << Rl.s << ",";
		refdata << Rl.x << ",";
		refdata << Rl.cp << ",";
		refdata << Rl.thc << ",";
		refdata << Rl.visc << ",";
		refdata << Rl.st << ",";
		refdata << Rl.cv << ",";
		refdata << endl;

		refdata << Rv.T << ",";
		refdata << Rv.P << ",";
		refdata << Rv.rho << ",";
		refdata << Rv.h << ",";
		refdata << Rv.s << ",";
		refdata << Rv.x << ",";
		refdata << Rv.cp << ",";
		refdata << Rv.thc << ",";
		refdata << Rv.visc << ",";
		refdata << Rv.st << ",";
		refdata << Rv.cv << ",";
		refdata << endl;

	}

	criticalPoint(Rc);
	step_t = ( (Rc.T-3.0) - start_t ) / double(rep_t);
	tt = 0;
	for( t = Rc.T-3.0 ; tt < rep_t ; t -= step_t , tt++ ){
		sat_t( t );

		refdata << Rl.T << ",";
		refdata << Rl.P << ",";
		refdata << Rl.rho << ",";
		refdata << Rl.h << ",";
		refdata << Rl.s << ",";
		refdata << Rl.x << ",";
		refdata << Rl.cp << ",";
		refdata << Rl.thc << ",";
		refdata << Rl.visc << ",";
		refdata << Rl.st << ",";
		refdata << Rl.cv << ",";
		refdata << endl;

		refdata << Rv.T << ",";
		refdata << Rv.P << ",";
		refdata << Rv.rho << ",";
		refdata << Rv.h << ",";
		refdata << Rv.s << ",";
		refdata << Rv.x << ",";
		refdata << Rv.cp << ",";
		refdata << Rv.thc << ",";
		refdata << Rv.visc << ",";
		refdata << Rv.st << ",";
		refdata << Rv.cv << ",";
		refdata << endl;

	}

	std::cout << "�e�[�u���쐬�I��" << std::endl;


}

void DXProperty::make_table_du(double start_d , double end_d , double start_u , double end_u ){
	int i;

	int d_num = 0;
	int u_num = 0;

	double step_d = ( end_d - start_d ) / double(rep_d);
	double step_u = ( end_u - start_u ) / double(rep_u);

	double d;
	double u;

	double P = 1200.0;
	double h = 300.0;

	int dd;
	int uu;

	refdata.open( "refdata_vu.csv" , ios::out );
	if( refdata.fail() ){
		std::cout << "file open error" << endl;
		exit(1);    // �t�@�C�����J���܂���ł����B
	}

	refdata << std::setprecision(20);


	double P_[3];
	double D_[3];

	uu = 0;
	for( u = start_u ; uu < rep_u ; u += step_u , uu++ ){
		dd = 0;
		for( d = start_d ; dd < rep_d ; d += step_d  ,dd++ ){

			P = 10.0;
			while(1){

				if( !state_ph2(P , u + P/d ) ){
					P = 10000.0;
					break;
				}
				if( d < Rc.rho ){
					break;
				}
				P += 10.0;
			}
	

			if( P < 9000 ){

				P_[0] = P-20;
				P_[2] = P+10;	
				P_[1] = ( P_[0] + P_[2] ) / 2.0 ;

				for( i = 0 ; i < 50 ; i++ ){

					state_ph2( P_[0] , u + P_[0] / d );
					D_[0] = Rc.rho;
					state_ph2( P_[1] , u + P_[1] / d );
					D_[1] = Rc.rho;
					state_ph2( P_[2] , u + P_[2] / d );
					D_[2] = Rc.rho;

					if( d > D_[1] ){
						P_[0] = P_[1];
					}else{
						P_[2] = P_[1];
					}
					P_[1] = ( P_[0] + P_[2] ) / 2.0 ;

				}
				Rc.u = u;

			}else{

				Rc.T = 0.0;
				Rc.P = 0.0;
				Rc.rho = d;
				Rc.h = 0.0;
				Rc.s = 0.0;
				Rc.x = 0.0;
				Rc.cp = 0.0;
				Rc.thc = 0.0;
				Rc.visc = 0.0;
				Rc.st = 0.0;
				Rc.cv = 0.0;
				Rc.u = u;

			}
			data_table_du[uu][dd] = Rc;


			refdata << data_table_du[uu][dd].T << ",";
			refdata << data_table_du[uu][dd].P << ",";
			refdata << data_table_du[uu][dd].rho << ",";
			refdata << data_table_du[uu][dd].h << ",";
			refdata << data_table_du[uu][dd].s << ",";
			refdata << data_table_du[uu][dd].x << ",";
			refdata << data_table_du[uu][dd].cp << ",";
			refdata << data_table_du[uu][dd].thc << ",";
			refdata << data_table_du[uu][dd].visc << ",";
			refdata << data_table_du[uu][dd].st << ",";
			refdata << data_table_du[uu][dd].cv << ",";
			refdata << data_table_du[uu][dd].u << ",";
			refdata << endl;
		}
		std::cout << "\r" << (double)uu/(double)rep_u ;
	}

	std::cout << std::endl;
	std::cout << "�e�[�u���쐬�I��" << std::endl;


}
/*!
*******************************************************************************
 state_ph2 -- ����Ԃ�p�����v�Z���@->state_ph�̌v�Z���xUp
 ����� ���
*******************************************************************************
*/
bool DXProperty::state_ph2(double pressure, double enthalpy)
{

	if( enthalpy < data_table[0][0].h || data_table[rep_h-1][0].h < enthalpy ||
		pressure < data_table[0][0].P || data_table[0][rep_p-1].P < pressure ){
//			cout<< name.str() << " p = " << pressure << " h = " << enthalpy << endl;
			if( _DLLsetup == 1 ){
				state_ph( pressure , enthalpy );
			}
			return false;
	}

	DXFluid ref;
	int i;
	int j;
	for( i = 0 ; i < rep_h ; i++ ){
		if( enthalpy < data_table[i][0].h ){
			break;
		}
	}
	for( j = 0 ; j < rep_p ; j++ ){
		if( pressure < data_table[0][j].P ){
			break;
		}
	}


	double h1;
	double h2;
	double p1;
	double p2;

	h2 = ( enthalpy - data_table[i-1][j].h ) / ( data_table[i][j].h - data_table[i-1][j].h );
	h1 = 1.0 - h2;
	p2 = ( pressure - data_table[i][j-1].P ) / ( data_table[i][j].P - data_table[i][j-1].P );
	p1 = 1.0 - p2;

	double a11,a12,a21,a22;
	a11 = h1 * p1;
	a12 = h1 * p2;
	a21 = h2 * p1;
	a22 = h2 * p2;

	ref.T    = data_table[i-1][j-1].T * a11 + data_table[i-1][j].T * a12 + data_table[i][j-1].T * a21 + data_table[i][j].T * a22 ;
	ref.P    = pressure;
//	ref.rho  = data_table[i-1][j-1].rho*a11 + data_table[i-1][j].rho*a12 + data_table[i][j-1].rho*a21 + data_table[i][j].rho*a22 ;
	ref.rho  = 1.0/data_table[i-1][j-1].rho*a11 + 1.0/data_table[i-1][j].rho*a12 + 1.0/data_table[i][j-1].rho*a21 + 1.0/data_table[i][j].rho*a22 ;
	ref.rho  = 1.0 / ref.rho;
	ref.h    = enthalpy;
	ref.s    = data_table[i-1][j-1].s * a11 + data_table[i-1][j].s * a12 + data_table[i][j-1].s * a21 + data_table[i][j].s * a22 ;
	ref.x    = data_table[i-1][j-1].x * a11 + data_table[i-1][j].x * a12 + data_table[i][j-1].x * a21 + data_table[i][j].x * a22 ;
	ref.cp   = data_table[i-1][j-1].cp * a11 + data_table[i-1][j].cp * a12 + data_table[i][j-1].cp * a21 + data_table[i][j].cp * a22 ;
	ref.thc  = data_table[i-1][j-1].thc * a11 + data_table[i-1][j].thc * a12 + data_table[i][j-1].thc * a21 + data_table[i][j].thc * a22 ;
	ref.visc = data_table[i-1][j-1].visc * a11 + data_table[i-1][j].visc * a12 + data_table[i][j-1].visc * a21 + data_table[i][j].visc * a22 ;
	ref.st   = data_table[i-1][j-1].st * a11 + data_table[i-1][j].st * a12 + data_table[i][j-1].st * a21 + data_table[i][j].st * a22 ;
	ref.cv   = data_table[i-1][j-1].cv * a11 + data_table[i-1][j].cv * a12 + data_table[i][j-1].cv * a21 + data_table[i][j].cv * a22 ;


	Rc = ref;

	return true;



}
/*!
*******************************************************************************
 state_ph3    ---   ���O�����W����Ԗ@��p�����v�Z���@->state_ph2�̌v�Z���xUp
 ����� ����
*******************************************************************************
*/
bool DXProperty::state_ph3(double pressure, double enthalpy)
{
	if( enthalpy < data_table[0][0].h || data_table[rep_h-2][0].h < enthalpy ||
		pressure < data_table[0][0].P || data_table[0][rep_p-2].P < pressure ){
	//		cout<< name.str() << " p = " << pressure << " h = " << enthalpy << endl;
			if( _DLLsetup == 1 ){
				state_ph( pressure , enthalpy );
			}
			return false;
	}

	DXFluid ref;
	int i;
	int j;
	for( i = 0 ; i < rep_h ; i++ ){
		if( enthalpy < data_table[i][0].h ){
			break;
		}
	}
	for( j = 0 ; j < rep_p ; j++ ){
		if( pressure < data_table[0][j].P ){
			break;
		}
	}

	double h1;
	double h2;
	double h3;
	double p1;
	double p2;
	double p3;

	h1 = ( enthalpy - data_table[i][j].h ) * ( enthalpy - data_table[i+1][j].h ) / ( data_table[i-1][j].h - data_table[i][j].h ) / ( data_table[i-1][j].h - data_table[i+1][j].h );
	h2 = ( enthalpy - data_table[i-1][j].h ) * ( enthalpy - data_table[i+1][j].h ) / ( data_table[i][j].h - data_table[i-1][j].h ) / ( data_table[i][j].h - data_table[i+1][j].h );
	h3 = ( enthalpy - data_table[i-1][j].h ) * ( enthalpy - data_table[i][j].h ) / ( data_table[i+1][j].h - data_table[i-1][j].h ) / ( data_table[i+1][j].h - data_table[i][j].h );

	p1 = ( pressure - data_table[i][j].P ) * ( pressure - data_table[i][j+1].P ) / ( data_table[i][j-1].P - data_table[i][j].P ) / ( data_table[i][j-1].P - data_table[i][j+1].P );
	p2 = ( pressure - data_table[i][j-1].P ) * ( pressure - data_table[i][j+1].P ) / ( data_table[i][j].P - data_table[i][j-1].P ) / ( data_table[i][j].P - data_table[i][j+1].P );
	p3 = ( pressure - data_table[i][j-1].P ) * ( pressure - data_table[i][j].P ) / ( data_table[i][j+1].P - data_table[i][j-1].P ) / ( data_table[i][j+1].P - data_table[i][j].P );


	double a11,a12,a13,a21,a22,a23,a31,a32,a33;
	a11 = h1 * p1;
	a12 = h1 * p2;
	a13 = h1 * p3;
	a21 = h2 * p1;
	a22 = h2 * p2;
	a23 = h2 * p3;
	a31 = h3 * p1;
	a32 = h3 * p2;
	a33 = h3 * p3;

	ref.T    =   data_table[i-1][j-1].T * a11 + data_table[i-1][j].T * a12 + data_table[i-1][j+1].T * a13
		       + data_table[i][j-1].T   * a21 + data_table[i][j].T   * a22 + data_table[i][j+1].T   * a23
			   + data_table[i+1][j-1].T * a31 + data_table[i+1][j].T * a32 + data_table[i+1][j+1].T * a33;
	ref.P    = pressure;
/*	ref.rho  =   data_table[i-1][j-1].rho * a11 + data_table[i-1][j].rho * a12 + data_table[i-1][j+1].rho * a13
		       + data_table[i][j-1].rho   * a21 + data_table[i][j].rho   * a22 + data_table[i][j+1].rho   * a23
			   + data_table[i+1][j-1].rho * a31 + data_table[i+1][j].rho * a32 + data_table[i+1][j+1].rho * a33;*/
	ref.rho  =   1.0/data_table[i-1][j-1].rho * a11 + 1.0/data_table[i-1][j].rho * a12 + 1.0/data_table[i-1][j+1].rho * a13
		       + 1.0/data_table[i][j-1].rho   * a21 + 1.0/data_table[i][j].rho   * a22 + 1.0/data_table[i][j+1].rho   * a23
			   + 1.0/data_table[i+1][j-1].rho * a31 + 1.0/data_table[i+1][j].rho * a32 + 1.0/data_table[i+1][j+1].rho * a33;
	ref.rho  = 1.0 / ref.rho;
	ref.h    = enthalpy;
	ref.s    =   data_table[i-1][j-1].s * a11 + data_table[i-1][j].s * a12 + data_table[i-1][j+1].s * a13
		       + data_table[i][j-1].s   * a21 + data_table[i][j].s   * a22 + data_table[i][j+1].s   * a23
			   + data_table[i+1][j-1].s * a31 + data_table[i+1][j].s * a32 + data_table[i+1][j+1].s * a33;
	ref.x    =   data_table[i-1][j-1].x * a11 + data_table[i-1][j].x * a12 + data_table[i-1][j+1].x * a13
		       + data_table[i][j-1].x   * a21 + data_table[i][j].x   * a22 + data_table[i][j+1].x   * a23
			   + data_table[i+1][j-1].x * a31 + data_table[i+1][j].x * a32 + data_table[i+1][j+1].x * a33;
	ref.cp   = data_table[i-1][j-1].cp * a11 + data_table[i-1][j].cp * a12 + data_table[i-1][j+1].cp * a13
		       + data_table[i][j-1].cp   * a21 + data_table[i][j].cp   * a22 + data_table[i][j+1].cp   * a23
			   + data_table[i+1][j-1].cp * a31 + data_table[i+1][j].cp * a32 + data_table[i+1][j+1].cp * a33;
	ref.thc  = data_table[i-1][j-1].thc * a11 + data_table[i-1][j].thc * a12 + data_table[i-1][j+1].thc * a13
		       + data_table[i][j-1].thc   * a21 + data_table[i][j].thc   * a22 + data_table[i][j+1].thc   * a23
			   + data_table[i+1][j-1].thc * a31 + data_table[i+1][j].thc * a32 + data_table[i+1][j+1].thc * a33;
	ref.visc = data_table[i-1][j-1].visc * a11 + data_table[i-1][j].visc * a12 + data_table[i-1][j+1].visc * a13
		       + data_table[i][j-1].visc   * a21 + data_table[i][j].visc   * a22 + data_table[i][j+1].visc   * a23
			   + data_table[i+1][j-1].visc * a31 + data_table[i+1][j].visc * a32 + data_table[i+1][j+1].visc * a33;
	ref.st   = data_table[i-1][j-1].st * a11 + data_table[i-1][j].st * a12 + data_table[i-1][j+1].st * a13
		       + data_table[i][j-1].st   * a21 + data_table[i][j].st   * a22 + data_table[i][j+1].st   * a23
			   + data_table[i+1][j-1].st * a31 + data_table[i+1][j].st * a32 + data_table[i+1][j+1].st * a33;
	ref.cv   = data_table[i-1][j-1].cv * a11 + data_table[i-1][j].cv * a12 + data_table[i-1][j+1].cv * a13
		       + data_table[i][j-1].cv   * a21 + data_table[i][j].cv   * a22 + data_table[i][j+1].cv   * a23
			   + data_table[i+1][j-1].cv * a31 + data_table[i+1][j].cv * a32 + data_table[i+1][j+1].cv * a33;

	Rc = ref;

	return true;
}
/*!
*******************************************************************************
 state_tp2 -- ����Ԃ�p�����v�Z���@->state_tp�̌v�Z���xUp
 ����� ���
*******************************************************************************
*/

bool DXProperty::state_tp2(double temperature, double pressure)
{
	if( temperature < data_table_tp[0][0].T || data_table_tp[rep_h-1][0].T < temperature ||
		pressure < data_table_tp[0][0].P || data_table_tp[0][rep_p-1].P < pressure ){
//			cout<< name.str() << " p = " << pressure << " h = " << enthalpy << endl;
			if( _DLLsetup == 1 ){
				state_tp( temperature , pressure );
			}
			return false;
	}

	DXFluid ref;
	int i;
	int j;
	for( i = 0 ; i < rep_t ; i++ ){
		if( temperature < data_table_tp[i][0].T ){
			break;
		}
	}
	for( j = 0 ; j < rep_p ; j++ ){
		if( pressure < data_table_tp[0][j].P ){
			break;
		}
	}


	double t1;
	double t2;
	double p1;
	double p2;

	t2 = ( temperature - data_table_tp[i-1][j].T ) / ( data_table_tp[i][j].T - data_table_tp[i-1][j].T );
	t1 = 1.0 - t2;
	p2 = ( pressure - data_table_tp[i][j-1].P ) / ( data_table_tp[i][j].P - data_table_tp[i][j-1].P );
	p1 = 1.0 - p2;

	double a11,a12,a21,a22;
	a11 = t1 * p1;
	a12 = t1 * p2;
	a21 = t2 * p1;
	a22 = t2 * p2;

	ref.T    = temperature;
	ref.P    = pressure;
	ref.rho  = data_table_tp[i-1][j-1].rho*a11 + data_table_tp[i-1][j].rho*a12 + data_table_tp[i][j-1].rho*a21 + data_table_tp[i][j].rho*a22 ;
	ref.h    = data_table_tp[i-1][j-1].h * a11 + data_table_tp[i-1][j].h * a12 + data_table_tp[i][j-1].h * a21 + data_table_tp[i][j].h * a22 ;
	ref.s    = data_table_tp[i-1][j-1].s * a11 + data_table_tp[i-1][j].s * a12 + data_table_tp[i][j-1].s * a21 + data_table_tp[i][j].s * a22 ;
	ref.x    = data_table_tp[i-1][j-1].x * a11 + data_table_tp[i-1][j].x * a12 + data_table_tp[i][j-1].x * a21 + data_table_tp[i][j].x * a22 ;
	ref.cp   = data_table_tp[i-1][j-1].cp * a11 + data_table_tp[i-1][j].cp * a12 + data_table_tp[i][j-1].cp * a21 + data_table_tp[i][j].cp * a22 ;
	ref.thc  = data_table_tp[i-1][j-1].thc * a11 + data_table_tp[i-1][j].thc * a12 + data_table_tp[i][j-1].thc * a21 + data_table_tp[i][j].thc * a22 ;
	ref.visc = data_table_tp[i-1][j-1].visc * a11 + data_table_tp[i-1][j].visc * a12 + data_table_tp[i][j-1].visc * a21 + data_table_tp[i][j].visc * a22 ;
	ref.st   = data_table_tp[i-1][j-1].st * a11 + data_table_tp[i-1][j].st * a12 + data_table_tp[i][j-1].st * a21 + data_table_tp[i][j].st * a22 ;
	ref.cv   = data_table_tp[i-1][j-1].cv * a11 + data_table_tp[i-1][j].cv * a12 + data_table_tp[i][j-1].cv * a21 + data_table_tp[i][j].cv * a22 ;


	Rc = ref;

	return true;



}

/*!
*******************************************************************************
 state_tp3 -- ���O�����W����Ԗ@��p�����v�Z���@->state_tp2�̌v�Z���xUp
 ����� ����
*******************************************************************************
*/

bool DXProperty::state_tp3(double temperature, double pressure)
{
	if( temperature < data_table_tp[0][0].T || data_table_tp[rep_h-2][0].T < temperature ||
		pressure < data_table_tp[0][0].P || data_table_tp[0][rep_p-2].P < pressure ){
//			cout<< name.str() << " p = " << pressure << " T = " << temperature << endl;
			if( _DLLsetup == 1 ){
				state_tp( temperature , pressure );
			}
			return false;
	}

	DXFluid ref;
	int i;
	int j;
	for( i = 0 ; i < rep_t ; i++ ){
		if( temperature < data_table_tp[i][0].T ){
			break;
		}
	}
	for( j = 0 ; j < rep_p ; j++ ){
		if( pressure < data_table_tp[0][j].P ){
			break;
		}
	}


	double t1;
	double t2;
	double t3;
	double p1;
	double p2;
	double p3;

	t1 = ( temperature - data_table_tp[i][j].T ) * ( temperature - data_table_tp[i+1][j].T ) / ( data_table_tp[i-1][j].T - data_table_tp[i][j].T ) / ( data_table_tp[i-1][j].T - data_table_tp[i+1][j].T );
	t2 = ( temperature - data_table_tp[i-1][j].T ) * ( temperature - data_table_tp[i+1][j].T ) / ( data_table_tp[i][j].T - data_table_tp[i-1][j].T ) / ( data_table_tp[i][j].T - data_table_tp[i+1][j].T );
	t3 = ( temperature - data_table_tp[i-1][j].T ) * ( temperature - data_table_tp[i][j].T ) / ( data_table_tp[i+1][j].T - data_table_tp[i-1][j].T ) / ( data_table_tp[i+1][j].T - data_table_tp[i][j].T );

	p1 = ( pressure - data_table_tp[i][j].P ) * ( pressure - data_table_tp[i][j+1].P ) / ( data_table_tp[i][j-1].P - data_table_tp[i][j].P ) / ( data_table_tp[i][j-1].P - data_table_tp[i][j+1].P );
	p2 = ( pressure - data_table_tp[i][j-1].P ) * ( pressure - data_table_tp[i][j+1].P ) / ( data_table_tp[i][j].P - data_table_tp[i][j-1].P ) / ( data_table_tp[i][j].P - data_table_tp[i][j+1].P );
	p3 = ( pressure - data_table_tp[i][j-1].P ) * ( pressure - data_table_tp[i][j].P ) / ( data_table_tp[i][j+1].P - data_table_tp[i][j-1].P ) / ( data_table_tp[i][j+1].P - data_table_tp[i][j].P );


	double a11,a12,a13,a21,a22,a23,a31,a32,a33;
	a11 = t1 * p1;
	a12 = t1 * p2;
	a13 = t1 * p3;
	a21 = t2 * p1;
	a22 = t2 * p2;
	a23 = t2 * p3;
	a31 = t3 * p1;
	a32 = t3 * p2;
	a33 = t3 * p3;

	ref.T    = temperature;
	ref.P    = pressure;
	ref.rho  =   data_table_tp[i-1][j-1].rho * a11 + data_table_tp[i-1][j].rho * a12 + data_table_tp[i-1][j+1].rho * a13
		       + data_table_tp[i][j-1].rho   * a21 + data_table_tp[i][j].rho   * a22 + data_table_tp[i][j+1].rho   * a23
			   + data_table_tp[i+1][j-1].rho * a31 + data_table_tp[i+1][j].rho * a32 + data_table_tp[i+1][j+1].rho * a33;
	ref.h    =   data_table_tp[i-1][j-1].h * a11 + data_table_tp[i-1][j].h * a12 + data_table_tp[i-1][j+1].h * a13
		       + data_table_tp[i][j-1].h   * a21 + data_table_tp[i][j].h   * a22 + data_table_tp[i][j+1].h   * a23
			   + data_table_tp[i+1][j-1].h * a31 + data_table_tp[i+1][j].h * a32 + data_table_tp[i+1][j+1].h * a33;
	ref.s    =   data_table_tp[i-1][j-1].s * a11 + data_table_tp[i-1][j].s * a12 + data_table_tp[i-1][j+1].s * a13
		       + data_table_tp[i][j-1].s   * a21 + data_table_tp[i][j].s   * a22 + data_table_tp[i][j+1].s   * a23
			   + data_table_tp[i+1][j-1].s * a31 + data_table_tp[i+1][j].s * a32 + data_table_tp[i+1][j+1].s * a33;
	ref.x    =   data_table_tp[i-1][j-1].x * a11 + data_table_tp[i-1][j].x * a12 + data_table_tp[i-1][j+1].x * a13
		       + data_table_tp[i][j-1].x   * a21 + data_table_tp[i][j].x   * a22 + data_table_tp[i][j+1].x   * a23
			   + data_table_tp[i+1][j-1].x * a31 + data_table_tp[i+1][j].x * a32 + data_table_tp[i+1][j+1].x * a33;
	ref.cp   = data_table_tp[i-1][j-1].cp * a11 + data_table_tp[i-1][j].cp * a12 + data_table_tp[i-1][j+1].cp * a13
		       + data_table_tp[i][j-1].cp   * a21 + data_table_tp[i][j].cp   * a22 + data_table_tp[i][j+1].cp   * a23
			   + data_table_tp[i+1][j-1].cp * a31 + data_table_tp[i+1][j].cp * a32 + data_table_tp[i+1][j+1].cp * a33;
	ref.thc  = data_table_tp[i-1][j-1].thc * a11 + data_table_tp[i-1][j].thc * a12 + data_table_tp[i-1][j+1].thc * a13
		       + data_table_tp[i][j-1].thc   * a21 + data_table_tp[i][j].thc   * a22 + data_table_tp[i][j+1].thc   * a23
			   + data_table_tp[i+1][j-1].thc * a31 + data_table_tp[i+1][j].thc * a32 + data_table_tp[i+1][j+1].thc * a33;
	ref.visc = data_table_tp[i-1][j-1].visc * a11 + data_table_tp[i-1][j].visc * a12 + data_table_tp[i-1][j+1].visc * a13
		       + data_table_tp[i][j-1].visc   * a21 + data_table_tp[i][j].visc   * a22 + data_table_tp[i][j+1].visc   * a23
			   + data_table_tp[i+1][j-1].visc * a31 + data_table_tp[i+1][j].visc * a32 + data_table_tp[i+1][j+1].visc * a33;
	ref.st   = data_table_tp[i-1][j-1].st * a11 + data_table_tp[i-1][j].st * a12 + data_table_tp[i-1][j+1].st * a13
		       + data_table_tp[i][j-1].st   * a21 + data_table_tp[i][j].st   * a22 + data_table_tp[i][j+1].st   * a23
			   + data_table_tp[i+1][j-1].st * a31 + data_table_tp[i+1][j].st * a32 + data_table_tp[i+1][j+1].st * a33;
	ref.cv   = data_table_tp[i-1][j-1].cv * a11 + data_table_tp[i-1][j].cv * a12 + data_table_tp[i-1][j+1].cv * a13
		       + data_table_tp[i][j-1].cv   * a21 + data_table_tp[i][j].cv   * a22 + data_table_tp[i][j+1].cv   * a23
			   + data_table_tp[i+1][j-1].cv * a31 + data_table_tp[i+1][j].cv * a32 + data_table_tp[i+1][j+1].cv * a33;

	Rc = ref;

	return true;
}
/*!
*******************************************************************************
 sat_p2 -- ����Ԃ�p�����v�Z���@->sat_p�̌v�Z���xUp
 ����� ���
*******************************************************************************
*/
bool DXProperty::sat_p2(double pressure)
{


	if(	pressure < data_table_p[0][rep_p-1].P || data_table_p[0][0].P < pressure ){
//			cout<< "p = " << pressure << endl;
			if( _DLLsetup == 1 ){
				sat_p( pressure );
			}
			return false;
	}


	DXFluid refl,refv;

	int j;
	for( j = 0 ; j < rep_p ; j++ ){
		if( pressure > data_table_p[0][j].P ){
			break;
		}
	}

	double p1;
	double p2;

	p2 = ( pressure - data_table_p[0][j-1].P ) / ( data_table_p[0][j].P - data_table_p[0][j-1].P );
	p1 = 1.0 - p2;
 
	refl.T    = data_table_p[0][j-1].T    * p1 + data_table_p[0][j].T    * p2;
	refl.P    = data_table_p[0][j-1].P    * p1 + data_table_p[0][j].P    * p2;
	refl.rho  = data_table_p[0][j-1].rho  * p1 + data_table_p[0][j].rho  * p2;
	refl.h    = data_table_p[0][j-1].h    * p1 + data_table_p[0][j].h    * p2;
	refl.s    = data_table_p[0][j-1].s    * p1 + data_table_p[0][j].s    * p2;
	refl.x    = data_table_p[0][j-1].x    * p1 + data_table_p[0][j].x    * p2;
	refl.cp   = data_table_p[0][j-1].cp   * p1 + data_table_p[0][j].cp   * p2;
	refl.thc  = data_table_p[0][j-1].thc  * p1 + data_table_p[0][j].thc  * p2;
	refl.visc = data_table_p[0][j-1].visc * p1 + data_table_p[0][j].visc * p2;
	refl.st   = data_table_p[0][j-1].st   * p1 + data_table_p[0][j].st   * p2;
	refl.cv   = data_table_p[0][j-1].cv   * p1 + data_table_p[0][j].cv   * p2;


	refv.T    = data_table_p[1][j-1].T    * p1 + data_table_p[1][j].T    * p2;
	refv.P    = data_table_p[1][j-1].P    * p1 + data_table_p[1][j].P    * p2;
	refv.rho  = data_table_p[1][j-1].rho  * p1 + data_table_p[1][j].rho  * p2;
	refv.h    = data_table_p[1][j-1].h    * p1 + data_table_p[1][j].h    * p2;
	refv.s    = data_table_p[1][j-1].s    * p1 + data_table_p[1][j].s    * p2;
	refv.x    = data_table_p[1][j-1].x    * p1 + data_table_p[1][j].x    * p2;
	refv.cp   = data_table_p[1][j-1].cp   * p1 + data_table_p[1][j].cp   * p2;
	refv.thc  = data_table_p[1][j-1].thc  * p1 + data_table_p[1][j].thc  * p2;
	refv.visc = data_table_p[1][j-1].visc * p1 + data_table_p[1][j].visc * p2;
	refv.st   = data_table_p[1][j-1].st   * p1 + data_table_p[1][j].st   * p2;
	refv.cv   = data_table_p[1][j-1].cv   * p1 + data_table_p[1][j].cv   * p2;

	Rl = refl;
	Rv = refv;


	return true;


}
/*!
*******************************************************************************
 sat_p3 -- ���O�����W����Ԗ@��p�����v�Z���@->sat_p2�̌v�Z���xUp
 ����� ����
*******************************************************************************
*/
bool DXProperty::sat_p3(double pressure)
{


	if(	pressure < data_table_p[0][rep_p-1].P || data_table_p[0][0].P < pressure ){
//			cout<< "p = " << pressure << endl;
			if( _DLLsetup == 1 ){
				sat_p( pressure );
			}
			return false;
	}


	DXFluid refl,refv;

	int j;
	for( j = 0 ; j < rep_p ; j++ ){
		if( pressure > data_table_p[0][j].P ){
			break;
		}
	}

	double p1;
	double p2;
	double p3;

	p1 = ( pressure - data_table_p[0][j-1].P ) * ( pressure - data_table_p[0][j].P ) / ( data_table_p[0][j-2].P - data_table_p[0][j-1].P ) / ( data_table_p[0][j-2].P - data_table_p[0][j].P );
	p2 = ( pressure - data_table_p[0][j-2].P ) * ( pressure - data_table_p[0][j].P ) / ( data_table_p[0][j-1].P - data_table_p[0][j-2].P ) / ( data_table_p[0][j-1].P - data_table_p[0][j].P );
	p3 = ( pressure - data_table_p[0][j-2].P ) * ( pressure - data_table_p[0][j-1].P ) / ( data_table_p[0][j].P - data_table_p[0][j-2].P ) / ( data_table_p[0][j].P - data_table_p[0][j-1].P );

	refl.T    = data_table_p[0][j-2].T    * p1 + data_table_p[0][j-1].T    * p2 + data_table_p[0][j].T    * p3;
	refl.P    = data_table_p[0][j-2].P    * p1 + data_table_p[0][j-1].P    * p2 + data_table_p[0][j].P    * p3;
	refl.rho  = data_table_p[0][j-2].rho  * p1 + data_table_p[0][j-1].rho  * p2 + data_table_p[0][j].rho  * p3;
	refl.h    = data_table_p[0][j-2].h    * p1 + data_table_p[0][j-1].h    * p2 + data_table_p[0][j].h    * p3;
	refl.s    = data_table_p[0][j-2].s    * p1 + data_table_p[0][j-1].s    * p2 + data_table_p[0][j].s    * p3;
	refl.x    = data_table_p[0][j-2].x    * p1 + data_table_p[0][j-1].x    * p2 + data_table_p[0][j].x    * p3;
	refl.cp   = data_table_p[0][j-2].cp   * p1 + data_table_p[0][j-1].cp   * p2 + data_table_p[0][j].cp   * p3;
	refl.thc  = data_table_p[0][j-2].thc  * p1 + data_table_p[0][j-1].thc  * p2 + data_table_p[0][j].thc  * p3;
	refl.visc = data_table_p[0][j-2].visc * p1 + data_table_p[0][j-1].visc * p2 + data_table_p[0][j].visc * p3;
	refl.st   = data_table_p[0][j-2].st   * p1 + data_table_p[0][j-1].st   * p2 + data_table_p[0][j].st   * p3;
	refl.cv   = data_table_p[0][j-2].cv   * p1 + data_table_p[0][j-1].cv   * p2 + data_table_p[0][j].cv   * p3;


	refv.T    = data_table_p[1][j-2].T    * p1 + data_table_p[1][j-1].T    * p2 + data_table_p[1][j].T    * p3;
	refv.P    = data_table_p[1][j-2].P    * p1 + data_table_p[1][j-1].P    * p2 + data_table_p[1][j].P    * p3;
	refv.rho  = data_table_p[1][j-2].rho  * p1 + data_table_p[1][j-1].rho  * p2 + data_table_p[1][j].rho  * p3;
	refv.h    = data_table_p[1][j-2].h    * p1 + data_table_p[1][j-1].h    * p2 + data_table_p[1][j].h    * p3;
	refv.s    = data_table_p[1][j-2].s    * p1 + data_table_p[1][j-1].s    * p2 + data_table_p[1][j].s    * p3;
	refv.x    = data_table_p[1][j-2].x    * p1 + data_table_p[1][j-1].x    * p2 + data_table_p[1][j].x    * p3;
	refv.cp   = data_table_p[1][j-2].cp   * p1 + data_table_p[1][j-1].cp   * p2 + data_table_p[1][j].cp   * p3;
	refv.thc  = data_table_p[1][j-2].thc  * p1 + data_table_p[1][j-1].thc  * p2 + data_table_p[1][j].thc  * p3;
	refv.visc = data_table_p[1][j-2].visc * p1 + data_table_p[1][j-1].visc * p2 + data_table_p[1][j].visc * p3;
	refv.st   = data_table_p[1][j-2].st   * p1 + data_table_p[1][j-1].st   * p2 + data_table_p[1][j].st   * p3;
	refv.cv   = data_table_p[1][j-2].cv   * p1 + data_table_p[1][j-1].cv   * p2 + data_table_p[1][j].cv   * p3;

	Rl = refl;
	Rv = refv;


	return true;


}
/*!
*******************************************************************************
 sat_t2 -- ����Ԃ�p�����v�Z���@->sat_t�̌v�Z���xUp
 ����� ���
*******************************************************************************
*/
bool DXProperty::sat_t2(double temperature)
{


	if(	temperature < data_table_t[0][rep_t-2].T || data_table_t[0][0].T < temperature ){
//			cout<< "p = " << pressure << endl;
			if( _DLLsetup == 1 ){
				sat_t( temperature );
			}
			return false;
	}


	DXFluid refl,refv;

	int j;
	for( j = 0 ; j < rep_t ; j++ ){
		if( temperature > data_table_t[0][j].T ){
			break;
		}
	}

	double t1;
	double t2;

	t2 = ( temperature - data_table_t[0][j-1].T ) / ( data_table_t[0][j].T - data_table_t[0][j-1].T );
	t1 = 1.0 - t2;
 
	refl.T    = data_table_t[0][j-1].T    * t1 + data_table_t[0][j].T    * t2;
	refl.P    = data_table_t[0][j-1].P    * t1 + data_table_t[0][j].P    * t2;
	refl.rho  = data_table_t[0][j-1].rho  * t1 + data_table_t[0][j].rho  * t2;
	refl.h    = data_table_t[0][j-1].h    * t1 + data_table_t[0][j].h    * t2;
	refl.s    = data_table_t[0][j-1].s    * t1 + data_table_t[0][j].s    * t2;
	refl.x    = data_table_t[0][j-1].x    * t1 + data_table_t[0][j].x    * t2;
	refl.cp   = data_table_t[0][j-1].cp   * t1 + data_table_t[0][j].cp   * t2;
	refl.thc  = data_table_t[0][j-1].thc  * t1 + data_table_t[0][j].thc  * t2;
	refl.visc = data_table_t[0][j-1].visc * t1 + data_table_t[0][j].visc * t2;
	refl.st   = data_table_t[0][j-1].st   * t1 + data_table_t[0][j].st   * t2;
	refl.cv   = data_table_t[0][j-1].cv   * t1 + data_table_t[0][j].cv   * t2;


	refv.T    = data_table_t[1][j-1].T    * t1 + data_table_t[1][j].T    * t2;
	refv.P    = data_table_t[1][j-1].P    * t1 + data_table_t[1][j].P    * t2;
	refv.rho  = data_table_t[1][j-1].rho  * t1 + data_table_t[1][j].rho  * t2;
	refv.h    = data_table_t[1][j-1].h    * t1 + data_table_t[1][j].h    * t2;
	refv.s    = data_table_t[1][j-1].s    * t1 + data_table_t[1][j].s    * t2;
	refv.x    = data_table_t[1][j-1].x    * t1 + data_table_t[1][j].x    * t2;
	refv.cp   = data_table_t[1][j-1].cp   * t1 + data_table_t[1][j].cp   * t2;
	refv.thc  = data_table_t[1][j-1].thc  * t1 + data_table_t[1][j].thc  * t2;
	refv.visc = data_table_t[1][j-1].visc * t1 + data_table_t[1][j].visc * t2;
	refv.st   = data_table_t[1][j-1].st   * t1 + data_table_t[1][j].st   * t2;
	refv.cv   = data_table_t[1][j-1].cv   * t1 + data_table_t[1][j].cv   * t2;

	Rl = refl;
	Rv = refv;


	return true;


}
/*!
*******************************************************************************
  sat_t3 -- ���O�����W����Ԗ@��p�����v�Z���@->sat_t2�̌v�Z���xUp
  ����� ����
*******************************************************************************
*/
bool DXProperty::sat_t3(double temperature)
{


	if(	temperature < data_table_t[0][rep_t-2].T || data_table_t[0][0].T < temperature ){
//			cout<< "t = " << temperature << endl;
			if( _DLLsetup == 1 ){
				sat_t( temperature );
			}
			return false;
	}


	DXFluid refl,refv;

	int j;
	for( j = 0 ; j < rep_t ; j++ ){
		if( temperature > data_table_t[0][j].T ){
			break;
		}
	}

	double t1;
	double t2;
	double t3;

	t1 = ( temperature - data_table_t[0][j].T ) * ( temperature - data_table_t[0][j+1].T ) / ( data_table_t[0][j-1].T - data_table_t[0][j].T ) / ( data_table_t[0][j-1].T - data_table_t[0][j+1].T );
	t2 = ( temperature - data_table_t[0][j-1].T ) * ( temperature - data_table_t[0][j+1].T ) / ( data_table_t[0][j].T - data_table_t[0][j-1].T ) / ( data_table_t[0][j].T - data_table_t[0][j+1].T );
	t3 = ( temperature - data_table_t[0][j-1].T ) * ( temperature - data_table_t[0][j].T ) / ( data_table_t[0][j+1].T - data_table_t[0][j-1].T ) / ( data_table_t[0][j+1].T - data_table_t[0][j].T );

	refl.T    = data_table_t[0][j-1].T    * t1 + data_table_t[0][j].T    * t2 + data_table_t[0][j+1].T    * t3;
	refl.P    = data_table_t[0][j-1].P    * t1 + data_table_t[0][j].P    * t2 + data_table_t[0][j+1].P    * t3;
	refl.rho  = data_table_t[0][j-1].rho  * t1 + data_table_t[0][j].rho  * t2 + data_table_t[0][j+1].rho  * t3;
	refl.h    = data_table_t[0][j-1].h    * t1 + data_table_t[0][j].h    * t2 + data_table_t[0][j+1].h    * t3;
	refl.s    = data_table_t[0][j-1].s    * t1 + data_table_t[0][j].s    * t2 + data_table_t[0][j+1].s    * t3;
	refl.x    = data_table_t[0][j-1].x    * t1 + data_table_t[0][j].x    * t2 + data_table_t[0][j+1].x    * t3;
	refl.cp   = data_table_t[0][j-1].cp   * t1 + data_table_t[0][j].cp   * t2 + data_table_t[0][j+1].cp   * t3;
	refl.thc  = data_table_t[0][j-1].thc  * t1 + data_table_t[0][j].thc  * t2 + data_table_t[0][j+1].thc  * t3;
	refl.visc = data_table_t[0][j-1].visc * t1 + data_table_t[0][j].visc * t2 + data_table_t[0][j+1].visc * t3;
	refl.st   = data_table_t[0][j-1].st   * t1 + data_table_t[0][j].st   * t2 + data_table_t[0][j+1].st   * t3;
	refl.cv   = data_table_t[0][j-1].cv   * t1 + data_table_t[0][j].cv   * t2 + data_table_t[0][j+1].cv   * t3;

	refv.T    = data_table_t[1][j-1].T    * t1 + data_table_t[1][j].T    * t2 + data_table_t[1][j+1].T    * t3;
	refv.P    = data_table_t[1][j-1].P    * t1 + data_table_t[1][j].P    * t2 + data_table_t[1][j+1].P    * t3;
	refv.rho  = data_table_t[1][j-1].rho  * t1 + data_table_t[1][j].rho  * t2 + data_table_t[1][j+1].rho  * t3;
	refv.h    = data_table_t[1][j-1].h    * t1 + data_table_t[1][j].h    * t2 + data_table_t[1][j+1].h    * t3;
	refv.s    = data_table_t[1][j-1].s    * t1 + data_table_t[1][j].s    * t2 + data_table_t[1][j+1].s    * t3;
	refv.x    = data_table_t[1][j-1].x    * t1 + data_table_t[1][j].x    * t2 + data_table_t[1][j+1].x    * t3;
	refv.cp   = data_table_t[1][j-1].cp   * t1 + data_table_t[1][j].cp   * t2 + data_table_t[1][j+1].cp   * t3;
	refv.thc  = data_table_t[1][j-1].thc  * t1 + data_table_t[1][j].thc  * t2 + data_table_t[1][j+1].thc  * t3;
	refv.visc = data_table_t[1][j-1].visc * t1 + data_table_t[1][j].visc * t2 + data_table_t[1][j+1].visc * t3;
	refv.st   = data_table_t[1][j-1].st   * t1 + data_table_t[1][j].st   * t2 + data_table_t[1][j+1].st   * t3;
	refv.cv   = data_table_t[1][j-1].cv   * t1 + data_table_t[1][j].cv   * t2 + data_table_t[1][j+1].cv   * t3;

	Rl = refl;
	Rv = refv;


	return true;
}



bool DXProperty::state_du2(double density, double energy)
{
	if( energy < data_table_du[0][0].h || data_table_du[rep_h-1][0].h < energy ||
		density < data_table_du[0][0].P || data_table_du[0][rep_p-1].P < density ){
//			cout<< name.str() << " p = " << pressure << " h = " << enthalpy << endl;
//			if( _DLLsetup == 1 ){
//				state_ph( pressure , enthalpy );
//			}
			return false;
	}

	DXFluid ref;
	int i;
	int j;
	for( i = 0 ; i < rep_u ; i++ ){
		if( energy < data_table_du[i][0].u ){
			break;
		}
	}
	for( j = 0 ; j < rep_d ; j++ ){
		if( density < data_table_du[0][j].rho ){
			break;
		}
	}


	double u1;
	double u2;
	double d1;
	double d2;

	u2 = ( energy - data_table_du[i-1][j].u ) / ( data_table_du[i][j].u - data_table_du[i-1][j].u );
	u1 = 1.0 - u2;
	d2 = ( density - data_table_du[i][j-1].rho ) / ( data_table_du[i][j].rho - data_table_du[i][j-1].rho );
	d1 = 1.0 - d2;

	double a11,a12,a21,a22;
	a11 = u1 * d1;
	a12 = u1 * d2;
	a21 = u2 * d1;
	a22 = u2 * d2;

	ref.T    = data_table_du[i-1][j-1].T * a11 + data_table_du[i-1][j].T * a12 + data_table_du[i][j-1].T * a21 + data_table_du[i][j].T * a22 ;
	ref.P    = data_table_du[i-1][j-1].P * a11 + data_table_du[i-1][j].P * a12 + data_table_du[i][j-1].P * a21 + data_table_du[i][j].P * a22 ;
	ref.rho  = density;
	ref.h    = data_table_du[i-1][j-1].h * a11 + data_table_du[i-1][j].h * a12 + data_table_du[i][j-1].h * a21 + data_table_du[i][j].h * a22 ;
	ref.s    = data_table_du[i-1][j-1].s * a11 + data_table_du[i-1][j].s * a12 + data_table_du[i][j-1].s * a21 + data_table_du[i][j].s * a22 ;
	ref.x    = data_table_du[i-1][j-1].x * a11 + data_table_du[i-1][j].x * a12 + data_table_du[i][j-1].x * a21 + data_table_du[i][j].x * a22 ;
	ref.cp   = data_table_du[i-1][j-1].cp * a11 + data_table_du[i-1][j].cp * a12 + data_table_du[i][j-1].cp * a21 + data_table_du[i][j].cp * a22 ;
	ref.thc  = data_table_du[i-1][j-1].thc * a11 + data_table_du[i-1][j].thc * a12 + data_table_du[i][j-1].thc * a21 + data_table_du[i][j].thc * a22 ;
	ref.visc = data_table_du[i-1][j-1].visc * a11 + data_table_du[i-1][j].visc * a12 + data_table_du[i][j-1].visc * a21 + data_table_du[i][j].visc * a22 ;
	ref.st   = data_table_du[i-1][j-1].st * a11 + data_table_du[i-1][j].st * a12 + data_table_du[i][j-1].st * a21 + data_table_du[i][j].st * a22 ;
	ref.cv   = data_table_du[i-1][j-1].cv * a11 + data_table_du[i-1][j].cv * a12 + data_table_du[i][j-1].cv * a21 + data_table_du[i][j].cv * a22 ;
	ref.u    = energy;


	Rc = ref;

	return true;



}

void DXProperty::CopyTable( DXProperty &ref ){

	data_table    = ref.data_table;
	data_table_p  = ref.data_table_p;
	data_table_tp = ref.data_table_tp;
	data_table_t  = ref.data_table_t;
	data_table_du = ref.data_table_du;


	min_P = ref.min_P;
	max_P = ref.max_P;
	min_h = ref.min_h;
	max_h = ref.max_h;
	min_T = ref.min_T;
	max_T = ref.max_T;


	return;
}
#ifndef __CNewtonRaphsonMethodPlus_H
#define __CNewtonRaphsonMethodPlus_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <vector>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif



#pragma warning( disable: 4996 )
#pragma warning( disable: 4819 )
#pragma warning( disable: 4018 )

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)






//Newton-Raphson�@�N���X
class CNewtonRaphsonMethodPlus
{
public:





	int Initialize();//�������֐�



	int SetValiable(int,double);//�ϐ����Z�b�g����֐�
	int SetValiable(int,double,double);//�ϐ����Z�b�g����֐�
	int SetValiable(int,double,double,double,double);//�ϐ����Z�b�g����֐�
	int SetValiableName(int, std::string);

	int SetMaxStep(int,double);//max�X�e�b�v����ݒ肷�郂�W���[��
	int SetMinMaxValue(int,double,double);//�ő�l�ƍŒ��l��ݒ肷�郂�W���[��
	
	int GetValiable(int,double&);//�ϐ����󂯎��֐�
	int GetValiableDirect(int,int,double&);//�ϐ����󂯎��֐�

	int SetResult(int,double,double);//�v�Z���ʂ�����֐�
	int SetResultDirect(int,int,double,double);//�v�Z���ʂ�����֐�
	int SetResultMatrix( int , int , double );//�v�Z���ʂ��i�[����}�g���b�N�X�ɒ��ڒl������֐�

	int SetError(double);//����������Z�b�g����֐�
	int SetDelta(double);//�����ړ��͂��Z�b�g����֐�

	int SetAcc(double);//�����W�����Z�b�g����֐�
	int SetAccDirect(int,double);//�����W����ϐ����ƂɃZ�b�g����֐�
	int SetAccAuto(bool _FlagAccAuto, double _AccAutoMin, double _AccAutoMax, double _AccAutoStep);//�����W���������I�ɕύX����

	int SetPrint( bool , bool );//error�̌��ʕ\�������邩�ǂ������Z�b�g����֐�
	int SetEndReport( bool );//�j���[�g���@�I����ɃG���[�l�Ȃǂ̌��ʂ��o�͂��邩�ǂ����̐ݒ�
	int SetMaxLoop( int );//�ő唽���񐔂���͂���֐�

	int SetName(std::string);//���̃j���[�g�����t�\���@�����ʂ��邽�߂̖��O����͂���֐�

	int StatusOut();//�s��̐��Ȃǂ��A�E�g�v�b�g����f�o�b�O�p�̊֐�

	int SetMatrixOut(bool);
	int MatrixOut();//���݂̃}�g���b�N�X���A�E�g�v�b�g����

	int ValueContainerResize( int );

	int MainLoopInit(void);		//���C�����[�v��������
	int MainLoopCheck(void);	//���C�����[�v�p��������
	int MainLoopReinit(void);	//���C�����[�v�ď�������(���C���ƂȂ�֐�)
	int SubLoopInit(void);		//�T�u���[�v��������
	int SubLoopCheck(void);		//�T�u���[�v�p��������
	int SubLoopReinit(void);	//�T�u���[�v�ď�������

	int GaussianElimination(void);

	int Print(void);
	int PrintAll(void);
	int EndReport(void);
	int SetJacobStep( int );

	void irekae(double &x,double &y);	//�s�����ւ��p
	void irekae(int    &x,   int &y);

	void NR_free(void);//2014�N8�����쌴�ǉ��@�j���[�g���@�̍s��v�f�̉��
	//�ϐ��̐��𐔂���J�E���^
	//���[�v�J�E��


	std::vector<std::vector<double>> Variable;//�C���l�̔z��
	std::vector<std::string> VariableName;//�C���l�̔z��
	std::vector<std::vector<double>> VariableMemory;//�C���l�̗���

	std::vector<double> RelError;//���Ό덷
	std::vector<double> AbsError;//��Ό덷

	std::vector<double> DeltaValue;//�ړ��ʂ̐�Βl


	std::vector<std::vector<double>> ResultMatrix;//�v�Z���ʂ̔z��[�s][��+1]
	std::vector<std::vector<double>> JacobMatrix;//���R�r�s��̔z��[�s][��+1]
	std::vector<std::vector<double>> UpperTriangularMatrix;//��O�p�s��[�s][��+1]
	std::vector<std::vector<double>> DiagonalMatrix;//�Ίi�s��[�s][��+1]

	std::vector<double> UpperC;//��O�p�s������Ƃ��̌W��

	std::vector<double> Acc;//�����W��
	std::vector<double> MaxValue;//�ϐ��̍ő�l
	std::vector<double> MinValue;//�ϐ��̍ŏ��l
	std::vector<double> MaxDisplacement;//�ϐ��̍ő�ړ���
	std::vector<double> Step;//�X�e�b�v
	std::vector<bool> StepControl;//�X�e�b�v����bool
	std::vector<bool> MinMaxControl;//MinMax����bool



	int NumVariable;//�ϐ��̐�=���̐�

	int MainLoopCounter;//�j���[�g���@���C�����[�v�̃J�E���^
	int SubLoopCounter;//�j���[�g���@�T�u���[�v�̃J�E���^

	int JacobSkipCounter;
	bool JacobReset;

	std::string Name;




	int MaxLoop;//�ő唽����

	int JacobStep;//


	double MaxError;
	double AveError;

	double Min;//
	double SparseMin;

	double keisu;
	double keisu2;


	double AccAcc;
	double AccMax;

	//AccAuto�֌W
	double AccAutoMax;
	double AccAutoMin;
	double AccAutoStep;
	bool FlagAccAuto;
	double ErrorMax0;
	double ErrorMax1;
	double ErrorDif0;
	double ErrorDif1;


	double Error;//��������
	double Delta;//���l�I�ɕΔ�������Ƃ��̔����ړ���

	bool FlagConbergence;//�������Ă��邩�ǂ����̃t���O.true�Ŏ������Ă���D
	bool FlagPrint;
	bool FlagPrintAll;
	bool FlagEndReport;
	bool FlagJacobSkip;
	bool FlagMatrixOut;//�}�g���b�N�X���A�E�g�v�b�g���邩�ǂ����̃t���O

	time_t TimeStart;
	time_t TimeNow;


	std::fstream MatrixFile;

	CNewtonRaphsonMethodPlus();
	~CNewtonRaphsonMethodPlus();


};
#endif//CNewtonRaphsonMethod

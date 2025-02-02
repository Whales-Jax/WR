#ifndef __CNewtonRaphsonMethod_H
#define __CNewtonRaphsonMethod_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>

#pragma warning( disable: 4996 )
//#include <boost/numeric/ublas/lu.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>


//�K�E�X�̏����@�i���S�s�{�b�e�B���O�j���s�����߂̊֐��D���L�Q�Ƃ̂��ƁD������Ȃ��Ă���肠��܂���


//���ϐ��j���[�g���@���\�b�h�t���[�����\�z���邽�߂̃N���X
class CNewtonRaphsonMethod
{
public:
	std::string name;
	//�ϐ�
	int i_sys,i_error;		//���[�v���Ŏg�p����J�E���^
	int count;				//�v�Z��
	int EscLoop;			//���������Ƃ݂Ȃ��Ĕ����o���v�Z��
	int count_max;			//�ő�v�Z��			
	int count_div;			//��������
	int num;				//�p�����[�^�̐�
	int _num;				//���p�����[�^�̐�
	int ArrayNum;			//�z��̐�
	int Solver;				//�A���������̉�@ 0:���S�s�{�b�g�I�� 1:�����s�{�b�g�I��
	int ErrorSolver;
	int ErrorLoopOver;
	bool FlagInitial;		//�������̊m�F�p
	bool FlagSetup;			//setup�̊m�F�p
	bool check;				//�����m�F�p
	bool check2;			//�����m�F���ꂪ�m�F�ł��������
	bool Bprt;				//�v�����g�\���p�ϐ�
	bool Bprt_sum;			//�v�����g�\���p�ϐ�
	double delta;			//�����ړ����̌W�� x*delta �������ړ��ʂɂȂ�
	double MinimumEPS;		//�ŏ���delta�ړ���
	double pibot;			//�ŏ��s�{�b�g�l
	double err;				//��Ό덷���e�n
	double acc;				//�ɘa�W��
	double err_sum;
	double err_ave;
	double err_max;
	double AccAcc;			//�����x���z�̉����x
	double AccMax;			//�����x���z�̍ő�l
	//�z��
	int *ErrorMethod;			//�G���[���� 0:���Ό덷�@1:��Ό덷
	double *Error;
	double *answer;				//���ݒn�_�Ǝ��̒n�_�Ƃ̍���
	double *set_value;			//�J��Ԃ��v�Z�ɂ����Ĉꎞ�I�Ɏg�p����l�i������T_vap_eva_en�݂����́j
	double *residual;			//���ݒn�_�ł̃G���[�l
	double *MaxDisplacement;	//�ő�ړ���(max displacement)
	double *_MaxDisplacement;	//�ő�ړ��ʂ̉��z��
	double *MaxValue;			//�ő�l(max value)
	double *MinValue;			//�ŏ��l(min value)
	double *assumed_value;		//���ݒn�_�̒l���X�g�b�N�i������T_vap_eva_en_nt�݂����́j
	double *_assumed_value;		//����l�̉��z��
	double *step_value;			//�����ړ��n�_�ł̒l�̃X�g�b�N(T_vap_eva_en_nt)
	double **jacob;				//�Δ����l�i������Jacob[NO][NO]�j
	double **ErrorAbs;			//���ݒn�_����є����ړ������_�ł̃G���[�l(������Error[NO][NO])��Βl
	double **ErrorRel;			//���ݒn�_����є����ړ������_�ł̃G���[�l(������Error[NO][NO])���Βl
	//matAvecX=vecB�Ƃ���vecX�����߂�
	//boost::numeric::ublas::matrix<double> matA;
	//boost::numeric::ublas::vector<double> vecX;
	//boost::numeric::ublas::vector<double> vecB;
	//�֐�
	CNewtonRaphsonMethod();//���p�����[�^�̐��C�ꕔ�������m��
	~CNewtonRaphsonMethod(void);				//���������
	void setup(int n = 1000,double _delta=1.0+1e-6,double _err2=1e-6, double _pibot=-1.0);
	void reset(void);
	void setErrorMethod(int);
	void main_loop_init(void);			//���C�����[�v��������
	bool main_loop_check(void);			//���C�����[�v�p��������
	void main_loop_reinit(void);		//���C�����[�v�ď�������
	void sub_loop_init(void);			//�T�u���[�v��������
	bool sub_loop_check(void);			//�T�u���[�v�p��������
	void sub_loop_reinit(void);			//�T�u���[�v�ď�������
	void initial(int n = 0);			//�p�����[�^������C�������m��
	void prt(double value=-1.0);		//�G���[�\��
	void prt_sum(void);					//�S�̃G���[�\��
	void prt2(double value=-1.0);		//�G���[�\��
	void prt_sum2(void);					//�S�̃G���[�\��
	void spcacc(double maxacc);				//acc��������

	double getValue(int i);				//�C���l��Ԃ��֐�

	void set_assumed_value(int i,double value,double _max = 0.0);//����l���
	void setValue(int i , double value , double _max = 0.0 , double _MinValue = 0.0 , double _MaxValue = 0.0 );//����l���
	void set_error(double,double x=1.0);//�G���[�Z�b�g
	void set_error2(int,double,double,double x=1.0);//�G���[�Z�b�g2
	void setError(int i ,double Value1 , double Value2 , int _ErrorMethod = 0 , double error = -1 );//�G���[�Z�b�g3
	void setAcc(double);
	void setSolver(int);
	void setEscLoop(int);
	void setDeltaError(double _delta=1.0+1e-6,double _err2=1e-6);


	std::fstream matrixFile;
	void matrixOut(void);
	void matrixFileMake(void);

	//�j���[�g���@�p
private:
	int gaussian_elimination(void);		//���S�s�{�b�g�I���i�x��
	int gaussian_elimination_2(void);	//�����s�{�b�g�I���i����
//	int lu_factorize_ublas(void);		//lu�����@���Ԃ񑬂�
	void irekae(double &x,double &y);
	void irekae(int    &x,   int &y);
};
#endif//CNewtonRaphsonMethod

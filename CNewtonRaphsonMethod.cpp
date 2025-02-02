/**************************************************************************************************
*
*	���ϐ��j���[�g���@ �F CNewtonRaphsonMethod 
*
*	1G03 ���c���
*	1G04 ���c�S
*
*	tool_nc_multi_newton_method�����W���[����͂���ыN����~���܂߂��������ɑς����Ȃ��̂�
*	�J�X�^�}�C�Y����multi_newton�̎g��������悭���邽�߂ɂ���ɃJ�X�^�}�C�Y���܂����D
*	�K�E�X�̏����@�����͂킩��Ȃ��̂őS�R�s�y�ł��D
*	[�ύX�_]
*	�E�����R���p�C���Ή��D
*	�E�ϐ��p�ꉻ�C�z��2�������ς݁D
*	�E�֐��Ƃ��ēƗ����Ă���gaussian_elimination,irekae���N���X�̒��Ɋi�[�D
*	�E���z���o���ۂ�2�_�Ԃ̍����ŏ��ł�DELTA-1.0�ɂȂ�悤�ɏC���D
*	�E�ω����C�ő�ω��ʂ̓����D
*	�E�ϐ��̐�������ŃJ�E���g���邱�Ƃŕϐ��𐔂��Ȃ��Ă������悤�ɉ��ǁD
*   �E�K�E�X�̏����@�ɕ����s�{�b�e�B���O��ǉ�
*   �Eboost���g���Ă��܂��D�C���X�g�[�����Ă��������D
*	
*	[�\��]
*	gaussian_elimination����__GEM�̔�}�N����
*
*	[�X�V����]
*	2008/1/15:	��ʂ�̃o�O���I���D
*   2008/11/20: �����s�{�b�g�I���̃K�E�X�̏����@�쐬
*
*   [boost�̃C���X�g�[�����@]
*   �P�Dboost����肷��D���L�t�H���_�ɂ���̂ŏꏊ�̓V�X�J���ɕ����܂��傤
*   �Q�Dc:\lib\boost�ɉ񓚂���D�t�H���_�\���͈ȉ��̂Ƃ���ɂȂ�͂��ł��D
*   c��[lib]������������[boost]��[boost]
*  �@��[windows]�@�@�@�@�@�@�@ ��[doc]
*  �@��[prigram files]�@�@�@�@ ��[libs]
*    �F�@�@�@�@�@�@�@�@�@�@�@�@�F
*   �R�D�u�c�[���v���u�I�v�V�����v���u�v���W�F�N�g����у\�����[�V�����v���uVC++�f�B���N�g���v
*        �� �u�f�B���N�g����\������v���W�F�N�g�v���u�C���N���[�h �t�@�C���v�� c:\lib\boost ��������B 
*   
*	[��{�I�Ȏg�p�@]�i���Z�����j
*


#include <iostream>
#include "CNewtonRaphsonMethod.h"

void main(void){

	//�����l
	double a = 5;
	double b = 5;
	double c = 5;
	
	//a * b = 10
	//b * c = 20
	//c * a = 30
	//�������v���O����
	
	CNewtonRaphsonMethod mnm;		//class��錾���܂��@��1
	mnm.setup(10);					//�Z�b�g�A�b�v�@�����̓j���[�g���@�̕ϐ��ȏ�ɂ��Ă��������@��2:

	mnm.setValue( 0 , a );
	mnm.setValue( 1 , b );
	mnm.setValue( 2 , c );

	//mnm.setAcc(0.8);				//�����x���z�̓��́@�Ȃ��Ă��@��3:

	mnm.initial();					//�v�Z���J�n����O�ɕK�����������Ă�������
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		//���܂��Ȃ�
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	//���܂��Ȃ�

			a = mnm.getValue(0);
			b = mnm.getValue(1);
			c = mnm.getValue(2);

			//�G���[�l
			mnm.setError( 0 , a * b , 10 );
			mnm.setError( 1 , b * c , 20 );
			mnm.setError( 2 , c * a , 30 );

			mnm.prt();				//�G���[�\��
			mnm.prt_sum();			//�G���[�̍��v��\��
		}
	}
	std::cout << " a = " << a << std::endl;
	std::cout << " b = " << b << std::endl;
	std::cout << " c = " << c << std::endl;

}

*	��1:	����mnm�Ƃ������O���قȂ�ΈقȂ������ϐ��j���[�g���@�ƔF������܂��D���ϐ��j���[�g���@�̒��ɑ��ϐ��j���[�g���@��
*			�\�z�������Ƃ��͂����̖��O�����ς������̂��g�̑��ϐ��j���[�g���@�̒��ɃR�s�y����΂n�j�ł��D
*	��2:	setup����ہC�K�v�ɉ����Ĕ����ړ���(������DELTA)�Ǝ����덷(������ERROR)�ƍŏ��s�{�b�g�l�i�K�E�X�̏����@�ɂ���
*			�ē��ʂȈӖ������l�D�ʏ�͂��̒l���C�ɂ���K�v�͂���܂���D�������Ă���l�������̒l���Βl�Őݒ肵�Ă��������j
*			���w��ł��܂��D�w�肵�����ꍇ�̓p�����[�^�̐��̂��ƂɁC�����ړ��ʁC�����덷�C�ŏ��s�{�b�g�l�̏��ŋL�q���܂��D��
*			��
*			mnm.setup(2,1+1e-3,1e-3);     //���̏ꍇ�����ړ���:1+1e-3�C�����덷:1e-3//////////////////////////////////////////////
*			mnm.setup(2,DELTA,ERROR,1e-9);//���̏ꍇ�����ړ���:DELTA�C�����덷:ERROR,�ŏ��s�{�b�g�l:1e-9
*			�ȂǂƂ����l�ɏ����܂��D���w��̏ꍇ�͔����ړ���=1.0+1e-6�C�����덷=1e-6�C�ŏ��s�{�b�g�l=-1.0�ƂȂ�܂��D��L�̎g
*			�p�@�ł̓p�����[�^�̐�=2�C�����ړ���=1.0+1e-6�C�����덷=1e-6�C�ŏ��s�{�b�g�l=1e-9�ƂȂ�܂��D
*			�s�{�b�g���ŏ��s�{�b�g�l��菬�����Ȃ�ƌv�Z�r���ł��I�����܂��D�i�ŏ��s�{�b�g�l���w��̏ꍇ�͓r���I�����邱�Ƃ͂�
*			��܂���j
*	��3:	�����x���z��ݒ肷�邱�ƂŃj���[�g���@�ɂ��ω��ʂ�}���邱�Ƃ��ł��܂��D�Ⴆ�Εω�����0.3�ɂ������ꍇ��
*			mnm.acc = 0.3;
*			�Ƃ��Ă��������D���̏ꍇ�͌��̕ω��ʂ�0.3�{���܂��D�f�t�H���g��1.0�ł��D
*	��4:	mnm.setValue�ɏ����l�������Cmnm.getValue�Œl�̌Ăяo�������܂��D�Ⴆ�Ώ����������P_vap_eva_in�ƋÏk�����
*			��P_vap_con_in�̂Q���p�����[�^�Ƃ��đ��ϐ��j���[�g���@���s���ꍇ�C���̊e�X�̏����l��1.0,5.0�Ƃ���ƁC
*			mnm.setValue( 0 , 1.0 );//P_vap_eva_in
*			mnm.setValue( 1 , 5.0 );//P_vap_con_in
*			mnm.initial();
*				�c
*			P_vap_eva_in = mnm.getValue( 0 );//�l�̑��
*			P_vap_con_in = mnm.getValue( 1 );
*			�ƂȂ�܂��D
*			�����l�ɑ����čő�ω��ʂ����邱�ƂŃj���[�g���@�ɂ��ω���}�����邱�Ƃ��ł��܂��D����͏�q�̕ω������D�悳
*			��܂��D�Ⴆ��
*			mnm.setValue( 0 , 1.0 , 0.05 );//G_vap_eva_in
*			�Ƃ���Ƃ��̃p�����[�^�͍ő�ł�0.05�����ϓ����܂���D
*	��5:	mnm.initial()�ɂ��z��̃T�C�Y�����������܂��D�z���setValue�֐��ɂ��J�E���g����̂ł������O�ɉ���l�����
*			���Ă��������D
*	��6:	���ϐ��j���[�g���@�̃G���[�l��setError�֐��ɂ��ݒ肵�܂��D���R�p�����[�^�̐�����setError�֐����L�q���܂��D�G��
*			�[�l�̐ݒ�̊�{��
*			mnm.setError( �G���[�ԍ� , �ϐ��@, �ϐ��@);
*			�ƂȂ�܂����C���̓�̕ϐ��̒l����v����悤�Ɏ����v�Z���s���܂��D
*	��7:	�G���[�\���͌ʂ̕\���ƑS�̂̕\����2��ނ��g���܂��D�v�f�������Ȃ��Ƃ��Ȃǂ�prt()���C�����Ƃ���prt_sum()�Ƃ��g����
*			���Ă��������D
*
*
*
**************************************************************************************************/
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include "CNewtonRaphsonMethod.h"

double Newton(double X1,double X2,double Y1,double Y2)
{
	return X1-Y1*(X2-X1)/(Y2-Y1);
}

CNewtonRaphsonMethod::CNewtonRaphsonMethod()
{
	FlagSetup = false;

}
CNewtonRaphsonMethod::~CNewtonRaphsonMethod(void)
{

	if( FlagSetup == false ){
		setup(2);
	}

	int i;
	//�z��p���������
	delete [] ErrorMethod;
	delete [] Error;
	delete [] answer;
	delete [] set_value;
	delete [] residual;
	delete [] MaxDisplacement;
	delete [] _MaxDisplacement;
	delete [] assumed_value;
	delete [] _assumed_value;
	delete [] step_value;
	for( i = 0 ; i < ArrayNum+1 ; i++ ){
		delete [] ErrorAbs[i];
		delete [] ErrorRel[i];
	}
	for( i = 0 ; i < ArrayNum ; i++ ){
		delete [] jacob[i];
	}
	delete [] ErrorAbs;
	delete [] ErrorRel;
	delete [] jacob;
}
//����l�̑��
void CNewtonRaphsonMethod::setup(int n,double _delta,double _err,double _pibot)
{

	FlagSetup = true;

	//���ϐ��̐��C�����ړ����C���e�G���[�l
	_num  = 0;
	delta = _delta;
	err	  = _err;
	pibot = _pibot;
	num   = 0;
	ArrayNum = n;
	FlagInitial  = false;
	acc    = 1.0;
	count_div = 10;
	count_max = 10000;
	EscLoop = 30;
	Solver = 0;
	err_sum = 100.0;
	err_max = 100000;
	ErrorLoopOver = 1;
	MinimumEPS = _delta - 1.0;
	Bprt = false;
	Bprt_sum = false;
	//�z��p�������m��

	int i ;
	ErrorMethod			= new int[n];
	Error				= new double[n];
	answer				= new double[n];
	residual			= new double[n];
	set_value			= new double[n];
	MaxDisplacement		= new double[n];
	_MaxDisplacement	= new double[n];
	MaxValue			= new double[n];
	MinValue			= new double[n];
	assumed_value		= new double[n];
	_assumed_value		= new double[n];
	step_value			= new double[n];
	jacob				= new double*[n];
	ErrorAbs			= new double*[n+1];
	ErrorRel			= new double*[n+1];
	for( i = 0 ; i < n+1 ; i++ ){
		ErrorAbs[i] = new double[n];
		ErrorRel[i] = new double[n];
	}
	for( i = 0 ; i < n ; i++ ){
		jacob[i] = new double[n];
	}

	for( i = 0 ; i < n ; i++ ){
		ErrorMethod[i] = 0;
		Error[i] = err;
	}

}
void CNewtonRaphsonMethod::setDeltaError( double _delta , double _err ){
	delta = _delta;
	err	  = _err;
	MinimumEPS = _delta - 1.0;
	for( int i = 0 ; i < ArrayNum ; i++ ){
		Error[i] = err;
	}
}

void CNewtonRaphsonMethod::setErrorMethod( int _ErrorMethod ){
	int i ;
	for( i = 0 ; i < ArrayNum ; i++ ){
		ErrorMethod[i] = _ErrorMethod;
	}
}
void CNewtonRaphsonMethod::setAcc( double _acc ){
	acc = _acc;
}
void CNewtonRaphsonMethod::setSolver( int _Solver ){
	Solver = _Solver;
}
void CNewtonRaphsonMethod::setEscLoop( int _EscLoop ){
	EscLoop = _EscLoop;
}
void CNewtonRaphsonMethod::reset(){
	FlagInitial = false;
	_num = 0;
}

void CNewtonRaphsonMethod::set_assumed_value(int i,double value,double _max)
{

	if( FlagInitial == true ){
		reset();
	}

	if( i+1 > _num ){
		_num = i+1;
	}
	_assumed_value[i] = value;
	_MaxDisplacement[i] = _max;

}
void CNewtonRaphsonMethod::setValue(int i,double value,double _max , double _MinValue , double _MaxValue )
{

	if( FlagInitial == true ){
		reset();
	}

	if( i+1 > _num ){
		_num = i+1;
	}
	_assumed_value[i] = value;
	_MaxDisplacement[i] = _max;

	MinValue[i] = _MinValue;
	MaxValue[i] = _MaxValue;


}
double CNewtonRaphsonMethod::getValue(int i){
	return set_value[i];
}

void CNewtonRaphsonMethod::initial(int n)
{
	int i;
	if( n == 0 ){
		num = _num;
	}else{
		num = n;
	}

//	matA = boost::numeric::ublas::matrix<double>( num , num );
//	vecX = boost::numeric::ublas::vector<double>( num );
//	vecB = boost::numeric::ublas::vector<double>( num );

	if( FlagInitial == false ){
		FlagInitial = true;
		for( i = 0 ; i < num ; i++ ){
			assumed_value[i] = _assumed_value[i];
			MaxDisplacement[i]	 = _MaxDisplacement[i];
		}
	}

	ErrorRel[num][0] = 100000000.0;
	ErrorAbs[num][0] = 100000000.0;
	err_sum = num * 100;
	err_ave = err_sum / num;
	ErrorLoopOver = 1;
}
void CNewtonRaphsonMethod::main_loop_init()
{
	//�l�̏�����
	if( FlagInitial == false ){
		std::cout<<"error : initialize"<<std::endl;
		getchar();
		exit(1);
	}
	check = false;
	check2 = false;
	count = 1;

	int i,j;
	for( i = 0 ; i < ArrayNum+1 ; i++ ){
		for( j = 0 ; j < ArrayNum ; j++ ){
			ErrorAbs[i][j] = 100.0;
			ErrorRel[i][j] = 100.0;
		}
	}

}
bool CNewtonRaphsonMethod::main_loop_check()//���ꂪ�Ō�̏���
{


	if( count > EscLoop ){
		std::cout<< name << "loop max skip "<< "error_max = " << err_max << std::endl;
		ErrorLoopOver = 0;
		return (false);
	}

	return (!check);

/*
	//���v�G���[�̌v�Z
	err_max = 0.0;
	err_sum = 0.0;
	for(int i=0;i<num;i++){
		if( ErrorMethod[i] == 0 ){
			err_sum  += fabs(ErrorRel[num][i]);
			if( fabs(ErrorRel[num][i]) > err_max ){
				err_max = fabs(ErrorRel[num][i]);
			}
		}else if( ErrorMethod[i] == 1 ){
			err_sum  += fabs(ErrorAbs[num][i]);
			if( fabs(ErrorAbs[num][i]) > err_max ){
				err_max = fabs(ErrorAbs[num][i]);
			}
		}
	}
	err_ave = err_sum / num;


	if( Bprt == true ){
		prt2();
	}
	if( Bprt_sum == true ){
		prt_sum2();
	}


	//�������f
	check = true;
	if( count > EscLoop ){
		std::cout<< name << "loop max skip "<< "error_max = " << err_max << std::endl;
		ErrorLoopOver = 0;
		return (!check);
	}

	for( int i=0 ; i <= num-1 ; i++ )
	{
		if( ErrorMethod[i] == 0 ){
			if( fabs(ErrorRel[num][i] ) >= Error[i] ){
				check = false;
				break;
			}
		}else if( ErrorMethod[i] == 1 ){
			if( fabs(ErrorAbs[num][i] ) >= Error[i] ){
				check = false;
				break;
			}
		}

	}
	return (!check);
	*/

}
void CNewtonRaphsonMethod::main_loop_reinit()
{
	int i;
	//���v�G���[�̌v�Z
	err_max = 0.0;
	err_sum = 0.0;
	for(i=0;i<num;i++){
		if( ErrorMethod[i] == 0 ){
			err_sum  += fabs(ErrorRel[num][i]);
			if( fabs(ErrorRel[num][i]) > err_max ){
				err_max = fabs(ErrorRel[num][i]);
			}
		}else if( ErrorMethod[i] == 1 ){
			err_sum  += fabs(ErrorAbs[num][i]);
			if( fabs(ErrorAbs[num][i]) > err_max ){
				err_max = fabs(ErrorAbs[num][i]);
			}
		}
	}
	err_ave = err_sum / num;

	if( Bprt == true ){
		prt2();
	}
	if( Bprt_sum == true ){
		prt_sum2();
	}

	//�������f
	check = true;
	for(i=0 ; i <= num-1 ; i++ )
	{
		if( ErrorMethod[i] == 0 ){
			if( fabs(ErrorRel[num][i] ) >= Error[i] ){
				check = false;
				break;
			}
		}else if( ErrorMethod[i] == 1 ){
			if( fabs(ErrorAbs[num][i] ) >= Error[i] ){
				check = false;
				break;
			}
		}
	}
	if( check == true ){
		return;//�������Ă����ꍇ�͏C���l���v�Z������main_loop_check()��
	}



	int gyou,retu;
	double tmp,sum;
	//�Δ����̎Z�o
	for( gyou = 0 ; gyou <= num-1 ; gyou++ )
	{
		residual[gyou] = ErrorAbs[num][gyou];

		for( retu = 0 ; retu <= num-1 ; retu++ )
		{
			jacob[gyou][retu] = ( ErrorAbs[num][gyou] - ErrorAbs[retu][gyou] )
				/ (assumed_value[retu] - step_value[retu] );
		}
	}

	//�A���ꎟ������
	if( Solver == 0 ){
		ErrorSolver = gaussian_elimination();
	}else if( Solver == 1 ){
		ErrorSolver = gaussian_elimination_2();
	}else if( Solver == 2 ){
//		ErrorSolver = lu_factorize_ublas();
	}else{
		ErrorSolver = gaussian_elimination();
	}


	if( ErrorSolver == 0 )
	{
		std::cout<<"error : solver"<<std::endl;
		getchar();
		exit(1);
	}
	//�l�̕ύX
	for( i=0; i<=num-1; i++ )
	{
		if(MaxDisplacement[i]==0.0){
			assumed_value[i] = set_value[i] - acc*answer[i];
		}
		else{
			double temp2;
			temp2 = acc * answer[i];
			tmp = fabs(temp2);
			if( tmp > MaxDisplacement[i] ){
				if( answer[i]>0 ){
					assumed_value[i] = set_value[i] - MaxDisplacement[i];
				}else{
					assumed_value[i] = set_value[i] + MaxDisplacement[i];
				}
			}else{
				assumed_value[i] = set_value[i] -  temp2;
			}
		}

		if( MinValue[i] != 0.0 && assumed_value[i] < MinValue[i] ){
			assumed_value[i] = MinValue[i];
		}
		if( MaxValue[i] != 0.0 && assumed_value[i] > MaxValue[i] ){
			assumed_value[i] = MaxValue[i];
		}



	}
	//�J�E���^��������
	if(count==count_max){
		sum = 0;
		for(i=0;i<num;i++){
			sum += fabs(ErrorAbs[num][i]);
		}
		if(sum<sqrt(err)){
			std::cout<<"loop over"<<std::endl;
			check=true;
		}
		else std::cout<<"error : loop over"<<std::endl;
		getchar();
		exit(1);
	}



	count++;
	
}
void CNewtonRaphsonMethod::spcacc(double maxacc){

	acc = maxacc * pow( 1.5 , count - 10.0 );
	if( acc > maxacc ){
		acc = maxacc;
	}
}
void CNewtonRaphsonMethod::sub_loop_init()
{
	//�l�̏�����
	i_sys   = 0;
	i_error = 0;
//	double MinimumEPS = 1.0;
	//�l�̑��
	for(int i=0; i<=num-1; i++ )
	{
		if(fabs(assumed_value[i])>MinimumEPS)	step_value[i] = assumed_value[i]*delta;
		else									step_value[i] = assumed_value[i]+(assumed_value[i]>0?MinimumEPS:-MinimumEPS)*(delta-1.0);
//		step_value[i] = assumed_value[i]*delta;

		set_value[i] = (i==i_sys)?step_value[i]:assumed_value[i];
	}
}
bool CNewtonRaphsonMethod::sub_loop_check()
{
	return ( i_sys <= num );
}
void CNewtonRaphsonMethod::sub_loop_reinit()
{
	i_sys++;
	//�l�̏�����
	i_error = 0;
	//�l�̑��
	if( i_sys != num + 1 )
	{
		for(int i=0; i<=num-1; i++ )
		{
			set_value[i] = (i==i_sys)?step_value[i]:assumed_value[i];
		}
	}

}

//�G���[�l�̑��
void CNewtonRaphsonMethod::set_error(double atai,double bairitu)
{
	ErrorAbs[i_sys][i_error] = atai*bairitu;
	ErrorRel[i_sys][i_error] = ErrorAbs[i_sys][i_error];
	i_error++;
}
void CNewtonRaphsonMethod::set_error2(int bangou,double atai1,double atai2,double bairitu)
{
	ErrorAbs[i_sys][bangou] = ( atai1 - atai2 ) * bairitu;
/*
	if( fabs(atai1) < 1e-3 || fabs(atai2) < 1e-3 ){
		if( atai1 > 0.0 ){
			atai1 += 1e-3;
			atai2 += 1e-3;
		}else{
			atai1 -= 1e-3;
			atai2 -= 1e-3;
		}
	}*/

	ErrorRel[i_sys][bangou] = ( 1.0 - ( atai1 / atai2 ) ) * bairitu;
}

void CNewtonRaphsonMethod::setError(int i,double Value1,double Value2, int _ErrorMethod , double error )
{
	if( error > 0.0 ){
		Error[i] = error;
	}else{
		Error[i] = err;
	}

	ErrorMethod[i] = _ErrorMethod;

	ErrorAbs[i_sys][i] = ( Value1 - Value2 );

	//���̏������������̂ŃR�����g���@by���
	if( fabs(Value1) < 1e-3 || fabs(Value2) < 1e-3 ){
		if( Value2 > 0.0 ){
			Value1 += 1e-3;
			Value2 += 1e-3;
		}else{
			Value1 -= 1e-3;
			Value2 -= 1e-3;
		}
	}
	
	ErrorRel[i_sys][i] = ( 1.0 - ( Value1 / Value2 ) );


//	if( ErrorRel[i_sys][i] < -100000.0 || 100000 < ErrorRel[i_sys][i] ){
//		std::cout << "asdf";
//	}


}


//�G���[�l�̕\��
void CNewtonRaphsonMethod::prt(double value)
{
	int i;
	if( i_sys == num )
	{
		std::cout<<"\tnum of loop = "<<count<<std::endl;
		for( i=0; i<=num-1; i++ )
		{
			if( fabs(ErrorAbs[i_sys][i]) > value ){
				std::cout << std::setprecision(6);
				std::cout << "ErrorRel[" << i+1 << "] = " << fabs(ErrorRel[i_sys][i]) << "\t";
				std::cout << "ErrorAbs[" << i+1 << "] = " << fabs(ErrorAbs[i_sys][i]) << "\t";
				std::cout << std::endl;
			}
		}
	}
}
void CNewtonRaphsonMethod::prt_sum()
{
	if(i_sys==num)
	{
		std::cout.precision(6);
		std::cout << "\tnum of loop = "<< count; 
		std::cout << "\tacc = "<< acc;
		std::cout << "\terr_sum = "<< err_sum;
		std::cout << "\terr_ave = " << err_ave;
		std::cout << "\terr_max = " << err_max;
		std::cout << std::endl;
	}
}
void CNewtonRaphsonMethod::prt2(double value)
{
	int i;
	std::cout<<"\tnum of loop = "<<count<<std::endl;
	for( i=0; i<=num-1; i++ )
	{
		if( fabs(ErrorAbs[i_sys][i]) > value ){
			std::cout << std::setprecision(6);
			std::cout << "ErrorRel[" << i+1 << "] = " << fabs(ErrorRel[i_sys][i]) << "\t";
			std::cout << "ErrorAbs[" << i+1 << "] = " << fabs(ErrorAbs[i_sys][i]) << "\t";
			std::cout << std::endl;
		}
	}

}
void CNewtonRaphsonMethod::prt_sum2()
{

	std::cout.precision(6);
	std::cout << "\tnum of loop = "<< count; 
	std::cout << "\tacc = "<< acc;
	std::cout << "\terr_sum = "<< err_sum;
	std::cout << "\terr_ave = " << err_ave;
	std::cout << "\terr_max = " << err_max;
	std::cout << std::endl;

}

void CNewtonRaphsonMethod::matrixFileMake(){
	
	matrixFile.open("matrix.csv" , std::ios::out );

}
void CNewtonRaphsonMethod::matrixOut(){
	


	matrixFile << "jacob and answer" << std::endl;
	for( int i = 0 ; i < num ; i++ ){
		for( int j = 0 ; j < num ; j++ ){
			matrixFile << jacob[i][j] << ",";
		}
		matrixFile << answer[i] << ",";
		matrixFile << std::endl;
	}
	matrixFile << std::endl;

	matrixFile << "REL error" << std::endl;
	for( int i = 0 ; i < num ; i++ ){
		for( int j = 0 ; j < num ; j++ ){
			matrixFile << ErrorRel[i][j] << ",";
		}
		matrixFile << std::endl;
	}
	matrixFile << std::endl;

	matrixFile << "ABS error" << std::endl;
	for( int i = 0 ; i < num ; i++ ){
		for( int j = 0 ; j < num ; j++ ){
			matrixFile << ErrorAbs[i][j] << ",";
		}
		matrixFile << std::endl;
	}
	matrixFile << std::endl;


}

//�K�E�X�̏����@�@�����s�{�b�g�I��ver
int CNewtonRaphsonMethod::gaussian_elimination_2(){
	int i,j,k;
	double max;
	double keisu;
	int si,sj,maxi;//�s�{�b�g�I�����Ɏg�p
	for( k = 0 ; k < num ; k++ ){
		//�s�{�b�g�I��
		max = 0.0;
		for( si = k ; si < num ; si++ ){
			if( max <= fabs( jacob[si][k] ) ){
				maxi = si;
				max  = fabs( jacob[si][k] );
			}
			if( max > 1.0 ){
				break;
			}
		}
		if( max <= pibot ){
			return(1);
		}
		//����ւ�
		if( maxi != k ){
			for( sj = 0 ; sj < num ; sj++ ){
				irekae( jacob[k][sj] , jacob[maxi][sj] );
			}
			irekae( residual[k] , residual[maxi] );
		}
		for( i = k+1 ; i < num ; i++ ){
			if( fabs(jacob[i][k]) > 1e-12 ){
				keisu = jacob[i][k] / jacob[k][k];
				for( j = k+1 ; j < num ; j++ ){
					jacob[i][j] -= jacob[k][j] * keisu; 
				}
				residual[i] -= residual[k] * keisu;
				jacob[i][k] = 0.0;
			}
		}
	}
	//��ޑ��
	for( k = num-1 ; k >= 1 ; k-- ){
		for( i = k-1 ; i >= 0 ; i-- ){
			residual[i] -= residual[k] * jacob[i][k] / jacob[k][k];
			jacob[i][k] = 0.0;
		}
	}
	//�����v�Z
	for( k = 0 ; k < num ; k++ ){
		answer[k] = residual[k] / jacob[k][k];
	}
	return 1;
}
/*
int CNewtonRaphsonMethod::lu_factorize_ublas(){
	int gyou;
	int retu;
	for( gyou = 0 ; gyou < num ; gyou++ ){
		for( retu = 0 ; retu < num ; retu++ ){
			matA(gyou,retu) = jacob[gyou][retu];
		}
	}
	for( gyou = 0 ; gyou < num ; gyou++ ){
		vecB(gyou) = residual[gyou];
	}
	boost::numeric::ublas::matrix<double> lhs(matA);
	boost::numeric::ublas::vector<double> rhs(vecB);
	boost::numeric::ublas::permutation_matrix<> pm(matA.size1());
	boost::numeric::ublas::lu_factorize(lhs, pm);
	boost::numeric::ublas::lu_substitute(lhs, pm, rhs);
	vecX.assign_temporary(rhs);
	for( gyou = 0 ; gyou < num ; gyou++ ){
		answer[gyou] = vecX(gyou);
	}

	return 1;
}
*/




/*
 *  Numerical Calculation : �K�E�X�̏����@�i���S�s�{�b�e�B���O�j ver5.0 G98�x���Y��
 *
 *  keisuu * x = y �ikeisuu,y:���m,x:���m�j�������D
 *
 *  ��������������������������������
 *  ����              ���^    �����l
 *  ��������������������������������
 *  �W���s��ikeisuu�j��type *��*1
 *  �萔�s��(y)       ��type *��*2
 *  ��(x)             ��type *��*2
 *  �W���s��̎���(n) ��int   ��
 *  �ŏ��s�{�b�g�l    ��type  ��
 *  ��������������������������������
 *  *1:�W���s��̎����Ɣz��̎����𓙂������Ă��������D
 *     �܂�2x2�s��Ȃ�a[2][2]�Ȃǂ������a[5][5]�Ȃǂ͓���Ȃ��ŉ�����
 */

void CNewtonRaphsonMethod::irekae(double &x,double &y ){double z=x;x=y;y=z;}
void CNewtonRaphsonMethod::irekae(int    &x,int    &y ){int    z=x;x=y;y=z;}

#define __GEM(i,j)  ( (j!=num) ? (jacob[i][j]) : (residual[i]) )
int CNewtonRaphsonMethod::gaussian_elimination()
{
	int i,j,k;
	double max;
	int si,sj,maxi,maxj;//�s�{�b�g�I�����Ɏg�p
	int *x_narabi = new int [num];//���̕��я�
	for( i=0; i<=num-1; i++ )
	{
		x_narabi[i]=i;
	}
	double *x_tmp = new double [num];//�����ꎞ�I�ɃX�g�b�N

	for( k=0; k<=num-1; k++ )
	{
		//���S�s�{�b�g�I��
		max = 0.0;
		for( si=k; si<=num-1; si++ )
		{
			for( sj=k; sj<=num-1; sj++ )
			{
				if( max <= fabs( __GEM(si,sj) ) )
				{
					max  = fabs( __GEM(si,sj) );
					maxi = si;
					maxj = sj;
				}
			}
		}
		if( max <= pibot )
		{
			return (1);
		}
		if( maxi != k )//�s���ꊷ��
		{
			for( sj=k; sj<=num; sj++ )
			{
				irekae( __GEM(k,sj), __GEM(maxi,sj) );
			}
		}
		if( maxj != k )//����ꊷ��
		{
			for( si=0; si<=num-1; si++ )
			{
				irekae( __GEM(si,k), __GEM(si,maxj) );
			}
			irekae( x_narabi[k], x_narabi[maxj] );
		}

	
		//�ΏۂƂȂ�s�𑀍�
		for( j=k+1; j<=num; j++ )
		{
			__GEM(k,j) /= __GEM(k,k);
		}
		//�ΏۂƂȂ�s��艺���ɂ���s�𑀍�
		for( i=k+1; i<=num-1; i++ )
		{
			if( __GEM(i,k) != 0.0 ){
				for( j = k+1 ; j < num+1 ; j++ )
				{
					__GEM(i,j) -= __GEM(i,k) *__GEM(k,j);
				}
			}
		}
	}
	//���̎Z�o
	for( k=num-1; k>=0; k-- )
	{
		x_tmp[k] = residual[k];
		for( j=k+1; j<=num-1; j++ )
		{
			x_tmp[k] -= x_tmp[j] * __GEM(k,j);
		}
	}
	for( i=0; i<=num-1; i++ )
	{
		answer[ x_narabi[i] ] = x_tmp[i];
	}
	//���������
	delete [] x_narabi;
	delete [] x_tmp;
	return 1;
}

#undef __GEM
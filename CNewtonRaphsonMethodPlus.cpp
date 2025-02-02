/**************************************************************************************************
*
*	���ϐ��j���[�g���@ �F CNewtonRaphsonMethodPlus
*
*	1G04 ���c�S
*
*	CNewtonRaphsonMethod������v�Z�⃄�R�r�s��쐬�X�L�b�v�@�\�ɑΉ����ɂ������߁C��蒼�����D
*   
*   [��ȕύX�_]
*   �Eopenmp�ɂ�����v�Z�ɑΉ�
*   �E���R�r�s��̍쐬���X�L�b�v����@�\�ɑΉ�
*   �E�����̍s���vector�^�Ő錾
*   
*	[��{�I�Ȏg�p�@]
*

#include <iostream>
#include "CNewtonRaphsonMethodPlus.h"

void main(){

	//�N���X��錾
	CNewtonRaphsonMethodPlus nrm;

	double a;
	double b;
	double c;

//------------------------�ȉ��m�[�}���Ȏg�p���@---------------------------//

	a = 1.0;
	b = 1.0;
	c = 1.0;

	//������
	nrm.Initialize();

	//�ϐ��̃Z�b�g
	nrm.SetValiable( 0 , a );
	nrm.SetValiable( 1 , b );
	nrm.SetValiable( 2 , c );

	//�e��Z�b�g�A�b�v

	//�����W���̐ݒ�
	nrm.SetAcc( 0.5 );

	//1���[�v���Ƃ̃G���[�l�̃f�B�X�v���C�ݒ�
	//nrm.SetPrint( �S�̃G���[�\�� , �ϐ����Ƃ̃G���[�\�� );
	nrm.SetPrint( true , true );

	//�����ړ��ʂ̐ݒ�
	nrm.SetDelta( 1.0001 );

	//��������̐ݒ�
	nrm.SetError( 0.0001 );

	for(nrm.MainLoopInit();nrm.MainLoopCheck();nrm.MainLoopReinit()){
		for(nrm.SubLoopInit();nrm.SubLoopCheck();nrm.SubLoopReinit()){
			
			//�C���l�̎󂯎��
			nrm.GetValiable( 0 , a );
			nrm.GetValiable( 1 , b );
			nrm.GetValiable( 2 , c );

			//�G���[�l�̑��
			nrm.SetResult( 0 , a * b , 10.0 );
			nrm.SetResult( 1 , b * c , 20.0 );
			nrm.SetResult( 2 , c * a , 30.0 );
		
		}
	}

	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << "c = " << c << std::endl;

//------------------------�ȏ�m�[�}���Ȏg�p���@---------------------------//


//------------------------�ȉ��A�h�o���X�h�Ȏg�p���@---------------------------//

	a = 1.0;
	b = 1.0;
	c = 1.0;

	//������
	nrm.Initialize();

	//�ϐ��̃Z�b�g
	nrm.SetValiable( 0 , a );
	nrm.SetValiable( 1 , b );
	nrm.SetValiable( 2 , c );

	//�e��Z�b�g�A�b�v

	//�����W���̐ݒ�
	nrm.SetAcc( 0.5 );

	//1���[�v���Ƃ̃G���[�l�̃f�B�X�v���C�ݒ�
	nrm.SetPrint( false , false );

	//�����ړ��ʂ̐ݒ�
	nrm.SetDelta( 1.0001 );

	//��������̐ݒ�
	nrm.SetError( 0.000001 );

	//���R�u�s��̍쐬���X�L�b�v����X�e�b�v����ݒ�D
	//���̏ꍇ�ł́C�j���[�g���@5���[�v�͓������R�r�s����g�p����
	nrm.SetJacobStep( 5 );

	//�j���[�g���@�̏I�����ɍŏI�G���[���f�B�X�v���C�ɕ\������
	nrm.SetEndReport( true ) ;


	for(nrm.MainLoopInit();nrm.MainLoopCheck();nrm.MainLoopReinit()){
		for(nrm.SubLoopInit();nrm.SubLoopCheck();nrm.SubLoopReinit()){
			
			//�C���l�̎󂯎��
			nrm.GetValiable( 0 , a );
			nrm.GetValiable( 1 , b );
			nrm.GetValiable( 2 , c );

			//�G���[�l�̑��
			nrm.SetResult( 0 , a * b , 10.0 );
			nrm.SetResult( 1 , b * c , 20.0 );
			nrm.SetResult( 2 , c * a , 30.0 );
		
		}
	}

	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << "c = " << c << std::endl;

//------------------------�Ȉȏ�A�h�o���X�h�Ȏg�p���@---------------------------//
	
}

*	�ȉ��̕��͂�CNewtonRaphsonMethod�̃R�s�y�ł����CCNewtonRaphsonMethodPlus�ł��قړ����@�\�����낦�Ă��܂��D
*
*	��1:	����nrm�Ƃ������O���قȂ�ΈقȂ������ϐ��j���[�g���@�ƔF������܂��D���ϐ��j���[�g���@�̒��ɑ��ϐ��j���[�g���@��
*			�\�z�������Ƃ��͂����̖��O�����ς������̂��g�̑��ϐ��j���[�g���@�̒��ɃR�s�y����΂n�j�ł��D
*	��2:	setup����ہC�K�v�ɉ����Ĕ����ړ���(������DELTA)�Ǝ����덷(������ERROR)�ƍŏ��s�{�b�g�l�i�K�E�X�̏����@�ɂ���
*			�ē��ʂȈӖ������l�D�ʏ�͂��̒l���C�ɂ���K�v�͂���܂���D�������Ă���l�������̒l���Βl�Őݒ肵�Ă��������j
*			���w��ł��܂��D�w�肵�����ꍇ�̓p�����[�^�̐��̂��ƂɁC�����ړ��ʁC�����덷�C�ŏ��s�{�b�g�l�̏��ŋL�q���܂��D��
*			��
*			mnm.setup(2,1+1e-3,1e-3);     //���̏ꍇ�����ړ���:1+1e-3�C�����덷:1e-3
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
#include "CNewtonRaphsonMethodPlus.h"


#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)

/*
memo

for( 1 init ; 2 check ; 3 reinit ){
	4
}

12/432/432/432/432/432/43......


1������
2�����o�����聨����l�̌v�Z
4���W���[���v�Z
3�������聨�K�E�X�̏����@���C���l�̌v�Z




������̃��X�g
�E����l�̗����̕ۑ�
�E������
�E�֐��̒���
�Ekeisu��z�񉻂���
�E���ƕϐ��̐�������Ȃ��������̃G���[
�E�T�u���[�v�񐔂��w�肷��ƕϐ����o���Ă����֐�

*/





#include "CNewtonRaphsonMethodPlus.h"



int CNewtonRaphsonMethodPlus::MatrixOut(){


	int i;
	int j;


	MatrixFile.open(Name + "matrix.csv", std::ios::out);

	if (!MatrixFile){

		std::cout << "file open error" << std::endl;
		std::cout << Name << "�̃}�g���b�N�X�t�@�C�����J���܂���" << std::endl;
		system("pause");
		exit(0);

	}


	//���R�r�s��̃A�E�g�v�b�g

	MatrixFile << "JacobMatrix" << std::endl;

	MatrixFile << "#" << ",";
	for (j = 0; j < NumVariable; j++){
		MatrixFile << "#eq" << j << ",";
	}
	MatrixFile << "Results"<<std::endl;

	for ( i = 0; i < NumVariable; i++){
		MatrixFile << "#Variable_" << i << "_" << VariableName[i] << ",";
		for (j = 0; j < NumVariable + 1; j++){
			MatrixFile << JacobMatrix[i][j] << ",";
		}
		MatrixFile << std::endl;
	}
	MatrixFile << std::endl;





	//��O�p�s��̃A�E�g�v�b�g

	MatrixFile << "UpperTriangularMatrix" << std::endl;

	MatrixFile << "#" << ",";
	for (j = 0; j < NumVariable; j++){
		MatrixFile << "#eq" << j << ",";
	}
	MatrixFile << "Results" << std::endl;

	for (i = 0; i < NumVariable; i++){
		MatrixFile << "#Variable_" << i << "_" <<VariableName[i] << ",";
		for ( j = 0; j < NumVariable + 1; j++){
			MatrixFile << UpperTriangularMatrix[i][j] << ",";
		}
		MatrixFile << std::endl;
	}
	MatrixFile << std::endl;



	//�Ίp�s��̃A�E�g�v�b�g

	MatrixFile << "DiagonalMatrix" << std::endl;

	MatrixFile << "#" << ",";
	for (j = 0; j < NumVariable; j++){
		MatrixFile << "#eq" << j << ",";
	}
	MatrixFile << "Results" << std::endl;

	for (i = 0; i < NumVariable; i++){
		MatrixFile << "#Variable_" << i << "_" <<VariableName[i] << ",";
		for ( j = 0; j < NumVariable + 1; j++){
			MatrixFile << DiagonalMatrix[i][j] << ",";
		}
		MatrixFile << std::endl;
	}
	MatrixFile << std::endl;






	MatrixFile.close();



	return 1;
}



int CNewtonRaphsonMethodPlus::SetName( std::string _Name ){
	Name = _Name;
	return 1;
}

/*
 * ���������s��
 *
 *
 */
int CNewtonRaphsonMethodPlus::Initialize(){
	NumVariable = 0;
	FlagJacobSkip = false;
	return 1;
}

int CNewtonRaphsonMethodPlus::SetAccAuto(bool _FlagAccAuto, double _AccAutoMin , double _AccAutoMax , double _AccAutoStep ){
	FlagAccAuto = _FlagAccAuto;
	AccAutoMin = _AccAutoMin;
	AccAutoMax = _AccAutoMax;
	AccAutoStep = _AccAutoStep;


	return 1;

}

int CNewtonRaphsonMethodPlus::SetAcc( double _Acc ){
	int i;
	for( i = 0 ; i < NumVariable ; i++ ){
		SetAccDirect( i , _Acc );
	}
	return 1;
}

int CNewtonRaphsonMethodPlus::SetAccDirect( int _Num , double _Acc ){
	Acc[ _Num ] = _Acc;

	return 1;
}

int CNewtonRaphsonMethodPlus::SetEndReport( bool _Flag ){
	FlagEndReport = _Flag;
	return 1;
}

int CNewtonRaphsonMethodPlus::ValueContainerResize( int _NumVariable ){
	Variable.resize( _NumVariable+1 );//�ϐ����ꎩ��
	VariableName.resize( _NumVariable + 1 );
	VariableMemory.resize( _NumVariable+1 );//�ϐ��̗���


	Acc.resize( _NumVariable+1 );//�ϐ��̉����W��
	MaxValue.resize( _NumVariable+1 );//�ϐ��̍ő�l
	MinValue.resize( _NumVariable+1 );//�ϐ��̍ŏ��l
	MaxDisplacement.resize( _NumVariable+1 );//�ϐ��̍ő�ړ���

	MaxValue[_NumVariable] = 0.0;
	MinValue[_NumVariable] = 0.0;
	MaxDisplacement[_NumVariable] = 0.0;

	StepControl.resize( _NumVariable+1 );
	MinMaxControl.resize( _NumVariable+1 );
	StepControl[_NumVariable] = false;
	MinMaxControl[_NumVariable] = false;

	Step.resize( _NumVariable+1 );//�ϐ��̍ő�ړ���


	for( int i = NumVariable ; i < _NumVariable+1 ; i++ ){
		Variable[ i ].resize(1);
		VariableMemory[ i ].resize(1);
		Acc[ i ] = 1.0;//Acc��1�ŏ���������D
	}
	return 1;
}

int CNewtonRaphsonMethodPlus::SetValiable( int _NumVariable , double _Variable ){

	if( NumVariable < _NumVariable+1 ){
		ValueContainerResize( _NumVariable );
		NumVariable = _NumVariable+1;
	}
	
	Variable[ _NumVariable ][ 0 ] = _Variable;
	VariableMemory[ _NumVariable ][ 0 ] = _Variable;
//	StepControl[_NumVariable] = false;
//	MinMaxControl[_NumVariable] = false;	
	return 1;
}


int CNewtonRaphsonMethodPlus::SetValiableName(int _NumVariable, std::string _VariableName){

	VariableName[_NumVariable]= _VariableName;

	return 1;
}


int CNewtonRaphsonMethodPlus::SetValiable( int _NumVariable , double _Variable , double _MaxStep ){
	SetValiable(  _NumVariable ,  _Variable );
	SetMaxStep(  _NumVariable ,  _MaxStep );
	return 1;
}

int CNewtonRaphsonMethodPlus::SetValiable( int _NumVariable , double _Variable , double _MaxStep , double _MinValue , double _MaxValue ){
	SetValiable(  _NumVariable ,  _Variable );
	SetMaxStep(  _NumVariable ,  _MaxStep );
	SetMinMaxValue(  _NumVariable ,  _MinValue , _MaxValue );
	return 1;
}

int CNewtonRaphsonMethodPlus::SetMaxStep( int _NumVariable , double _MaxStep ){
	MaxDisplacement[ _NumVariable ] = _MaxStep;
//	StepControl[_NumVariable] = true;

	return 1;
}
int CNewtonRaphsonMethodPlus::SetMinMaxValue( int _NumVariable , double _MinValue , double _MaxValue ){
	MinValue[ _NumVariable ] = _MinValue;
	MaxValue[ _NumVariable ] = _MaxValue;
//	MinMaxControl[_NumVariable] = true;
	return 1;
}


int CNewtonRaphsonMethodPlus::GetValiable( int _NumVariable , double &_Variable ){
	_Variable = Variable[ _NumVariable ][ SubLoopCounter ];
	return 1;
}
int CNewtonRaphsonMethodPlus::GetValiableDirect( int _NumVariable , int _SubLoopCounter, double &_Variable ){
	_Variable = Variable[ _NumVariable ][ _SubLoopCounter ];
	return 1;
}

int CNewtonRaphsonMethodPlus::SetResult( int _NumVariable , double _VariableA , double _VariableB ){
	ResultMatrix[ _NumVariable ][ SubLoopCounter ] = _VariableA - _VariableB;
	return 1;
}
int CNewtonRaphsonMethodPlus::SetResultDirect( int _NumVariable , int _SubLoopCounter , double _VariableA , double _VariableB ){
	ResultMatrix[ _NumVariable ][ _SubLoopCounter ] = _VariableA - _VariableB;
	return 1;
}


int CNewtonRaphsonMethodPlus::SetPrint( bool _Print , bool _PrintAll ){
	FlagPrint = _Print;
	FlagPrintAll = _PrintAll;
	return 1;
}

int CNewtonRaphsonMethodPlus::SetMatrixOut(bool _MatrixOut){
	FlagMatrixOut = _MatrixOut;
	return 1;
}



int CNewtonRaphsonMethodPlus::MainLoopInit(){
	int i;


	time(&TimeStart);

	FlagConbergence = false;
	MainLoopCounter = 0;
	JacobSkipCounter = 0;
	JacobReset = false;



	//�z��̃T�C�Y�����߂�

	for( i = 0 ; i < NumVariable ; i++ ){
		Variable[i].resize( NumVariable+1 );
	}


	ResultMatrix.resize( NumVariable );
	JacobMatrix.resize( NumVariable );
	UpperTriangularMatrix.resize( NumVariable );
	DiagonalMatrix.resize( NumVariable );
	UpperC.resize( NumVariable );

	RelError.resize( NumVariable );
	AbsError.resize( NumVariable );
	DeltaValue.resize( NumVariable );

	for( i = 0 ; i < NumVariable ; i++ ){
		ResultMatrix[i].resize( NumVariable+1 );
		JacobMatrix[i].resize( NumVariable+1 );
		UpperTriangularMatrix[i].resize( NumVariable+1 );
		DiagonalMatrix[i].resize( NumVariable+1 );
	}


	//�ŏ��̎�������Ŕ�΂���Ȃ��悤�ɃG���[�l��100�����Ă���
	ResultMatrix[0][NumVariable] = 100.0;

	for( i = 0 ; i < NumVariable ; i++ ){
		Variable[i][ NumVariable ] = Variable[i][0];
		DiagonalMatrix[i][ NumVariable ] = 0.0;
	}

	return 1;
}

int CNewtonRaphsonMethodPlus::MainLoopCheck(){

	//Error�̒l������
	MaxError = 0.0;
	AveError = 0.0;
	for( int i = 0 ; i < NumVariable ; i++ ){
		if( MaxError < fabs( ResultMatrix[i][NumVariable] ) ){
			MaxError = fabs( ResultMatrix[i][NumVariable] );
		}
		AveError += fabs( ResultMatrix[i][NumVariable] );
	}
	AveError /= NumVariable;


	//error�ɂ��Acc�̎�������
	ErrorMax0 = ErrorMax1;
	ErrorMax1 = MaxError;

	ErrorDif0 = ErrorDif0;
	ErrorDif1 = ErrorMax1 - ErrorMax0;

	if (FlagAccAuto){
		if (ErrorDif1 < 0.0){
		}else{
			double _Acc = Acc[0];
			if (_Acc / AccAutoStep < AccAutoMin){
				SetAcc( AccAutoMin);
			}
			else{
				SetAcc(_Acc / AccAutoStep);
			}
		}
	}


	//����������s��
	FlagConbergence = true;
	if( Error < MaxError ){
		FlagConbergence = false;
	}else{
		FlagConbergence = true;
	}


	//�G���[�f�B�X�v���C
	if( FlagPrint == true ){
		Print();
	}
	if( FlagPrintAll == true ){
		PrintAll();
	}


	//�������Ă����ꍇ��0��Ԃ��D
	if( FlagConbergence == true ){
		if( FlagEndReport == true ){
			EndReport();
		}
		return 0;
	}


	//MaxLoop�𒴂��Ă����甲���o���D
	if( MaxLoop < MainLoopCounter ){
		if( FlagEndReport == true ){
			EndReport();
		}
		//std::cout << "MaxLoop Over" << std::endl;
		return 0;
	}


	//�C���l������Ƃ���
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for( int i = 0 ; i < NumVariable ; i++ ){

		Step[i] = Acc[i] * DiagonalMatrix[i][ NumVariable ];

		if( MaxDisplacement[i] != 0.0 ){
			if( MaxDisplacement[i] < Step[i]  ){
				Step[i] = MaxDisplacement[i];
			}else if( Step[i] < -MaxDisplacement[i]){
				Step[i] = -MaxDisplacement[i];
			}
		}

		//���ƂŎg����������Ȃ�����c���Ă���
		//if ((Variable[i][NumVariable] - Step[i])*(Variable[i][NumVariable]) < 0.0){
		//	JacobReset = true;
		//}
		//else if ( 0.9 < fabs(Step[i]/Variable[i][NumVariable]) ){
		//	JacobReset = true;
		//}
		//
		Variable[i][NumVariable] = Variable[i][NumVariable] - Step[i];

		if( MinValue[i] != 0.0 ){
			if( Variable[i][NumVariable] < MinValue[i] ){
				Variable[i][NumVariable] = MinValue[i];
			}
		}
		if( MaxValue[i] != 0.0 ){
			if( Variable[i][NumVariable] > MaxValue[i] ){
				Variable[i][NumVariable] = MaxValue[i];
			}
		}

		for( int k = 0 ; k < NumVariable ; k++ ){
			Variable[i][k] = Variable[i][NumVariable];
		}

		int j;
		j = i;
		//�����ړ��ʂ̌v�Z�C����l��(Delta - 1.0)�ȉ��̎��͏�Ɂidelta-1(Delta - 1.0)�����������悤�ɂ���
		//if (fabs(Variable[i][NumVariable])  < (Delta - 1.0)){
		//		if (0 < Variable[i][NumVariable]){
		//		DeltaValue[i] = (Delta - 1.0)*(Delta - 1.0);
		//		Variable[i][j] = Variable[i][NumVariable] + (Delta - 1.0)*(Delta - 1.0);
		//	}
		//	else{
		//		DeltaValue[i] = -(Delta - 1.0)*(Delta - 1.0);
		//		Variable[i][j] = Variable[i][NumVariable] - (Delta - 1.0)*(Delta - 1.0);
		//	}
		//}
		//�����ړ��ʂ̌v�Z�C����l��1�ȉ��̎��͏�Ɂidelta-1�j�����������悤�ɂ���
		if (fabs(Variable[i][NumVariable])  < 1.0){
			if (0 < Variable[i][NumVariable]){
				DeltaValue[i] = (Delta - 1.0);
				Variable[i][j] = Variable[i][NumVariable] + (Delta - 1.0);
			}
			else{
				DeltaValue[i] = -(Delta - 1.0);
				Variable[i][j] = Variable[i][NumVariable] - (Delta - 1.0);
			}
		}
		else{
			DeltaValue[i] = Variable[i][NumVariable] * (Delta - 1.0);
			Variable[i][j] = Variable[i][NumVariable] * Delta;
		}


	}

	return 1;
}

int CNewtonRaphsonMethodPlus::MainLoopReinit(){
	double max;
	int maxi;

	//���R�r�s��\�z���X�L�b�v����ꍇ�͂��̃Z�N�V�������X�L�b�v
	if( FlagJacobSkip == false ){

		//���R�r�s������
		//�����͕��񉻂��Ȃ���������
		for( int i = 0 ; i < NumVariable ; i++ ){
			for( int j = 0 ; j < NumVariable ; j++ ){
				JacobMatrix[i][j] = ( ResultMatrix[i][j] - ResultMatrix[i][NumVariable] ) / DeltaValue[j];
			}
		}
	}

	//���ʂ̍s��͖��񏑂�������
	for( int i = 0 ; i < NumVariable ; i++ ){
		JacobMatrix[i][NumVariable] = ResultMatrix[i][NumVariable];
	}


	//�K�E�X�̏����@�͂��̂����֐�������
	//�������Ă��Ȃ��ꍇ�͂����ŃK�E�X�̏����@�ŏC���l�����߂�D

	//JacobMatrix ==> UpperTriangularMatrix �ւ̃R�s�[
	//�����͕��񉻂��Ȃ���������
	for (int i = 0; i < NumVariable; i++){
		for( int j = 0 ; j < NumVariable+1 ; j++ ){
			UpperTriangularMatrix[i][j] = JacobMatrix[i][j];
		}
	}


	//�O�i����
	for( int k = 0 ; k < NumVariable ; k++ ){

		//�����s�{�b�g�I��
		if( UpperTriangularMatrix[k][k] < Min ){
	
			max = 0.0;
			for( int si = k ; si < NumVariable ; si++ ){
				if (1.0 < max){
					break;
				}
				if( max  < fabs( UpperTriangularMatrix[si][k] ) ){
					maxi = si;
					max  = fabs( UpperTriangularMatrix[si][k] );
				}
			}

			//�s�{�b�g�I���G���[(0����)
			if( max == 0.0 ){
				std::cout << "error : pivot(" << k << ")" << std::endl;
				std::cout << Name << "��NewtonRaphsonMethodPlus�Ńs�{�b�g�I����0�����������܂����D�v���O�������I�����܂�" << std::endl;
				system("pause");
				exit(0);
//				return 0;
			}

			//����ւ�
			if( maxi != k ){
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for( int sj = k ; sj < NumVariable+1 ; sj++ ){
					double temp;
					temp = UpperTriangularMatrix[k][sj];
					UpperTriangularMatrix[k][sj] = UpperTriangularMatrix[maxi][sj];
					UpperTriangularMatrix[maxi][sj] = temp;
				}
			}
		}

		//�O�i��������
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for (int i = k + 1; i < NumVariable; i++){
			UpperC[i] = UpperTriangularMatrix[i][k] / UpperTriangularMatrix[k][k];
			if (fabs(UpperC[i]) > SparseMin ){
				for (int j = k; j < NumVariable + 1; j++){
					UpperTriangularMatrix[i][j] = UpperTriangularMatrix[i][j] - (UpperTriangularMatrix[k][j] * UpperC[i]);
				}
			}
		}		
	}

	//UpperTriangularMatrix ==> DiagonalMatrix �ւ̃R�s�[
	//�����͕��񉻂��Ȃ���������
	for (int i = 0; i < NumVariable; i++){
		for( int j = 0 ; j < NumVariable+1 ; j++ ){
			DiagonalMatrix[i][j] = UpperTriangularMatrix[i][j];
		}
	}

	//��ޑ��
	for( int k = NumVariable-1 ; k >= 1 ; k-- ){
		for( int i = k-1 ; i >= 0 ; i-- ){
			DiagonalMatrix[i][NumVariable] -= DiagonalMatrix[k][NumVariable] * DiagonalMatrix[i][k] / DiagonalMatrix[k][k];
			DiagonalMatrix[i][k] = 0.0;
		}
	}

	//�����v�Z
	//�����͕��񉻂��Ȃ���������
	for (int k = 0; k < NumVariable; k++){
		DiagonalMatrix[k][NumVariable] = DiagonalMatrix[k][NumVariable] / DiagonalMatrix[k][k];
		DiagonalMatrix[k][k] = 1.0;
	}
	

	//�����ŃJ�E���^�[���グ��D
	MainLoopCounter++;
	JacobSkipCounter++;

	//�X�L�b�v�����ɍ��v�����
	if (JacobStep == 1){
		FlagJacobSkip = false;
	}
	else if (JacobReset == true){
		FlagJacobSkip = false;
		JacobReset = false;
		JacobSkipCounter = 0;
	}
	else if (JacobStep <= JacobSkipCounter){
		FlagJacobSkip = false;
		JacobSkipCounter = 0;
	}
	else if (JacobSkipCounter < JacobStep){
		FlagJacobSkip = true;
	}


	if (FlagMatrixOut == true) {
		MatrixOut();
	}

	return 1;
}



//�T�u���[�v�͂���܂��邱�ƂȂ�
int CNewtonRaphsonMethodPlus::SubLoopInit(){
	if( FlagJacobSkip == true ){
		SubLoopCounter = NumVariable;//jacob�s��̍\�z���X�L�b�v����ꍇ�͂����Ȃ�Ō�̌v�Z����X�^�[�g
	}else{
		SubLoopCounter = 0;
	}
	return 1;
}
int CNewtonRaphsonMethodPlus::SubLoopCheck(){
	//���[�v�J�E���^�����Ĕ����o��
	if( NumVariable < SubLoopCounter ){
		return 0;
	}
	return 1;
}
int CNewtonRaphsonMethodPlus::SubLoopReinit(){
	//�T�u�J�E���^�[���グ��D
	SubLoopCounter++;
	return 1;
}


int CNewtonRaphsonMethodPlus::SetDelta( double _Delta ){
	Delta = _Delta;
	return 1;
}
int CNewtonRaphsonMethodPlus::SetError( double _Error ){
	Error = _Error;
	return 1;
}
int CNewtonRaphsonMethodPlus::SetMaxLoop( int _MaxLoop ){
	MaxLoop = _MaxLoop;
	return 1;
}
int CNewtonRaphsonMethodPlus::SetJacobStep( int _Step ){
	JacobStep = _Step;
	return 1;
}
int CNewtonRaphsonMethodPlus::EndReport( void ){
	time(&TimeNow);

	std::cout << "MaxError = " << MaxError << std::endl;
	std::cout << "Time = " << (difftime(TimeNow,TimeStart)) << std::endl;

	return 1;
}



CNewtonRaphsonMethodPlus::CNewtonRaphsonMethodPlus(){

	NumVariable = 0;
	FlagConbergence = false;
	FlagEndReport = false;
	FlagAccAuto = false;
	FlagMatrixOut = false;
	MaxLoop = 5;

	Variable.resize(1);
	Variable[0].resize(1);

	VariableMemory.resize(1);
	VariableMemory[0].resize(1);

	ResultMatrix.resize(1);
	ResultMatrix[0].resize(1);

	JacobMatrix.resize(1);
	JacobMatrix[0].resize(1);

	UpperTriangularMatrix.resize(1);
	UpperTriangularMatrix[0].resize(1);

	DiagonalMatrix.resize(1);
	DiagonalMatrix[0].resize(1);


	SetDelta( 1.00001 );
	SetError( 0.000001 );

	SetJacobStep( 1 );
	FlagJacobSkip = false;

	Min = 1.0;
	SparseMin = 1e-8;

	SetPrint( false , false );

	ErrorMax0 = 100001.0;
	ErrorMax1 = 100000.0;
	ErrorDif0 = -1.0;
	ErrorDif1 = -1.0;


}


CNewtonRaphsonMethodPlus::~CNewtonRaphsonMethodPlus(void){

}


int CNewtonRaphsonMethodPlus::StatusOut(){
	std::cout << "Variable[ " << Variable.size() << " ][ " << Variable[0].size() << " ]" << std::endl;
	std::cout << "VariableMemory[ " << VariableMemory.size() << " ][ " << VariableMemory[0].size() << " ]" << std::endl;
	std::cout << "JacobMatrix[ " << JacobMatrix.size() << " ][ " << JacobMatrix[0].size() << " ]" << std::endl;
	std::cout << "UpperTriangularMatrix[ " << UpperTriangularMatrix.size() << " ][ " << UpperTriangularMatrix[0].size() << " ]" << std::endl;
	std::cout << "DiagonalMatrix[ " << DiagonalMatrix.size() << " ][ " << DiagonalMatrix[0].size() << " ]" << std::endl;
	return 1;
}

int CNewtonRaphsonMethodPlus::GaussianElimination(){

	return 1;
}



int CNewtonRaphsonMethodPlus::Print(){
	std::cout << " Main Loop = " << MainLoopCounter << " Acc = " << Acc[0] << " MaxError = " << MaxError << " AveError = " << AveError << std::endl;
	return 1;
}

int CNewtonRaphsonMethodPlus::PrintAll(){
	int i;
	std::cout << " Main Loop = " << MainLoopCounter << std::endl;
	for( i = 0 ; i < NumVariable ; i++ ){
		std::cout << " Variable[" << i << "] = " << Variable[i][NumVariable] << " Error[" << i << "] = " << ResultMatrix[i][NumVariable] << std::endl;
	}
	return 1;
}

void CNewtonRaphsonMethodPlus::irekae(double &x,double &y ){double z=x;x=y;y=z;}
void CNewtonRaphsonMethodPlus::irekae(int    &x,int    &y ){int    z=x;x=y;y=z;}

void CNewtonRaphsonMethodPlus::NR_free(void){
	if(Variable.size()>0){
		std::vector<std::vector<double>>().swap(Variable);
	}
	if(VariableMemory.size()>0){
		std::vector<std::vector<double>>().swap(VariableMemory);
	}
	if(RelError.size()>0){
		std::vector<double>().swap(RelError); 
	}
	if(AbsError.size()>0){
		std::vector<double>().swap(AbsError); 
	}
	if(DeltaValue.size()>0){
	std::vector<double>().swap(DeltaValue); 
	}
	if(ResultMatrix.size()>0){
		std::vector<std::vector<double>>().swap(ResultMatrix);
	}
	if(JacobMatrix.size()>0){	
	std::vector<std::vector<double>>().swap(JacobMatrix);
	}
	if(UpperTriangularMatrix.size()>0){
		std::vector<std::vector<double>>().swap(UpperTriangularMatrix);
	}
	if(DiagonalMatrix.size()>0){
		std::vector<std::vector<double>>().swap(DiagonalMatrix);
	}
	if(UpperC.size()>0){
		std::vector<double>().swap(UpperC);
	}
	if(Acc.size()>0){
		std::vector<double>().swap(Acc);
	}
	if(MaxValue.size()>0){	
		std::vector<double>().swap(MaxValue);
	}
	if(MinValue.size()>0){
		std::vector<double>().swap(MinValue);
	}
	if(MaxDisplacement.size()){	
		std::vector<double>().swap(MaxDisplacement);
	}
	if(Step.size()){
		std::vector<double>().swap(Step);
	}
	if(StepControl.size()){
		std::vector<bool>().swap(StepControl);
	}
	if(MinMaxControl.size()){	
	std::vector<bool>().swap(MinMaxControl);
	}
//		std::cout<<"Newton�̃f�X�g���N�^"<<std::endl;

}
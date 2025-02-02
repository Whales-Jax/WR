//�M������
#include"DXProperty_ver06.h"
#include<math.h>
using namespace std;
#include "CNewtonRaphsonMethod.h"
#include"CFluidParameter.h"
void main()
{
	int n;
	double P;	// ����[kPa]
	double h;	// �G���^���s[kJ/kg]
	double T;	// ���x[��]
	double x;   //�����x
	double rho1;
	double hl2;
	double rho2;
	double A;double UA;
	double mc;
	double mh;
    double ans;
	double Ti,To,Q,cpi,cpo,ho,dt;
	DXProperty ref;					// DXProperty�N���X�̐錾

	CFluidParameter H_t1;

	mc=0.3;
	mh=0.2;
    A=2.0*3.14*0.01*10.0;
    UA=A*1.0;
	// �Z�b�g�A�b�v
	ref.reibai = "R410A.PPF";		// ��}�̐ݒ�
	ref.joutai = "IIR";				// ���Ԃ̐ݒ� ( ���ʂ�IIR, ����NBP )
	ref.LoadDLL("refprop.DLL");		// �ǂݍ���DLL�t�@�C���̐ݒ� ( 2��ވȏ�̗�}�������K�v�ȏꍇ�C���ꂼ��̗�}�ɑ΂���DLL�t�@�C�����K�v )
	ref.setup();					// �Z�b�g�A�b�v

	//�����l
	Ti =10.0;
	To = 10.0;
	Q = 0;
	
	CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethod�N���X�̐錾
	mnm.setup(10);					//�Z�b�g�A�b�v�@�������̓j���[�g���@�̕ϐ��ȏ�ɂ���

	// �����l�̑��
	mnm.setValue( 1 , To );
	mnm.setValue( 2 , Q );

	//mnm.setAcc(0.8);				//�����x���z�̓��́@���Ȃ��Ă���

	mnm.initial();					//�v�Z���J�n����O�ɕK�����������Ă�������
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// ���܂��Ȃ�
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// ���܂��Ȃ�

			//�l�̐ݒ�
			To = mnm.getValue(0);
			Q = mnm.getValue(1);
            
	// �����l�̌v�Z
	//ref.state_ph( 1000, 350 );	// �C�ӏ�Ԋ֐�( ���́C�G���^���s ) -> �����l��Rc�N���X�Ɋi�[�����
	ref.state_tp( 40, 4000 );	// �P����Ԋ֐�( ���x, ���� )		-> �����l��Rc�N���X�Ɋi�[�����
	//ref.sat_p( P );			// �O�a��Ԋ֐�( ���� ) -> �O�a�t�̕����l��Rl�N���X��, �O�a���C�̕����l��Rv�N���X�Ɋi�[�����
	//ref.sat_t( T );			// �O�a��Ԋ֐�( ���x ) -> �O�a�t�̕����l��Rl�N���X��, �O�a���C�̕����l��Rv�N���X�Ɋi�[�����

	// �����l�̌Ăяo��
	hl2 = ref.Rc.h;			//Rc�N���X�̉��x���Ăяo���Ƃ�
	rho1 = ref.Rc.rho;
	cout<< hl2<<endl;
	cout<< rho1<<endl;



	// �����l�̌v�Z
	//ref.state_ph( 1000, 350 );	// �C�ӏ�Ԋ֐�( ���́C�G���^���s ) -> �����l��Rc�N���X�Ɋi�[�����
	ref.state_tp( Ti, 4000 );	// �P����Ԋ֐�( ���x, ���� )		-> �����l��Rc�N���X�Ɋi�[�����
	//ref.sat_p( P );			// �O�a��Ԋ֐�( ���� ) -> �O�a�t�̕����l��Rl�N���X��, �O�a���C�̕����l��Rv�N���X�Ɋi�[�����
	//ref.sat_t( T );			// �O�a��Ԋ֐�( ���x ) -> �O�a�t�̕����l��Rl�N���X��, �O�a���C�̕����l��Rv�N���X�Ɋi�[�����
		        ho = ref.Rv.h;
	// �����l�̌Ăяo��
	hl2 = ref.Rc.h;			//Rc�N���X�̉��x���Ăяo���Ƃ�
	rho2 = ref.Rc.rho;
	cout<< hl2<<endl;
	cout<< rho2<<endl;

			ref.state_ph( 4000, hl2 );
			Ti= ref.Rc.T;
            ref.state_ph( 4000, ho );
			To= ref.Rc.T;
			//�G���[�l(�����Ɏ�����͂���)
			mnm.setError( 0 , Q/mc+hl2 , ho );
		    dt=(40-10)/(40-To);
			mnm.setError( 1 , UA*(To-10)/(log10(dt)) , Q );
			mnm.prt();				//�G���[�\��
			mnm.prt_sum();			//�G���[�̍��v��\��
		}
	}

	cout << "ho = " << ho << endl;
}

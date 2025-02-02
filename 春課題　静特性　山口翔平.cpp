#include <iostream>
#include <fstream>
#include <cmath>
#include "multi_newton.h"
#include "Property.h"
using namespace std;

	
//****************************************************************************
//void�֐��@�e���u�ł̌v�Z����ݒ�
//****************************************************************************


//*********************
//�@������iEVAPORATER)
//*********************


//****
//����
//****

//�����E�E�E���ʁA���́A��G���^���s
//�o���E�E�E���ʁA���́A��G���^���s�i�|�C���^�j
//�[�U��M


 void Eva(double Eva_Rin_G, double Eva_Rin_P, double Eva_Rin_h,
		 double *Eva_Rout_G, double *Eva_Rout_P, double *Eva_Rout_h, double *M){


//*******
//�ϐ��ݒ�
//*******

	//�M������
	double Eva_Ref_Q;

	//���̔�G���^���s
	double Eva_Ref_h_if;  
	
	//�o���ɂ������}�[�U��	
	double Eva_Rout_M = 0;  


	//�z��ݒ�@��������20�̕��z�萔�n��p���邽�߁A�z��̗v�f�̐���20�Ƃ���.
	//����
	double Eva_Ref_G[21];
	
	//��G���^���s
	double Eva_Ref_h[21];
	
	//���x
	double Eva_Ref_T[21];
	
	//���x
	double Eva_Ref_p[21];

	//�����������C���x�i12��=285K)
	double Eva_Air_T = 12.0;


//******
//�݌v�l
//******

	//����
//	double Eva_ele_t = 0.0;

	//�~������
	double Eva_ele_pi = atan(1.0) * 4.0;
	
	//�`�M�ǒ���
	double Eva_ele_L = 1.0;
	
	//�`�M�ǂ�20�����������̒���
	double Eva_ele_dx = Eva_ele_L/20;

	//�`�M�ǒ��a
	double Eva_ele_D = 0.01;

	//������M�ʉߗ�
	double Eva_ele_K = 13.2;

	//�`�M�ǒf�ʐ�
	double Eva_ele_a = Eva_ele_pi * Eva_ele_D * Eva_ele_D / 4;	


//******	
//�����l
//******

	//Fluid�֐����FEva_r
	Fluid Eva_r;

	//R134A�̕����l�Ăяo��
	Property ref("refprop00.DLL");  
	ref.setup("R134A.FLD","IIR");  


//******
//�v�Z��
//******
	
	//���������
	//0�Ԗڂ̔�G���^���s��void�֐��̈���Eva_Rin_h����.
	Eva_Ref_h[0] = Eva_Rin_h;

	//0�Ԗڂ̗��ʂ�void�֐��̈���Eva_Rin_G����.
	Eva_Ref_G[0] = Eva_Rin_G;
	

	//���ϐ��j���[�g���@
	for(int i=0 ; i<20 ; i++){
				multi_newton mnm;
				mnm.setup(10);

				//�����l�̐ݒ�  �����̔�G���^���s�ɑ�����鎖�ƂȂ�.
				mnm.set_assumed_value( 0 , 300 );

				mnm.initial();

				for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){
					for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){

					//�����x���z
					//mnm.acc = 0.8;

					//�����x���z�����ݒ�i�������₷�����x���j
					//mnm.spcacc();
					//�ϐ��ɏ����l����

					//���̔�G���^���s�ɏ����l����.
					Eva_Ref_h_if = mnm.set_value[0];

					//���ʈ��̎�
					//(i�Ԗڂ̗���)��(i+1�Ԗڂ̗���)
					Eva_Ref_G[i+1] = Eva_Ref_G[i];

					//Eva_r�֐��̌Ăяo��
					//���́F�������́Ai+1�Ԗڂ̉��̔�G���^���s
					ref.state_ph(Eva_Rin_P, Eva_Ref_h_if, Eva_r);
					
					//Eva_r�֐��̏o�͌��ʂ��Ai+1�Ԗڂ̉��x�Ɩ��x�Ƃ���.
					Eva_Ref_T[i+1] = Eva_r.T;
					Eva_Ref_p[i+1] = Eva_r.rho;

					//i+1�Ԗڂ̏ꏊ�ɂ�����A��C�Ɨ�}�̔M������
					Eva_Ref_Q = Eva_ele_K * Eva_ele_pi * Eva_ele_D * 
								Eva_ele_dx*(Eva_Ref_T[i+1] - Eva_Air_T); //<0

					//i+1�Ԗڂ̔�G���^���s
					Eva_Ref_h[i+1]=(Eva_Ref_G[i]*Eva_Ref_h[i]-Eva_Ref_Q)/Eva_Ref_G[i+1];

					//�G���[�l
					//���̔�G���^���s��i+1�Ԗڂ̃G���^���s���r
					mnm.set_error2( 0 , Eva_Ref_h[i+1], Eva_Ref_h_if );
					
					//�G���[�\��
					//mnm.prt();

					//�G���[�̍��v��\��
					//mnm.prt_sum();
					}
				}

				//������o���ɂ������}�[�U��
				Eva_Rout_M = Eva_Rout_M + Eva_Ref_p[i+1] * Eva_ele_pi * Eva_ele_D * Eva_ele_D * Eva_ele_dx/4;
	}

	Eva_Rin_G = Eva_Ref_G[0];
	Eva_Rin_h = Eva_Ref_h[0];

	//�|�C���^
	//*Eva_Rout_G��Eva_Ref_G[20]�����[.
	*Eva_Rout_G = Eva_Ref_G[20];

	//*Eva_Rout_P��Eva_Rin_P�����[.
	*Eva_Rout_P = Eva_Rin_P;

	//*Eva_Rout_h��Eva_Ref_h[20]�����[.
	*Eva_Rout_h = Eva_Ref_h[20];

	//*M��Eva_Rout_M�����[.
	*M = Eva_Rout_M;

	//������o����G�Ƃ�����@�����o���ɂ�����T,P,p���Ăяo��.
	cout << "Eva has caliculated" << endl;

	ref.state_ph(Eva_Rin_P , Eva_Rin_h , Eva_r);
	cout << "����T = " << Eva_r.T<<endl;

//�o��P,h����T�����߂�.
	ref.state_ph(*Eva_Rout_P , *Eva_Rout_h , Eva_r);
	cout << "�o��T = " << Eva_r.T << " " << "�o��P = " << Eva_r.P << " " << "�o��rho = " << Eva_r.rho << endl;
	cout << "�o��h = " << Eva_r.h << " " << "�o��s = " << Eva_r.s << endl;
	cout << "�o��G = " << *Eva_Rout_G << " " << "�o��M = " << Eva_Rout_M << endl;
}


//*********************
//�A���k�@�iCOMPRESSOR)
//*********************

void Com(double Com_Rin_G, double Com_Rin_P, double Com_Rin_h,
		 double *Com_Rout_h , double *Com_Rout_P, double *Com_Rout_G, double *Com_Rout_M ){

//*******
//�ϐ��ݒ�
//*******

	//���� ��G���g���s�[
	double Com_Rin_s ;
	
	//�o�� ��G���g���s�[
	double Com_Rout_s ;
	
	//�o�� ��G���^���s�i�f�M�ߒ��j
	double Com_Rout_had ;

	//���ʈ��
	*Com_Rout_G = Com_Rin_G ;
	
	//�݌v�l
	double Com_ele_W = 0.45;
	double Com_ele_V = 0.00002;
	double Com_ele_eff = 1.0;
	Fluid Com_r;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");
	
	multi_newton mnm;
	mnm.setup(10);

	//�o�����͂����ɐݒ�
	mnm.set_assumed_value( 0 , 1000 );

	mnm.initial();
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){

		//�����x���z
		//mnm.acc = 0.8;

		//�����x���z�����ݒ�i�������₷�����x���j
		//mnm.spcacc();

		//*�o�����͂ɏ����l���
		*Com_Rout_P = mnm.set_value[0];
		
		//Com_r�֐��̌Ăяo��
		//P��h����s�����߂�.(at�����j
		ref.state_ph(Com_Rin_P , Com_Rin_h , Com_r);
		Com_Rin_s = Com_r.s;

		//*�o���G���^���s
		*Com_Rout_h = (Com_Rin_G * Com_Rin_h + Com_ele_W)/(*Com_Rout_G);
		
		//�f�M�������l�������o���G���^���s
		Com_Rout_had = Com_ele_eff * (*Com_Rout_h - Com_Rin_h) +Com_Rin_h;

		//P��had����s�����߂�.(at�o��)
		ref.state_ph(*Com_Rout_P , Com_Rout_had , Com_r);
		Com_Rout_s = Com_r.s;

		//�G���[�l
		//�����Əo���̃G���g���s�[��r.
		mnm.set_error2(0 , Com_Rin_s , Com_Rout_s );
		}
	} 

	//*Rout_P��*Rout_h����T,P,rho,h,s�����߂�.
	ref.state_ph(*Com_Rout_P , *Com_Rout_h , Com_r);
	cout << "Com has caliculated�i�S�ďo���̒l�j" << endl;
	cout << "T = " << Com_r.T << " " << "P = " << Com_r.P << " " << "rho = " << Com_r.rho << endl;
	cout << "h = " << Com_r.h << " " << "s = " << Com_r.s << endl;

	//*�o����}�[�U��
	*Com_Rout_M = Com_ele_V * Com_r.rho;
	cout << "G = " << *Com_Rout_G << " " << "M = " << Com_ele_V*Com_r.rho << endl;
}



//*********************
//�B�Ïk��iCONDENSER)
//*********************

void Con(double Con_Rin_h, double Con_Rin_P, double Con_Rin_G,
		 double *Con_Rout_h, double *Con_Rout_P, double *Con_Rout_G, double *M1){

	double Con_Ref_Q;
	double Con_Ref_h_if;
	double Con_Rout_M = 0;

	double Con_Ref_G[21];
	double Con_Ref_h[21];
	double Con_Ref_T[21];
	double Con_Ref_p[21];

	//�Ïk�������C���x(32.0��=305.0K)
	double Con_Air_T = 32.0;

	//�݌v�l
//	double Con_ele_t = 0.0;
	double Con_ele_pi = atan(1.0) * 4.0;
	double Con_ele_L = 1.0;
	double Con_ele_dx = Con_ele_L/20;
	double Con_ele_D = 0.01;
	double Con_ele_K = 13.3;
	double Con_ele_a = Con_ele_pi * Con_ele_D * Con_ele_D / 4;
	
	//�����l
	Fluid Con_r;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");

	//h,G�̐ݒ�at�Ïk�����
	Con_Ref_h[0] = Con_Rin_h;
	Con_Ref_G[0] = Con_Rin_G;

	Con_Rout_M = 0;

			for(int i=0 ; i<20 ; i++){
				multi_newton mnm;
				mnm.setup(10);
				//�����l�̐ݒ�
				mnm.set_assumed_value( 0 , 300 );

				mnm.initial();

				for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){
					for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){

					//�����x���z
					//mnm.acc = 0.8;

					//�����x���z�����ݒ�i�������₷�����x���j
					//mnm.spcacc();

					//�ϐ��ɏ����l����
					//���̃G���^���s�ɏ����l���.
					Con_Ref_h_if = mnm.set_value[0];
					
					//���ʈ��
					Con_Ref_G[i+1] = Con_Ref_G[i];
					
					//����P�Ɖ��G���^���s����i+1�Ԗڂ�T��p�����߂�.
					ref.state_ph(Con_Rin_P, Con_Ref_h_if, Con_r);
					Con_Ref_T[i+1] = Con_r.T;
					Con_Ref_p[i+1] = Con_r.rho;

					//�����M��
					Con_Ref_Q = Con_ele_K * Con_ele_pi * Con_ele_D * 
								Con_ele_dx*(Con_Ref_T[i+1] - Con_Air_T);	//<0

					//��}�G���^���sat(i+1)�Ԗ�
					Con_Ref_h[i+1]=(Con_Ref_G[i]*Con_Ref_h[i]-Con_Ref_Q)/Con_Ref_G[i+1];

					//�G���[�l
					//i+1�ԖڃG���^���s�Ɖ��G���^���s��r
					mnm.set_error2( 0 , Con_Ref_h[i+1], Con_Ref_h_if );
					
					//�G���[�\��
					//mnm.prt();

					//�G���[�̍��v��\��
					//mnm.prt_sum();
					}
				}
				
				//�o����}�[�U��
				Con_Rout_M = Con_Rout_M + Con_Ref_p[i+1] * Con_ele_pi * Con_ele_D * 
								Con_ele_D * Con_ele_dx/4;
			}

	Con_Rin_G = Con_Ref_G[0];
	Con_Rin_h = Con_Ref_h[0];

	//�o���ɂ������ԗʂ�20�Ԗڂ̒l�����ꂼ����.
	*Con_Rout_G = Con_Ref_G[20];
	*Con_Rout_P = Con_Rin_P;
	*Con_Rout_h = Con_Ref_h[20];
	*M1 = Con_Rout_M;

	//�o��P��h����T,rho,s,G,M�����߂�.


	//�v�Z���ʏo��
	cout << "Con has caliculated" << endl;
	ref.state_ph(Con_Rin_P , Con_Rin_h , Con_r);
	cout << "����T = " << Con_r.T << endl;
	
	ref.state_ph(*Con_Rout_P , *Con_Rout_h , Con_r);
	cout << "�o��T = " << Con_r.T << " " << "�o��P = " << Con_r.P << " " << "�o��rho = " << Con_r.rho << endl;
	cout << "�o��h = " << Con_r.h << " " << "�o��s = " << Con_r.s << endl;
	cout << "�o��G = " << *Con_Rout_G << " " << "�o��M = " << Con_Rout_M << endl;	
}


//*******************************
//�C�A�L�������[�^�[�iACCUMULATOR)
//*******************************

void Acc(double Acc_Vin_P, double Acc_Vin_G, double Acc_Lin_L, double *Acc_Vin_h,
		 double *Acc_Vout_G, double *Acc_Vout_P, double *Acc_Vout_h, double *Acc_Lout_M, double *Acc_Vout_M){

	//����
//	double Acc_Vin_G_int;
//	double Acc_Vin_P_int;
	double Acc_Vin_h_int;

	double Acc_Lout_rhol = 0;
	double Acc_Vout_rhov = 0;

	//�݌v�l
	double Acc_ele_a = 0.01 ;
	double Acc_ele_L = 1.0;
	double Acc_ele_g = 9.80665 ;

	Fluid Acc_rv;
	Fluid Acc_rl;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");


	ref.sat_p(Acc_Vin_P , Acc_rl , Acc_rv);

	//����P�������h���Z�o.
	*Acc_Vin_h = Acc_rv.h;

	//�������C��G���^���s
	Acc_Vin_h_int = *Acc_Vin_h;

	//�o��h=����h
	*Acc_Vout_h = Acc_Vin_h_int;

	//���ʈ��
	*Acc_Vout_G = Acc_Vin_G;

	ref.sat_p(Acc_Vin_P , Acc_rl , Acc_rv);
	Acc_Vout_rhov = Acc_rv.rho;
	Acc_Lout_rhol = Acc_rl.rho;

//	*Acc_Lout_P = Acc_Vin_P_int + (Acc_Vout_rhol) / 1000 * Acc_ele_g * Acc_Lin_L;
	*Acc_Lout_M = Acc_ele_a * Acc_Lin_L * Acc_Lout_rhol;
	*Acc_Vout_M = Acc_ele_a * (Acc_ele_L -  Acc_Lin_L) * Acc_Vout_rhov;
	*Acc_Vout_P = Acc_Vin_P;

	cout << "Acc has caliculated�i�S�ďo���̒l�j" << endl;
	cout << "G = " << *Acc_Vout_G << " " << "P = " << *Acc_Vout_P << endl;
	cout << "h = " << *Acc_Vin_h << " " << "M = " << *Acc_Lout_M + *Acc_Vout_M << endl;
}


//**************************
//�D�c���فiEXPANSION VALVE)
//**************************

void Exp(double Exp_Rin_h, double Exp_Rin_P, double Exp_Rin_G,
		 double *Exp_Rout_h, double *Exp_Rout_P, double *Exp_Rout_G){
	
	double Exp_Rin_p ;

	//�݌v�l
	double Exp_ele_Cv = 0.6;
	double Exp_ele_a = 0.0000269;     //2.69 * 1e-5;

	Fluid Exp_r;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");

	//�v�Z
	ref.state_ph(Exp_Rin_P , Exp_Rin_h , Exp_r);
	Exp_Rin_p = Exp_r.rho;

	//���ʈ��
	*Exp_Rout_G = Exp_Rin_G;
	
	*Exp_Rout_h = Exp_Rin_h;
	*Exp_Rout_P = Exp_Rin_P-(Exp_Rin_G*Exp_Rin_G)/(2 * Exp_Rin_p * Exp_ele_Cv * Exp_ele_Cv * Exp_ele_a * Exp_ele_a); 

	ref.state_ph(*Exp_Rout_P , *Exp_Rout_h , Exp_r);
	cout << "Exp has caliculated�i�o���̒l�j" << endl;
	cout << "T = " << Exp_r.T << " " << "P = " << Exp_r.P << " " << "rho = " << Exp_r.rho << endl;
	cout << "h = " << Exp_r.h << " " << "s = " << Exp_r.s << endl;
	cout << "G = " << *Exp_Rout_G << endl;
}


//****************************************************************************
//main�֐��@void�֐��ō�����e���u�̌v�Z�������s.�i�����ŏ����l����͒l�Ƃ��Ďg�p�j
//****************************************************************************

int main(){

	Fluid Eva_r;
	Fluid Con_r;
	Property ref("refprop00.DLL");
	ref.setup("R134A.FLD","IIR");

	double Eva_Rin_G, Eva_Rin_P, Eva_Rin_h ;
	double Com_Rin_G, Com_Rin_P, Com_Rin_h ;
	double Con_Rin_G, Con_Rin_P, Con_Rin_h ;
	double Acc_Vin_G, Acc_Vin_P, Acc_Lin_L ;
	double Exp_Rin_G, Exp_Rin_P, Exp_Rin_h ;
	

	double ALL_M;
	double Eva_Rout_G = 0;
	double Eva_Rout_P = 0;
	double Eva_Rout_h = 0;
	double Eva_Rout_M = 0;
	double Com_Rout_G = 0;
	double Com_Rout_P = 0;
	double Com_Rout_h = 0;
	double Com_Rout_M = 0;
	double Con_Rout_G = 0;
	double Con_Rout_P = 0;
	double Con_Rout_h = 0;
	double Con_Rout_M = 0;
	double Acc_Vout_G = 0;
	double Acc_Vout_P = 0;
	double Acc_Vout_h = 0;
	double Acc_Vin_h = 0;
	double Acc_Lout_M = 0;
	double Acc_Vout_M = 0;
	double Exp_Rout_G = 0;
	double Exp_Rout_P = 0;
	double Exp_Rout_h = 0;
				
			cout << "SYSTEM �Ó���" << endl;

			//�Ó���
			multi_newton mnm;
			mnm.setup(20);
			//�����l�̐ݒ�
			mnm.set_assumed_value( 0  , 0.02 );
			mnm.set_assumed_value( 1  , 350.0 );
			mnm.set_assumed_value( 2  , 250.0 );
			mnm.set_assumed_value( 3  , 0.02 );
			mnm.set_assumed_value( 4  , 1000.0 );
			mnm.set_assumed_value( 5  , 425.0 );
			mnm.set_assumed_value( 6  , 0.02 );
			mnm.set_assumed_value( 7  , 1000.0 );
			mnm.set_assumed_value( 8  , 250.0 );
			mnm.set_assumed_value( 9  , 0.02 );
			mnm.set_assumed_value( 10 , 350.0 );
			mnm.set_assumed_value( 11 , 400.0 );
			mnm.set_assumed_value( 12 , 0.02 );
			mnm.set_assumed_value( 13 , 350.0 );
			mnm.set_assumed_value( 14 , 0.1 );
			mnm.initial();

			for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){
				for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){

				//�����x���z
				//mnm2.acc = 0.8;

				//�����x���z�����ݒ�i�������₷�����x���j
				//mnm2.spcacc();

				//�ϐ��ɏ����l����
				Eva_Rin_G = mnm.set_value[0];
				Eva_Rin_P = mnm.set_value[1];
				Eva_Rin_h = mnm.set_value[2];
				Con_Rin_G = mnm.set_value[3];
				Con_Rin_P = mnm.set_value[4];
				Con_Rin_h = mnm.set_value[5];
				Exp_Rin_G = mnm.set_value[6];
				Exp_Rin_P = mnm.set_value[7];
				Exp_Rin_h = mnm.set_value[8];
				Com_Rin_G = mnm.set_value[9];
				Com_Rin_P = mnm.set_value[10];
				Com_Rin_h = mnm.set_value[11];
				Acc_Vin_G = mnm.set_value[12];
				Acc_Vin_P = mnm.set_value[13];
				Acc_Lin_L = mnm.set_value[14];


				//�v�Z�@��X�ŕϐ�X�̃A�h���X�������A&X��ǂݍ��ގ��ŕϐ�X���Ăяo����.

				Eva(Eva_Rin_G, Eva_Rin_P, Eva_Rin_h, &Eva_Rout_G, &Eva_Rout_P, 
						&Eva_Rout_h, &Eva_Rout_M);
				Acc(Acc_Vin_P, Acc_Vin_G, Acc_Lin_L, &Acc_Vin_h, &Acc_Vout_G, 
						&Acc_Vout_P, &Acc_Vout_h, &Acc_Lout_M, &Acc_Vout_M);
				Com(Com_Rin_G, Com_Rin_P, Com_Rin_h, &Com_Rout_h, &Com_Rout_P, 
						&Com_Rout_G, &Com_Rout_M);
				Con(Con_Rin_h, Con_Rin_P, Con_Rin_G, &Con_Rout_h, &Con_Rout_P, 
						&Con_Rout_G, &Con_Rout_M);
				Exp(Exp_Rin_h, Exp_Rin_P, Exp_Rin_G, &Exp_Rout_h, &Exp_Rout_P, 
						&Exp_Rout_G);
				ALL_M = Eva_Rout_M + Con_Rout_M + Com_Rout_M + Acc_Lout_M + Acc_Vout_M;

				cout << "ALL_M =0." << ALL_M << endl;
	cout << "," << endl;

				//�G���[�l
				mnm.set_error2( 0  , ALL_M , 1.5);
				mnm.set_error2(  1 , Con_Rout_P , Exp_Rin_P);
				mnm.set_error2(  2 , Con_Rout_h , Exp_Rin_h);
				mnm.set_error2(  3 , Exp_Rout_G , Eva_Rin_G);
				mnm.set_error2(  4 , Exp_Rout_P , Eva_Rin_P);
				mnm.set_error2(  5 , Exp_Rout_h , Eva_Rin_h);
				mnm.set_error2(  6 , Eva_Rout_G , Acc_Vin_G);
				mnm.set_error2(  7 , Eva_Rout_P , Acc_Vin_P);
				mnm.set_error2(  8 , Eva_Rout_h , Acc_Vin_h);
				mnm.set_error2(  9 , Acc_Vout_G , Com_Rin_G);
				mnm.set_error2( 10 , Acc_Vout_P , Com_Rin_P);
				mnm.set_error2( 11 , Acc_Vout_h , Com_Rin_h);
				mnm.set_error2( 12 , Com_Rout_G , Con_Rin_G);
				mnm.set_error2( 13 , Com_Rout_P , Con_Rin_P);
				mnm.set_error2( 14 , Com_Rout_h , Con_Rin_h);
				
				//�G���[�\��
				//mnm.prt();
				//�G���[�̍��v��\��
				mnm.prt_sum();
				}
			}

//**************
//�v�Z���ʂ̏o��
//**************

			cout << "\t" << "G" << "\t" << "\t" << "P" << "\t" << "\t" << "h" << "\t" << endl;

			cout << "EVA" << "\t" << Eva_Rin_G << "\t" << Eva_Rin_P << "\t" <<Eva_Rin_h << "\t" << endl;

			cout << "COM" << "\t" << Com_Rin_G << "\t" << Com_Rin_P << "\t" << Com_Rin_h << "\t" << endl;

			cout << "CON" << "\t" << Con_Rin_G << "\t" << Con_Rin_P << "\t" << Con_Rin_h << "\t" << endl;

			cout << "Acc" << "\t" << Acc_Vin_G << "\t" << Acc_Vin_P << "\t" <<Acc_Vin_h << "\t" << endl;

			cout << "EXP" << "\t" << Exp_Rin_G << "\t" << Exp_Rin_P << "\t" <<Exp_Rin_h << "\t" << "\t" << endl;
			cout << "," <<endl;
			cout << "EVA" << "\t" << Eva_Rout_G << "\t" << Eva_Rout_P << "\t" <<Eva_Rout_h << "\t" << endl;

			cout << "COM" << "\t" << Com_Rout_G << "\t" << Com_Rout_P << "\t" <<Com_Rout_h << "\t" << endl;

			cout << "CON" << "\t" << Con_Rout_G << "\t" << Con_Rout_P << "\t" <<Con_Rout_h << "\t" << endl;

			cout << "Acc" << "\t" << Acc_Vout_G << "\t" << Acc_Vout_P << "\t" <<Acc_Vout_h << "\t" << endl;

			cout << "EXP" << "\t" << Exp_Rout_G << "\t" << Exp_Rout_P << "\t" <<Exp_Rout_h << "\t" << endl;
			cout << "ALL_M" << "\t" << ALL_M << endl;
			cout << "Acc_ele.L = " << Acc_Lin_L << endl;

			ref.state_ph(Eva_Rin_P , Eva_Rin_h , Eva_r);
			ref.state_ph(Con_Rin_P , Con_Rin_h , Con_r);

			cout << "T_Eva_In = " << Eva_r.T << "\t" << "T_Con_In = " << Con_r.T <<endl;

			ref.state_ph(Eva_Rout_P , Eva_Rout_h , Eva_r);
			ref.state_ph(Con_Rout_P , Con_Rout_h , Con_r);

			cout << "T_Eva_Out = " << Eva_r.T << "\t" << "T_Con_Out = " << Con_r.T <<endl;
		
		

}
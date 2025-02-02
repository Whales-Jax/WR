#include <iostream>
#include <fstream>
#include<math.h>
using namespace std;
#include "CNewtonRaphsonMethod.h"
#define M 15		//��������
#define N 15		//�ƕ����i90�x�܂Łj
#define L 1		//Re�J�n
#define P 60	//Re�I���

int main(){
	

	//�Z�x56%	���x100��(20120926_3)(�e�함����56%,Tv����)
	double u[M][2*N];
	double uo[M][N];
	double T[M][2*N];
	double Tv=97.507325859;		//�t�����C�����x(�Z�x�C���͂��)
	double Ti=83.35 ;		//�`�M�ǉ��x
	double Ta=99.90 ;		//�n�t�������x
	double k=0.476361222;	//�M�`����
	double p=1603.612594 ;		//���x
	double q=0.001645237;	//�S���W��
	double cp=2050.162898;			//��M
	double sig=0.082086053;		//�\�ʒ���
	double c;				//��/��
	c=q/p;
	double a;				//k/��c
	a=k/(p*cp);
	int i=0;
	int j=0;
	int t=0;
	double g=9.80665;
	double r=1.8e-2;	//�ǔ��a
	double rr[N];		//film thickness[m]
	double dr[2*N];		//film thickness per a cell[m]
	double dro[2*N];		//film thickness per a cell[m](ho)
	double uave[N];	
	double uoave[N];
	double sum=0;
	double h[2*N];
	double b[2*N];
	double bave;		//�X������
	double have[P]={0};		//�M�`�B��
	double rrave;		//�t���������ρi�S�p�x�j
	double uuave;		//�������ρi�S�Z���j
	double width[N];

	//contact angle 29.7[deg]
	double hoo=0.5285;	//minimum thickness[-]
	double ho[N]={0};			//
	double WR[N]={0};			//wetting ratio[-]
	double WRave[P]={0};			//average of WR[-]
	double sine[N]={0};
	double aa;
	aa=pow(p,3)*pow(g,2)/(15*pow(q,2)*sig);
	double aaa;

	for(j=0 ; j<N ; j++){			//calculation of ho
		sine[j]=sin(3.14159*(2*j+1)/(4*N));
		aaa=sine[j]*sine[j]*aa;
		ho[j]=hoo/pow(aaa,0.2);
		dro[j]=ho[j]/M;
	}

	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<=N ; j++){
			uo[i][j]=0.03;
		}
	}
	for(j=0 ; j<=N ; j++){
		uoave[j]=0.03;
	}

	//*******************************�@����ho�̎��̑��x�v�Z�@**************************************************
	CNewtonRaphsonMethod mnm;		//**************************���x���z******************************************************
	mnm.setup((M+1)*(N+1));					//�Z�b�g�A�b�v�@�����̓j���[�g���@�̕ϐ��ȏ�ɂ��Ă��������@��2:
	t=0;
	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			mnm.setValue( t , uo[i][j] );
			t++;
		}
	}
	//mnm.setAcc(0.8);				//�����x���z�̓��́@�Ȃ��Ă��@��3:
	mnm.initial();					//�v�Z���J�n����O�ɕK�����������Ă�������
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		//���܂��Ȃ�
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	//���܂��Ȃ�
			
			t=0;
			for(i=0 ; i<=M ; i++){
				for(j=0 ; j<N ; j++){
					uo[i][j] = mnm.getValue(t);
					t++;
				}
			}
		
			sum=0;
			for(j=0 ; j<N ; j++){
				for(i=1 ; i<=M ; i++){
					sum=sum+uo[i][j];
				}
				uoave[j]=sum/M;
				sum=0;
			}

			//�G���[�l
			t=0;
			for(j=0 ; j<N ; j++){		//���C��
				mnm.setError( t , uo[0][j] , uo[1][j] );
				t++;
			}

			for(i=1 ; i<M ; i++){
				for(j=0 ; j<N ; j++){
					mnm.setError( t , g*sin(3.14159*(2*j+1)/(4*N))+c*(uo[i+1][j]-2*uo[i][j]+uo[i-1][j])/(dro[j]*dro[j]) , 0 );
					t++;
				}
			}

			for(j=0 ; j<N ; j++){		//�Ǒ�
				mnm.setError( t , uo[M][j] , 0 );
				t++;
			}	



			//mnm.prt();				//�G���[�\��
			mnm.prt_sum();			//�G���[�̍��v��\��
		}
	}


	//�����l
	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<=2*N ; j++){
			u[i][j] = 0.02;
		}
	}
	for(j=0 ; j<=N ; j++){
		dr[j] = 0.001/N;
		uave[j]=0.02;
	}

		//�����l
	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<2*N ; j++){
			T[i][j] = 90.0;	
		}
	}

	for(int z=L; z<=P; z++){	////////////////////////////////�@Re���[�v�J�n///////////////////////
	double m=z*q/4 ;
	//double m=(16+0.1*z)*q/4;
	cout<<"Re="<<z<<"-------------------------------------------"<<endl;

	
	for(j=0 ; j<N ; j++){			//calculation of ho
		sine[j]=sin(3.14159*(2*j+1)/(4*N));
		aaa=sine[j]*sine[j]*aa;
		ho[j]=hoo/pow(aaa,0.2);
		dro[j]=ho[j]/M;
	}

	
	CNewtonRaphsonMethod mnm2;		//**************************���x���z******************************************************
	mnm2.setup((M+2)*(N+2));					//�Z�b�g�A�b�v�@�����̓j���[�g���@�̕ϐ��ȏ�ɂ��Ă��������@��2:
	t=0;
	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			mnm2.setValue( t , u[i][j] );
			t++;
		}
	}


	//mnm.setAcc(0.8);				//�����x���z�̓��́@�Ȃ��Ă��@��3:

	mnm2.initial();					//�v�Z���J�n����O�ɕK�����������Ă�������
	for(mnm2.main_loop_init();mnm2.main_loop_check();mnm2.main_loop_reinit()){		//���܂��Ȃ�
		for(mnm2.sub_loop_init();mnm2.sub_loop_check();mnm2.sub_loop_reinit()){	//���܂��Ȃ�

			t=0;
			for(i=0 ; i<=M ; i++){
				for(j=0 ; j<N ; j++){
					u[i][j] = mnm2.getValue(t);
					t++;
				}
			}
			
			//�G���[�l
			
			for(j=0 ; j<N ; j++){		//���ϗ���
				for(i=1 ; i<=M ; i++){
					sum=sum+u[i][j];
				}
				uave[j]=sum/M;
				sum=0;
			}
			for(j=0 ; j<N ; j++){		//�t������
				dr[j] = m/(M*p* uave[j]);
			}


			t=0;
			for(j=0 ; j<N ; j++){		//���C��
				mnm2.setError( t , u[0][j] , u[1][j] );
				t++;
			}

			for(i=1 ; i<M ; i++){
				for(j=0 ; j<N ; j++){
					mnm2.setError( t , g*sin(3.14159*(2*j+1)/(4*N))+c*(u[i+1][j]-2*u[i][j]+u[i-1][j])/(dr[j]*dr[j]) , 0 );
					t++;
				}
			}

			for(j=0 ; j<N ; j++){		//�Ǒ�
				mnm2.setError( t , u[M][j] , 0.0 );
				t++;
			}		



			
			

			//mnm2.prt();				//�G���[�\��
			mnm2.prt_sum();			//�G���[�̍��v��\��
		}
	}//***********************************************************************************************************************
	
	for(i=0; i<=M ; i++){		//90-180�x��0-90�x�̑��x������
		for(j=N ; j<=2*N ; j++){
			u[i][j]=u[i][2*N-1-j];
		}
	}



	
	


	
	CNewtonRaphsonMethod mnm3;		//*******************************���x���z**************************************************
	mnm3.setup((M+1)*(2*N+1));					//�Z�b�g�A�b�v�@�����̓j���[�g���@�̕ϐ��ȏ�ɂ��Ă��������@��2:
	t=0;
	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<2*N ; j++){
			mnm3.setValue( t , T[i][j] );
			t++;
		}
	}
	//mnm3.setAcc(0.8);				//�����x���z�̓��́@�Ȃ��Ă��@��3:

	mnm3.initial();					//�v�Z���J�n����O�ɕK�����������Ă�������
	for(mnm3.main_loop_init();mnm3.main_loop_check();mnm3.main_loop_reinit()){		//���܂��Ȃ�
		for(mnm3.sub_loop_init();mnm3.sub_loop_check();mnm3.sub_loop_reinit()){	//���܂��Ȃ�
		
			t=0;
			for(i=0 ; i<=M ; i++){
				for(j=0 ; j<2*N ; j++){
					T[i][j] = mnm3.getValue(t);
					t++;
				}
			}
		
			t=0;
			//�G���[�l
			for(j=0 ; j<2*N ; j++){		//���C��
				mnm3.setError( t , T[0][j] , Tv );
				t++;
			}
			
			for(i=1 ; i<M ; i++){		
				mnm3.setError( t , T[i][0] , Ta );
				t++;
			}
		
			for(i=1 ; i<M ; i++){
				for(j=1 ; j<2*N ; j++){
					mnm3.setError( t ,(u[i][j]*(T[i][j]-T[i][j-1]))/(r*(3.14159/M)), a*((T[i+1][j]-2*T[i][j]+T[i-1][j])/(dr[j]*dr[j])+(T[i][j+1]-2*T[i][j]+T[i][j-1])/(r*r*(3.14159/M)*(3.14159/M))) );
					t++;
				}
			}
			
			for(j=0 ; j<2*N ; j++){		//�Ǒ�
				mnm3.setError( t , T[M][j] , Ti );
				t++;
			}
			
			//for(i=1 ; i<M ; i++){		//
			//	mnm3.setError( t , T[i][N-1] , Ti+(Tv-Ti)*(M-i)/M);
			//	t++;
			//}

			//mnm3.prt();				//�G���[�\��
			mnm3.prt_sum();			//�G���[�̍��v��\��
		}
	}
	
	/*
	sum=0;
	//fout << endl << "����" <<endl;
	for(j=0 ; j<N ; j++){
		sum=sum+uave[j];
	}
	uuave=sum/N;
	//fout << endl << "�t������" << endl;
	for(j=0 ; j<N ; j++){
		rr[j]=dr[j]*M;
		//fout << rr[j] << "," << endl;
	}*/


	//*******************************�@WR�v�Z�@**************************************************
	for(j=0 ; j<N ; j++){
		//cout<<"ho"<<ho[j]<<endl;
		//cout<<"dro"<<dro[j]<<endl;
		if(dr[j]>=dro[j]){
			WR[j]=1.0;
		}else{
			//cout<<"uave"<<uave[j]<<endl;
			//cout<<"uoave"<<uoave[j]<<endl;
			width[j]=m/(uoave[j]*p*ho[j]);
			//cout<<"width"<<width[j]<<endl;
			WR[j]=0.4*width[j];
			dr[j]=dro[j];
			//cout<<"width"<<width[j]<<endl;
			//cout<<"uoave"<<uoave[j]<<endl;
			//cout<<"ho"<<ho[j]<<endl;
		}
		//cout<<"WR,"<<WR[j]<<endl;
	}	
	for(j=N ; j<=2*N ; j++){	//90-180�x��0-90�x��WR,dr,������
		WR[j]=WR[2*N-1-j];
		dr[j]=dr[2*N-1-j];
	}


	sum=0;
	for(j=0 ; j<N ; j++){
		sum=sum+WR[j];
	}
	WRave[z]=sum/N;
	//cout<<"WRave["<<z<<"]"<<WRave[z]<<endl;

	//fout << endl << "�X��" << endl;
	for(j=0 ; j<2*N ; j++){
		b[j]=(T[M-1][j]-T[M][j])/dr[j];
		//cout<<"b"<<b[j]<<endl;
		h[j]=k/(Tv-Ti)*b[j]*WR[j];
		//cout<<T[M-1][j]<<"-"<<T[M][j]<<"/"<<dr[j]<<"="<<b[j]<<endl;
		//cout<<k<<"*"<<b[j]<<"*"<<WR[j]<<"/"<<Tv-Ti<<"="<<h[j]<<endl;
		//fout<< b[j] << endl;
		//cout<<"h"<<h[j]<<endl;
		//cout<<h[j]<<"*"<<WR[j]<<endl<<endl;
	}

	sum=0;
	for(j=0 ; j<2*N ; j++){
		sum=sum+h[j];
	}
	have[z]=sum/N;
	//cout<<"WRave"<<WRave[z]<<endl;
	//cout<<"have"<<have[z]<<endl;




	}	//*******************************************Re���[�v�����********************************

	ofstream fout("a.csv");
	if(!fout){
	cout << "�t�@�C�����I�[�v���ł��܂���ł���"<< endl;
    return 1;
	}
	else
	cout<< "�t�@�C�����I�[�v�����܂���" << endl;
	
	
	fout << "���x���z" << endl <<",";
	for(j=1 ; j<=2*N ; j++){
		fout<< j << "," ;
	}
	fout << endl;
	for(i=0 ; i<M ; i++){
		fout << i <<",";
		for(j=0 ; j<2*N ; j++){
			fout << u[i][j] << ",";
		}
		fout << endl;
	}

	fout << endl << "���x���z" << endl<< ",";
	for(j=0 ; j<2*N ; j++){
		fout<< j << "," ;
	}
	fout << endl;
	for(i=0 ; i<=M ; i++){
		fout << i <<",";
		for(j=0 ; j<2*N ; j++){
			fout << T[i][j] << ",";
		}
		fout << endl;
	}

	/*
	sum=0;
	for(j=0 ; j<N ; j++){
		sum=sum+b[j];
	}
	bave=sum/N;
	sum=0;
	for(j=0 ; j<N ; j++){
		sum=sum+rr[j];
	}
	rrave=sum/N;*/

	//fout << endl << "�X������,"<<bave;
	//fout << endl << "�t����������,"<<rrave;
	//fout << endl << "��������," <<uuave;
	//fout << endl << "�M�`�B��,"<<have;
	

	/*
	fout<<"WR"<<endl;
	for(j=0 ; j<N ; j++){
		fout<<WR[j]<<endl;
	}*/

	fout << "WR����,�M�`�B��"<<endl;	
		for(int l=L; l<=P; l++){
			fout<<WRave[l]<<","<<have[l]<<endl;
		}

	fout.close();

	return 0;
}
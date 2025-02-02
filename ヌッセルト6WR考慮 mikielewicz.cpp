#include <iostream>
#include <fstream>
#include<math.h>
using namespace std;
#include "CNewtonRaphsonMethod.h"
#define M 10
#define N 10

int main(){

				ofstream fout("�k�b�Z���g6.csv");
	if(!fout){
	cout << "�t�@�C�����I�[�v���ł��܂���ł���"<< endl;
    return 1;
	}
	else
		cout<< "�t�@�C�����I�[�v�����܂���" << endl;
	



	//�Z�x60%	���x��
	double u[M][N];
	double T[M][N];
	double Tv=96.44432489;		//�t�����C�����x
	double Ti=83.35 ;		//�`�M�ǉ��x
	double Ta=99.90 ;		//�n�t�������x
	double k=0.479938777;	//�K���C�ق�Ƃ͔Z�x�C���x�ɂ���ĕς��
	double p=1711.790009 ;		//���x
	double m=0.3 ;
	double q=0.00650693170434067;	//�S���W��
	double cp=2069.194255;			//��M
	double sig=0.090283664901;		//�\�ʒ���
	double c;				//��/��
	c=q/p;
	double a;				//k/��c
	a=k/(p*cp);
	int i,j;
	double g=9.80665;
	double r=1.8e-2;	//�ǔ��a
	double rr[N];		//film thickness[m]
	double dr[N];		//film thickness per a cell[m]
	double uave[N];		
	double sum=0;
	double h[N];
	double b[N];
	double bave;		//�X������
	double have;		//�M�`�B��
	double rrave;		//�t���������ρi�S�p�x�j
	double uuave;		//�������ρi�S�Z���j
	double width[N];

	//contact angle 29.7[deg]
	double hoo=0.5285;	//minimum thickness[-]
	double ho[N];			//
	double WR[N];			//wetting ratio[-]
	double WRave;			//average of WR[-]
	double sine[N];
	double aa;
	aa=pow(p,3)*pow(g,2)/(15*pow(q,2)*sig);
	double aaa;

	for(j=0 ; j<=N ; j++){			//calculation of ho
		sine[j]=sin(3.14159*(2*j+1)/(2*N));
		aaa=sine[j]*sine[j]*aa;
		ho[j]=hoo/pow(aaa,0.2);
	}
	
//�����l
	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<=N ; j++){
			u[i][j] = 0.03;
		}
	}
	for(j=0 ; j<=N ; j++){
		dr[j] = 0.001/N;
		uave[j]=0.03;
	}
	
	CNewtonRaphsonMethod mnm;		//**************************���x���z******************************************************
	mnm.setup((M+2)*(N+2));					//�Z�b�g�A�b�v�@�����̓j���[�g���@�̕ϐ��ȏ�ɂ��Ă��������@��2:

	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			mnm.setValue( j+N*i , u[i][j] );
		}
	}


	//mnm.setAcc(0.8);				//�����x���z�̓��́@�Ȃ��Ă��@��3:

	mnm.initial();					//�v�Z���J�n����O�ɕK�����������Ă�������
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		//���܂��Ȃ�
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	//���܂��Ȃ�

			for(i=0 ; i<=M ; i++){
				for(j=0 ; j<N ; j++){
					u[i][j] = mnm.getValue(j+N*i);
				}
			}
		
			
			//�G���[�l
			for(j=0 ; j<N ; j++){
				for(i=0 ; i<=M ; i++){
					sum=sum+u[i][j];
				}
				uave[j]=sum/(M+1);
				sum=0;
			}
			for(j=0 ; j<N ; j++){		//�t������
				dr[j] = m/(M*p* uave[j]);
			}
			
			
			
			
			for(j=0 ; j<N ; j++){		//���C��
				mnm.setError( j , u[0][j] , u[1][j] );
			}

			for(i=1 ; i<M ; i++){
				for(j=0 ; j<N ; j++){
					mnm.setError( j+N*i , g*sin(3.14159*(2*j+1)/(2*N))+c*(u[i+1][j]-2*u[i][j]+u[i-1][j])/(dr[j]*dr[j]) , 0 );
				}
			}

			for(j=0 ; j<N ; j++){		//�Ǒ�
				mnm.setError( j+N*N , u[M][j] , 0.0 );
			}			

			mnm.prt();				//�G���[�\��
			mnm.prt_sum();			//�G���[�̍��v��\��
		}
	}//***********************************************************************************************************************
	





	fout << "���x���z" << endl <<",";
	for(j=0 ; j<N ; j++){
		fout<< j << "," ;
	}


	fout << endl;
	for(i=0 ; i<=M ; i++){
		fout << i <<",";
		for(j=0 ; j<N ; j++){
			fout << u[i][j] << ",";
		}
		fout << endl;
	}
	


	
	
	//�����l
	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			T[i][j] = 40.0;	
		}
	}

	CNewtonRaphsonMethod mnm2;		//*******************************���x���z**************************************************
	mnm2.setup((M+1)*(N+1));					//�Z�b�g�A�b�v�@�����̓j���[�g���@�̕ϐ��ȏ�ɂ��Ă��������@��2:

	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			mnm2.setValue( j+N*i , T[i][j] );
		}
	}
	
	
	//mnm.setAcc(0.8);				//�����x���z�̓��́@�Ȃ��Ă��@��3:

	mnm2.initial();					//�v�Z���J�n����O�ɕK�����������Ă�������
	for(mnm2.main_loop_init();mnm2.main_loop_check();mnm2.main_loop_reinit()){		//���܂��Ȃ�
		for(mnm2.sub_loop_init();mnm2.sub_loop_check();mnm2.sub_loop_reinit()){	//���܂��Ȃ�
		
			for(i=0 ; i<=M ; i++){
				for(j=0 ; j<N ; j++){
					T[i][j] = mnm2.getValue(j+N*i);
				}
			}
			int t=0;
			//�G���[�l
			for(j=0 ; j<N ; j++){		//���C��
				mnm2.setError( t , T[0][j] , Tv );
				t++;
			}
			
			for(i=1 ; i<M ; i++){		
				mnm2.setError( t , T[i][0] , Ta );
				t++;
			}
		
			for(i=1 ; i<M ; i++){
				for(j=1 ; j<N-1 ; j++){
					mnm2.setError( t ,(uave[j]*(T[i][j+1]-T[i][j-1]))/(r*(3.14159/M)), a*((T[i+1][j]-2*T[i][j]+T[i-1][j])/(dr[j]*dr[j])+(T[i][j+1]-2*T[i][j]+T[i][j-1])/(r*r*(3.14159/M)*(3.14159/M))) );
					t++;
				}
			}
		
			for(j=0 ; j<N ; j++){		//�Ǒ�
				mnm2.setError( t , T[M][j] , Ti );
				t++;
			}
			for(i=1 ; i<M ; i++){		//
				mnm2.setError( t , T[i][N-1] , Ti+(Tv-Ti)*(M-i)/M);
				t++;
			}

			mnm2.prt();				//�G���[�\��
			mnm2.prt_sum();			//�G���[�̍��v��\��
		}
	}



	

	fout << endl << "���x���z" << endl<< ",";
	for(j=0 ; j<N ; j++){
		fout<< j << "," ;
	}
	fout << endl;
	for(i=0 ; i<=M ; i++){
		fout << i <<",";
		for(j=0 ; j<N ; j++){
			fout << T[i][j] << ",";
		}
		fout << endl;
	}

	sum=0;
	fout << endl << "����" <<endl;
	for(j=0 ; j<N ; j++){
		sum=sum+uave[j];
	}
	uuave=sum/N;

	fout << endl << "�t������" << endl;
	for(j=0 ; j<N ; j++){
		rr[j]=dr[j]*M;
		fout << rr[j] << "," << endl;
	}

	


	fout << endl << "�X��" << endl;
	for(j=0 ; j<N ; j++){
		b[j]=(T[M-1][j]-T[M][j])/dr[j];
		h[j]=k/(Tv-Ti)*b[j];
		fout<<T[M-1][j]<<"-"<<T[M][j]<<"/"<<dr[j]<<"=";
		fout<< b[j] << endl;
	}




	sum=0;
	for(j=1 ; j<N ; j++){
		sum=sum+b[j];
	}
	cout<< sum;
	bave=sum/(N-1);

	sum=0;
	for(j=1 ; j<N ; j++){
		sum=sum+h[j];
	}
	cout<< sum;
	have=sum/(N-1);

	sum=0;
	for(j=1 ; j<N ; j++){
		sum=sum+rr[j];
	}
	cout<< sum;
	rrave=sum/(N-1);





	//*******************************�@����ho�̎��̑��x�v�Z�@**************************************************
	CNewtonRaphsonMethod mnm3;		//**************************���x���z******************************************************
	mnm3.setup((M+2)*(N+2));					//�Z�b�g�A�b�v�@�����̓j���[�g���@�̕ϐ��ȏ�ɂ��Ă��������@��2:

	for(i=0 ; i<=M ; i++){
		for(j=0 ; j<N ; j++){
			mnm3.setValue( j+N*i , u[i][j] );
		}
	}
	//mnm.setAcc(0.8);				//�����x���z�̓��́@�Ȃ��Ă��@��3:
	mnm3.initial();					//�v�Z���J�n����O�ɕK�����������Ă�������
	for(mnm3.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		//���܂��Ȃ�
		for(mnm3.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	//���܂��Ȃ�

			for(i=0 ; i<=M ; i++){
				for(j=0 ; j<N ; j++){
					u[i][j] = mnm3.getValue(j+N*i);
				}
			}
		
			
			//�G���[�l
			for(j=0 ; j<N ; j++){
				for(i=0 ; i<=M ; i++){
					sum=sum+u[i][j];
				}
				uave[j]=sum/(M+1);
				sum=0;
			}
			for(j=0 ; j<N ; j++){		//�t������
				dr[j] = ho[j]/M;
			}
			
			
			
			
			for(j=0 ; j<N ; j++){		//���C��
				mnm3.setError( j , u[0][j] , u[1][j] );
			}

			for(i=1 ; i<M ; i++){
				for(j=0 ; j<N ; j++){
					mnm3.setError( j+N*i , g*sin(3.14159*(2*j+1)/(2*N))+c*(u[i+1][j]-2*u[i][j]+u[i-1][j])/(dr[j]*dr[j]) , 0 );
				}
			}

			for(j=0 ; j<N ; j++){		//�Ǒ�
				mnm3.setError( j+N*N , u[M][j] , 0 );
			}			

			mnm3.prt();				//�G���[�\��
			mnm3.prt_sum();			//�G���[�̍��v��\��
		}
	}


	//*******************************�@WR�v�Z�@**************************************************
	for(j=0 ; j<N ; j++){
		cout<<"ho"<<ho[j]<<endl;
		cout<<"r"<<rr[j]<<endl;

		if(rr[j]>=ho[j]){
			WR[j]=1.0;
		}else{
			width[j]=m/(uave[j]*p*ho[j]);
			cout<<"width"<<width[j]<<endl;
			WR[j]=0.4*width[j];
		}
		fout<<"WR,"<<WR[j]<<endl;
	}

	sum=0;
	for(j=1 ; j<N ; j++){
		sum=sum+WR[j];
	}
	WRave=sum/(N-1);

	fout << endl << "�X������,"<<bave;
	fout << endl << "�t����������,"<<rrave;
	fout << endl << "��������," <<uuave;
	fout << endl << "�M�`�B��,"<<have;
	fout << endl << "WR����,"<<WRave;
	fout.close();

	return 0;
}
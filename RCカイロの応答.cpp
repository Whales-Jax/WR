/*
#include<iostream>
#include<fstream>
using namespace std;


int 
{
double x[10000],T,R,C,Ef,dt,dd;
int i,k,t;
x[0]=0;
	ofstream fout("sabunhouback.csv");
	if(!fout){
	cout << "�t�@�C�����I�[�v���ł��܂���ł���"<< endl;
    return 1;
	}
	else
		cout<< "�t�@�C�����I�[�v�����܂���" << endl;

	cout <<"��R�ƃR���f���T,�X�e�b�v���͒l����͂��Ă�������"<< endl;
	cin >> R;
	fout<<"R=,"<< R << endl;
	cin >> C;
	fout<<"C=,"<< C << endl;
	fout<<"T=,"<< R*C << endl;
	cin >> Ef;
		fout<<"Ef=,"<< Ef << endl;
	cout <<"�����@�ł̃�t����͂��Ă�������"<< endl;
	cin >> dt;
	    fout<<"��t=,"<< dt << endl;

		T=R*C;
t=0;
k=0;

for(i=0; i<9000; i++){
	
	x[i+1]=(dt/T*Ef+x[i])/(1+dt/T);
	dd=x[i+1]-x[i];


	if(dd<0.00000000000000001){
		if(k==0){
		t=i+1;
		k=1;
		}
	}
}


for(i=1; i<t; i++){
    C=dt*(i-1)/T;
	cout << "t=," << C <<",x["<< i <<"}=," << x[i]/Ef << endl;
	fout << "t=," << C <<",x["<< i <<"}=," << x[i]/Ef << endl;
}


return 0;
}*/
#include <iostream>
#include <fstream>
using namespace std;

int main(void) {
	double T[10][20];
	int i, j;
	double dt, t;

	t=0;
	dt=0.1;
	
	ofstream fout("cond2.csv");
	if(!fout){
		cout << "ファイルをオープンできませんでした\n";
		return 1;
	}
	else
		cout << "ファイルをオープンしました\n";
	
	for(j=0;j<20;j++) {
		for(i=0;i<=10;i++) {
			T[i][j]=1.0/10.0*(double)i;
			if(j==0)
				fout<<T[i][j] <<",";
		}
	}
	fout <<t <<"\n";

	for(j=1;j<20;j++) {
		T[0][j]=1;
	}

	for(j=1;j<20;j++) {
		fout <<T[0][j] <<",";
		for(i=1;i<=9;i++){
			T[i][j]=0.05*dt*(T[i+1][j-1]-2.0*T[i][j-1]+T[i-1][j-1])*100.0+T[i][j-1];
			fout <<T[i][j] <<",";
		}
		t += dt;
		fout <<T[10][j] <<"," <<t <<"\n";
	}
	


	cout << "ファイルに書き込みました\n";
	fout.close();
	cout << "ファイルに書き込みました\n";

	return 0;
}
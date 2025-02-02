/*///熱交換器
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "CNewtonRaphsonMethod.h"
#include"CFluidParameter.h"
#include"PropertyLithiumBromide_ver05.h"

PropertyLithiumBromide libr;

int main()
{
	int i,j,N,M,O;
	N=60;/////Angle's mesh
	M=40;/////Thickness's mesh
	O=10;/////Mass flow rate
	double st, L,rho,myu,beta,hmin,hminp,Re,h[100],g,CA[100],G,Ti,Xi,PI,X[100],f,m,gl,k,n,ma,Xt[10],a[100],A[100],aave,Reb,aavet,dp,b[100];

	ofstream fout("DynNuS.csv");
	if(!fout){
	cout << "ファイルをオープンできませんでした"<< endl;
    return 1;
	}
	else
		cout<< "ファイルをオープンしました" << endl;


////////////////////Experiment condition////////////////////
	Ti=180.0;
	Xi=61.0/100.0;

	st=libr.sc_st_XT(Xi,Ti);
	rho=libr.sc_rho_XT(Xi,Ti);
	myu=libr.sc_visc_XT(Xi,Ti);	
	gl=9.81;//[m/s2]
	k=libr.sc_thc_XT(Xi,Ti);
    L=0.878;//tube length[m]
	PI=(6*asin(0.5));//π
	beta=PI/N;//Angle of tube//////for 90degree
	dp=(rho*pow(st,3.0))/(pow(myu,4.0)*gl);
	hmin=pow(3.0*(0.2314*pow(dp,0.1761))*pow(myu/rho,2.0)/9.81,1.0/3.0);//minimum thickness [mm] from Maron1982

	Reb=2.0*(0.2314*pow(dp,0.1761));


	ma=0.0005;//[kg/s]

	for(n=1;n<200;n++){///////////////Mass flow rate loop
		m=ma*n;
//		fout<<"mass flow rate,"<<m<<endl;

		aave=0;
		aavet=0;

///////////////Angle loop

	
    Re=2.0*m/L/myu;

	for(i=1;i<N-3;i++){
		if(i<N/2){
		h[i]=pow(3.0*m/L*myu/(rho*rho*gl*sin(beta*i)),1.0/3.0);
//		fout<<h[i]<<"h,";
		aavet=0;
		if(h[i]>=hmin){
			X[i]=1.0;
//			fout<<"alpha for angle,"<<k*X[i]/h[i]<<",";
			aave += k*X[i]/h[i]/N;/////calculation of average of alpha

		}
		else{
			hminp=pow(((rho*rho*rho)*gl*gl*sin(beta*i)*sin(beta*i))/(15.0*myu*myu*st),0.2)*hmin;
//			cout<<hminp<<"hminp"<<endl;

			CA[i]= 1665.4*pow(hminp,6.0) - 3446.8*pow(hminp,5.0) + 2399.3*pow(hminp,4.0) - 459.73*pow(hminp,3.0) + 49.538*pow(hminp,2.0) + 2.9741*hminp + 0.8616;//////////////////////////////////////////////////////
			CA[i]=CA[i]/180*PI;
				f=-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])
		-3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i];//f(θo)
	g=(CA[i]*(5.0/16.0+15.0/4.0*pow(cos(CA[i]),2.0) +5.0/2.0*pow(cos(CA[i]),4.0))
	  -sin(CA[i])*(113.0/48.0*cos(CA[i]) +97.0/24.0*pow(cos(CA[i]),3.0) +1.0/6.0*pow(cos(CA[i]),5.0)));//Ψ(θo))

				X[i]=pow((h[i]),3.0)*sin(CA[i])/f*pow(2.0/45.0*(pow(rho,3.0)
					*pow((9.81*sin(beta*i)),2.0)/(st*pow(myu,2.0)))*g/sin(CA[i])*pow((CA[i]/sin(CA[i])-cos(CA[i])),-1),3.0/5.0);
				b[i]=(h[i]*(CA[i]/2.0-1.0/2.0*cos(CA[i])*sin(CA[i])))/(pow(X[i]*f,1.0/3.0)*pow(sin(CA[i]),2.0/3.0));

				A[i]=(log(X[i]+(1-X[i])*pow(Re/Reb,2.0))-log(X[i]))/9.0;//exponential coefficient describing wetting ratio increasing for previous tubes

//				fout<<"alpha for angle,"<<k*X[i]/h[i]<<",";
				aave += k*X[i]/b[i]/N;/////calculation of average of alpha

		}
		
	}
	
	else{
		h[i]=pow(3.0*m/L*myu/(rho*rho*gl*sin(beta*i)),1.0/3.0);
	//second half of the tube
	if(h[N/2-1]>=hmin){
		X[i]=1.0;
//			fout<<"alpha for angle,"<<k*X[i]/h[i]<<",";
			aave += k*X[i]/h[i]/N;/////calculation of average of alpha

		}
		else{

		
//		fout<<h[i]<<"h,";
		
		
			hminp=pow(((rho*rho*rho)*gl*gl*sin(beta*i)*sin(beta*i))/(15.0*myu*myu*st),0.2)*hmin;
//			cout<<hminp<<"hminp"<<endl;

			CA[i]= 1665.4*pow(hminp,6.0) - 3446.8*pow(hminp,5.0) + 2399.3*pow(hminp,4.0) - 459.73*pow(hminp,3.0) + 49.538*pow(hminp,2.0) + 2.9741*hminp + 0.8616;//////////////////////////////////////////////////////
			CA[i]=CA[i]/180*PI;
				f=-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])
		-3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i];//f(θo)
	g=(CA[i]*(5.0/16.0+15.0/4.0*pow(cos(CA[i]),2.0) +5.0/2.0*pow(cos(CA[i]),4.0))
	  -sin(CA[i])*(113.0/48.0*cos(CA[i]) +97.0/24.0*pow(cos(CA[i]),3.0) +1.0/6.0*pow(cos(CA[i]),5.0)));//Ψ(θo))

				X[i]=pow((h[i]),3.0)*sin(CA[i])/f*pow(2.0/45.0*(pow(rho,3.0)
					*pow((9.81*sin(beta*i)),2.0)/(st*pow(myu,2.0)))*g/sin(CA[i])*pow((CA[i]/sin(CA[i])-cos(CA[i])),-1),3.0/5.0);
				if(X[i]>=1){
					X[i]=1.0;
					aave += k*X[i]/h[i]/N;
				}
				else{
				b[i]=(h[i]*(CA[i]/2.0-1.0/2.0*cos(CA[i])*sin(CA[i])))/(pow(X[i]*f,1.0/3.0)*pow(sin(CA[i]),2.0/3.0));

				A[i]=(log(X[i]+(1-X[i])*pow(Re/Reb,2.0))-log(X[i]))/9.0;//exponential coefficient describing wetting ratio increasing for previous tubes

//				fout<<"alpha for angle,"<<k*X[i]/h[i]<<",";
				aave += k*X[i]/b[i]/N;/////calculation of average of alpha
				}
	}
		
	}
	}

	for(i=N-3;i<N;i++){
		
				X[i]=1.0;
					
//				fout<<"alpha for angle,"<<k*X[i]/h[i]<<",";
				aave += k*X[i]/b[i]/N;/////calculation of average of alpha

		}

	aavet+=aave/10;

	for(j=0;j<9;j++){
	aave=0;
		for(i=1;i<N;i++){
	
			if(X[i]==1){
				a[i]=k/h[i];
				aave+=a[i]/N;
			}
			else{
				b[i]=(h[i]*(CA[i]/2.0-1.0/2.0*cos(CA[i])*sin(CA[i])))/(pow((X[i]+(1-X[i])*pow(Re/Reb,2.0))*exp(-A[i]*(j))*(-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])
		-3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i]),1.0/3.0)*pow(sin(CA[i]),2.0/3.0));

			    a[i]=k/b[i]*(X[i]+(1-X[i])*pow(Re/Reb,2.0))*exp(-A[i]*(j));//local HTC calculated with exponential law-increasing wetting ratio
				aave+=a[i]/N;
			}

		}
		
		aavet+=aave/10;
	}
	fout<<aavet<<endl;

	// fout<<endl;

	}




}*/
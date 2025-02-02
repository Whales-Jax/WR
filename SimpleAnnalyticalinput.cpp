/*///熱交換器
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "CNewtonRaphsonMethod.h"
#include"CFluidParameter.h"
#include"PropertyLithiumBromide_ver05.h"
#include"PropertyWater_ver05.h"

PropertyLithiumBromide libr;
PropertyWater wat;

int main()
{
	int p,i,j,N,M,o,O,E,n,r,s,S,R,d;
	N=80;/////Angle's mesh
	M=80;/////Thickness's mesh
	E=6;//////Eigenvalues
	O=1;/////Mass flow rate
	R=1;//////Parameter
	S=220;/////Eigenfunction terms

	double DF0,Ti,Tsat,Tw,Te,Teti,Xi,Xe,P,st,rho,myu,ni,g,k,cp,D,a,L,ro,PI,beta,de,dh,hv,Le,Aa,Lc,Sc,Pr,da,m,Re,Ga,Reb,WR,lnc,SDF,SF,SDG,SG,Fnc,FN,GN,Fj,Gj,FJ,GJ,Ie,Dhtc,Dhtcu,Save,Suave,Gave,mtcave,mtcuave,htcave,htcuave,Nuave,Nuuave,Qave,Quave,qave,quave,teta;
	double Lamb[100],h[100],Tavu[100],Xavu[100],Tav[100],Xav[100],q[100],qu[100],Sh[100],Shu[100],Gv[100],mtc[100],mtcu[100],ln[100],Nu[100],Nuu[100],htc[100],htcu[100],ani[250],bni[250],A[100],B[100],NB[100],DB[100],DF[100],DG[100],IF[100];
	double u[100][100],X[100][100],Xu[100][100],T[100][100],Tu[100][100],an[100][250],bn[100][250],F[100][100],G[100][100];

	ofstream fout("AnalyticHMTw.csv");
	if(!fout){
	cout << "ファイルをオープンできませんでした"<< endl;
    return 1;
	}
	else
		cout<< "ファイルをオープンしました" << endl;

//////////////////////////////////////////////////////////////////////////////////////Parameter loop/////////////////

	for(r=1;r<=R;r++){

/////////////////////////////////////////////////////////////////////////////operative conditions////////////////////

	Ti=50.0;
	fout<<Ti<<",";
	Xi=60.0/100.0;
	fout<<Xi<<",";
	Tw=30.0;
	fout<<Tw<<",";

	
	P=1.5;//[kPa]
	fout<<P<<",";
	Tsat=wat.p_t( P );

	st=libr.sc_st_XT(Xi,Ti);
	rho=libr.sc_rho_XT(Xi,Ti);
	myu=libr.sc_visc_XT(Xi,Ti);	
	ni=myu/rho;
	g=9.81;//[m/s2]
	k=libr.sc_thc_XT(Xi,Ti);
	cp=libr.sc_cp_XT(Xi,Ti)*1000.0;////[J/kgK]
	D=libr.sc_d_XT(Xi,Ti);//////[m2/s]
    a=k/rho/cp;
	L=0.878;//tube length[m]
	ro=0.008;
	fout<<ro<<",";
	PI=(6*asin(0.5));//π
	teta=PI*50.0/180.0;
	beta=PI/double(N);//Angle of tube//////for 90degree
	de=beta/double(100);
	dh=1.0/double(100*M);
	Te=libr.sc_T_XTsat(Xi,Tsat);
	fout<<Te<<",";
	Xe=libr.sc_X_TTsat(Ti,Tsat);
	fout<<Xe<<",";
	Teti=(Ti-Tw)/(Te-Tw);
	hv = wat.sat_hv( Tsat )*1000.0;//-libr.sc_h_XT(Xi[i][M],Ti[i][M])*1000.0;/////[J/kg]?????check calculation formula

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Le=a/D;
	fout<<Le<<",";
	Aa=-hv*(Xe-Xi)/(cp*(Te-Ti)*Xe); ////check!!!!
	fout<<Aa<<",";
	Lc=pow(ni*ni/g,1.0/3.0);
	Sc=myu/(rho*D);
	fout<<Sc<<",";
	Pr=myu*cp/k;
	fout<<Pr<<",";
	da=2.0*PI*ro/Lc;
	fout<<da<<",";
	Ga=(rho*st*st*st)/(myu*myu*myu*myu*g);
//	Reb=34.1*pow(Ga,0.051);
	Reb=7.03*pow(Ga,0.2)*pow(0.2*log(teta)+0.69,3.0);

//////////////////////////////////////////////////////////////////////////////////////massflowrate loop/////////////////
	
	for(o=1;o<=O;o++){
	m=0.20712;//*o;//[kg/s]
	fout<<m/L<<",";
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	Re=4.0*m/L/myu;
	fout<<Re<<",";	
	
	if (Re<Reb) {
		WR=Re/Reb+0.1*(o-1);
		fout<<WR<<",";
	}

	else {
		WR=1;
		fout<<WR<<",";
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(j=0;j<=M;j++){

		T[0][j]=Ti;
		Tu[0][j]=Ti;
		X[0][j]=Xi;
		Xu[0][j]=Xi;
		u[0][j]=0;
		

	}

	for(j=0;j<=M;j++){
	for(i=1;i<=N;i++){
		T[i][j]=Tw;
		Tu[i][j]=Tw;
		X[i][j]=Xe;
		Xu[i][j]=Xe;
		Nu[i]=0;
		Nuu[i]=0;
		Sh[i]=0;
		Shu[i]=0;
		htc[i]=0;
		htcu[i]=0;
		Gv[i]=0;
		mtc[i]=0;
		mtcu[i]=0;
		q[i]=0;
		qu[i]=0;
		Tavu[i]=0;
		Xavu[i]=0;
		Tav[i]=0;
		Xav[i]=0;
		h[i]=pow(3.0*myu*m/L/(2.0*WR*pow(rho,2.0)*g*sin(beta*i)),1.0/3.0);
		Lamb[i]=0;

						}
						}

//////////////////////////////////////////////////////////////////////eigenfunctions initialization////////////////////

	

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////eigenvalues//////////////////////////////////

		ln[1]=0.06017;
		ln[2]=0.167391;
		ln[3]=0.2759875;
		ln[4]=0.3855765;
		ln[5]=0.4936935;
		ln[6]=0.599232;
		//ln[7]=0.702878;
		ln[7]=4.3799465;
		ln[8]=5.266974;
		ln[9]=6.1524248;
		//ln[10]=1.43276;
		//ln[10]=2.26341;
		//ln[11]=3.05577;
		//ln[12]=3.91461;
		//ln[13]=4.7213;
		//ln[14]=5.5314;
		//ln[15]=6.29;

	for(n=1;n<=E;n++){

		an[n][0]=0;
		bn[n][0]=1;
		an[n][1]=1;
		bn[n][1]=0;
		an[n][2]=0;
		bn[n][2]=0;
		an[n][3]=0;
		bn[n][3]=-ln[n]*ln[n]*Sc/3.0;

		for(s=4;s<=S;s++){

			an[n][s]=Pr*ln[n]*ln[n]*(an[n][s-4]-2.0*an[n][s-3])/double((s*(s-1)));
			bn[n][s]=Sc*ln[n]*ln[n]*(bn[n][s-4]-2.0*bn[n][s-3])/double((s*(s-1)));

		}

	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////Coefficients_An/Bn////////////////

	for(n=1;n<=E;n++){

		DF0=an[n][1];
		GN=0;
		FN=0;
		DF[n]=0;
		IF[n]=0;
		DG[n]=0;
		
		F[n][0]=0;
		G[n][0]=bn[n][0];
		NB[n]=0;
		DB[n]=0;

		for(s=0;s<=S;s++){

			GN +=bn[n][s];
			FN +=an[n][s];
			DF[n] +=s*an[n][s];
			DG[n] +=s*bn[n][s];
			IF[n] +=an[n][s]/double(s+1);
			

		}

			Gj=G[n][0];
			Fj=0;

		for(j=1;j<=100*M;j++){


			FJ=0;
			GJ=0;

			for(s=0;s<=S;s++){
			
			Fj +=an[n][s]*pow(dh*j,double(s));
			Gj +=bn[n][s]*pow(dh*j,double(s));

			}
				
				NB[n] +=dh/2.0*((2.0*dh*j-pow(dh*j,2.0))*(Le*Pr*Teti*GN/FN*FJ+Aa*Sc*GJ)+(2.0*dh*(j-1)-pow(dh*(j-1),2.0))*(Le*Pr*Teti*GN/FN*Fj+Aa*Sc*Gj));
				DB[n] +=dh/2.0*((2.0*dh*j-pow(dh*j,2.0))*(Le*Pr*pow(GN,2.0)/pow(FN,2.0)*pow(FJ,2.0)+Aa*Sc*pow(GJ,2.0))+(2.0*dh*(j-1)-pow(dh*(j-1),2.0))*(Le*Pr*pow(GN,2.0)/pow(FN,2.0)*pow(Fj,2.0)+Aa*Sc*pow(Gj,2.0)));
				
				Gj=GJ;
				Fj=FJ;
		
		}

		B[n]=NB[n]/DB[n];
		A[n]=B[n]*GN/FN;

		for(j=1;j<=M;j++){

			F[n][j]=0;
			G[n][j]=0;
			

			for(s=0;s<=S;s++){
			
			F[n][j] +=an[n][s]*pow(1.0/double(M)*double(j),double(s));
			G[n][j] +=bn[n][s]*pow(1.0/double(M)*double(j),double(s));

			}
		}

	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////angle loop//////////////////////



	Ie=0;
	htcave=0;
	htcuave=0;
	Nuave=0;
	Nuuave=0;
	Save=0;
	Suave=0;
	Gave=0;
	mtcave=0;
	mtcuave=0;
	Qave=0;
	Quave=0;
	qave=0;
	quave=0;
	
	for(i=1;i<N;i++){	

/////////////////////////////////////////////////////////////temperature and concentration distributions//////////////

		for(d=1;d<=100;d++){

		Ie +=de/2.0*(pow(sin(beta*(i-1)+de*(d-1)),1.0/3.0)+pow(sin(beta*(i-1)+de*(d)),1.0/3.0));

		}
		for(j=1;j<=M;j++){
	u[i][j]=(rho*g*h[i]*h[i]/2.0/myu*sin(beta*i)*(2.0/(double(M))*double(j)-pow(double(j)/(double(M)),2.0)));
			Lamb[i] +=rho*u[i][j]*h[i]/(double(M))*WR;
	}
		for(j=1;j<=M;j++){
		
			for(n=1;n<=E;n++){	

			T[i][j] +=-(Te-Tw)*A[n]*F[n][j]*exp(-ln[n]*ln[n]*da*pow(2*8.0*WR/(3.0*Re),4.0/3.0)*Ie);
			Tu[i][j] +=-(Te-Tw)*A[n]*F[n][j]*exp(-ln[n]*ln[n]*da*pow(2*8.0/(3.0*Re),4.0/3.0)*Ie);
			X[i][j] +=(Xi-Xe)*B[n]*G[n][j]*exp(-ln[n]*ln[n]*da*pow(2*8.0*WR/(3.0*Re),4.0/3.0)*Ie);
			Xu[i][j] +=(Xi-Xe)*B[n]*G[n][j]*exp(-ln[n]*ln[n]*da*pow(2*8.0/(3.0*Re),4.0/3.0)*Ie);

			}
			
			Tav[i] +=T[i][j]*rho*h[i]/(double(M))*WR/Lamb[i]*u[i][j];
			Xav[i] +=X[i][j]*rho*h[i]/(double(M))*WR/Lamb[i]*u[i][j];
			Tavu[i] +=Tu[i][j]*rho*h[i]/(double(M))*WR/Lamb[i]*u[i][j];
			Xavu[i] +=Xu[i][j]*rho*h[i]/(double(M))*WR/Lamb[i]*u[i][j];
			
				
			
		}

//////////////////////////////////////////////////////////Heat and Mass transfer///////////////////////////////////////

			Dhtc=0;
			Dhtcu=0;

			for(n=1;n<=E;n++){

			Dhtc +=-(Te-Tw)*A[n]*IF[n]*exp(-ln[n]*ln[n]*da*pow(2.0*8.0*WR/(3.0*Re),4.0/3.0)*Ie);
			Dhtcu +=-(Te-Tw)*A[n]*IF[n]*exp(-ln[n]*ln[n]*da*pow(2.0*8.0/(3.0*Re),4.0/3.0)*Ie);

			}

		for(n=1;n<=E;n++){

		htc[i] +=-(Te-Tw)*pow(WR,4.0/3.0)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*k/Lc*A[n]*DF0*exp(-ln[n]*ln[n]*da*pow(2*8.0*WR/(3.0*Re),4.0/3.0)*Ie)/Dhtc;
		htcu[i] +=-(Te-Tw)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*k/Lc*A[n]*DF0*exp(-ln[n]*ln[n]*da*pow(2*8.0/(3.0*Re),4.0/3.0)*Ie)/Dhtcu;
		Nu[i] +=-(Te-Tw)*pow(WR,4.0/3.0)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*A[n]*DF0*exp(-ln[n]*ln[n]*da*pow(2*8.0*WR/(3.0*Re),4.0/3.0)*Ie)/Dhtc;
		Nuu[i] +=-(Te-Tw)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*A[n]*DF0*exp(-ln[n]*ln[n]*da*pow(2*8.0/(3.0*Re),4.0/3.0)*Ie)/Dhtcu;
		mtc[i] +=(Xe-Xi)*pow(WR,4.0/3.0)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*D/Lc*B[n]*DG[n]*exp(-ln[n]*ln[n]*da*pow(2*8.0*WR/(3.0*Re),4.0/3.0)*Ie)/(X[i][M]*(X[i][1]-X[i][M]));
		mtcu[i] +=(Xe-Xi)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*D/Lc*B[n]*DG[n]*exp(-ln[n]*ln[n]*da*pow(2*8.0/(3.0*Re),4.0/3.0)*Ie)/(Xu[i][M]*(Xu[i][1]-Xu[i][M]));
		Sh[i] +=(Xe-Xi)*pow(WR,4.0/3.0)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*B[n]*DG[n]*exp(-ln[n]*ln[n]*da*pow(2*8.0*WR/(3.0*Re),4.0/3.0)*Ie)/(X[i][M]*(X[i][1]-X[i][M]));
		Shu[i] +=(Xe-Xi)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*B[n]*DG[n]*exp(-ln[n]*ln[n]*da*pow(2*8.0/(3.0*Re),4.0/3.0)*Ie)/(Xu[i][M]*(Xu[i][1]-Xu[i][M]));
		Gv[i] +=(Xe-Xi)*pow(WR,4.0/3.0)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*D*rho/Lc*B[n]*DG[n]*exp(-ln[n]*ln[n]*da*pow(2*8.0*WR/(3.0*Re),4.0/3.0)*Ie)/X[i][M];
		q[i] +=-(Ti-Tw)*pow(WR,4.0/3.0)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*k/Lc*A[n]*DF0*exp(-ln[n]*ln[n]*da*pow(2*8.0*WR/(3.0*Re),4.0/3.0)*Ie);
		qu[i] +=-(Ti-Tw)*pow(2*8.0*sin(PI*double(i)/double(N))/(3.0*Re),1.0/3.0)*k/Lc*A[n]*DF0*exp(-ln[n]*ln[n]*da*pow(2*8.0/(3.0*Re),4.0/3.0)*Ie);

		}

		htcave +=htc[i]/(double(N-1));
		htcuave +=htcu[i]/(double(N-1));
		Nuave +=Nu[i]/(double(N-1));
		Nuuave +=Nuu[i]/(double(N-1));
		Save +=Sh[i]/(double(N-1));
		Gave +=Gv[i]/(double(N-1));
		Suave +=Shu[i]/(double(N-1));
		mtcave +=mtc[i]*3600/(double(N-1));
		mtcuave +=mtcu[i]*3600/(double(N-1));
		qave +=q[i]/(double(N-1));//per tube unit area
		quave +=qu[i]/(double(N-1));//per tube unit area
		Qave +=2.0*PI*ro*q[i]/(double(N-1));//per tube unit length
		Quave +=2.0*PI*ro*qu[i]/(double(N-1));//per tube unit length
	
	}

for(j=0;j<=M;j++){
	for(i=0;i<=N;i++){
		fout<<T[i][j]<<",";
						}
	fout<<0<<endl;
						}

for(j=0;j<=M;j++){
	for(i=0;i<=N;i++){
		fout<<X[i][j]<<",";
						}
	fout<<0<<endl;
						}
for(j=0;j<=M;j++){
	for(i=0;i<=N;i++){
		fout<<u[i][j]<<",";
						}
	fout<<0<<endl;
						}

	for(i=0;i<=N;i++){
		fout<<Lamb[i]<<",";
		fout<<h[i]<<",";
	fout<<Tav[i]<<",";
		fout<<Xav[i]<<",";
		fout<<Tavu[i]<<",";
		fout<<Xavu[i]<<",";
		fout<<htc[i]<<",";
		fout<<htcu[i]<<",";
		fout<<Nu[i]<<",";
		fout<<Nuu[i]<<",";
		fout<<Sh[i]<<",";
		fout<<Shu[i]<<",";
		fout<<Gv[i]<<",";
		fout<<mtc[i]*3600<<",";
		fout<<mtcu[i]*3600<<",";
		fout<<q[i]<<",";
		fout<<qu[i]<<endl;

						}
	fout<<Lamb[N-1]<<",";
	fout<<Tav[N-1]<<",";
		fout<<Xav[N-1]<<",";
		fout<<Tavu[N-1]<<",";
		fout<<Xavu[N-1]<<",";
		fout<<htcave<<",";
	fout<<htcuave<<",";
	fout<<Nuave<<",";
	fout<<Nuuave<<",";
	fout<<mtcave<<",";
	fout<<mtcuave<<",";
	fout<<Gave<<",";
	fout<<Suave<<",";
	fout<<Save<<",";
	fout<<qave<<",";
	fout<<quave<<",";
	fout<<Qave<<",";
	fout<<Quave<<endl;
						
}
/////////////////////////////////////////////////////////////////////////////Massflowrate loop/////////////////////////

	//for(n=1;n<=E;n++){	

	//	fout<<ln[n]<<",";
		//fout<<A[n]<<",";
		//fout<<B[n]<<endl;

		//}

}
/////////////////////////////////////////////////////////////////////////////Parameter loop/////////////////////////

}*/
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
	int p,i,j,N,M,Z;
	N=1;/////
	M=60;/////
	Z=60;

	double	A,B,D,Dh,d,st,slp,Sgla,Sgls,sgla,sgls,ug,ugf,ul,ulf,ugla,ugls,Ug,Ul,Up,Ugla,Ugls,tet,L,G,rhog,rhol,rhop,myug,myul;
	double  Mi,r,rg,rl,Rel,Reg,Re,g,CA,Ti,X,Xi,P,X0,f,E,W,WD,WDD,We,Weg,Wel,Alf,Alff,PI,m,eps;
//	double	uave[450],EM[450],EMM[100],ZZ[450],h[450],t[450],CA[450],X[450],XX[100],R[450],u[100][100],RR[100];

	ofstream fout("MiniChAlfaX.csv");
	if(!fout){
	cout << "ファイルをオープンできませんでした"<< endl;
    return 1;
	}
	else
		cout<< "ファイルをオープンしました" << endl;

	PI=(6*asin(0.5));//π

////////////////////Experiment condition////////////////////
	WD=0.1;
	Alf=0.6; ///Initial value
	Alff=0.97;

	m=0.015;
	fout<<m<<",";
	X=0.5;//+0.05*double(i-1);
	fout<<X<<",";
	Xi=0.0/100.0;
	fout<<Xi<<",";

	P=101.325;
	fout<<P<<",";
	Ti=wat.p_t( P );
	fout<<Ti<<",";

	L=0.3;//Channel length[m]
	fout<<L<<",";
	D=0.005;//Channel diameter[m]
	fout<<D<<",";

	G=m*4.0/(PI*D*D);
	
////////////////////Fluid Properties//////////////////////////////////

	rhog=wat.sat_rouv( Ti );//
	myug=wat.sat_myuv( Ti );;	
	st=libr.sc_st_XT(Xi,Ti);
	rhol=libr.sc_rho_XT(Xi,Ti);
	myul=libr.sc_visc_XT(Xi,Ti);	
	g=9.81;//[m/s2]

//	k=libr.sc_thc_XT(Xi,Ti);
//	d=libr.sc_d_XT(Xi,Ti);//////[m2/s]
	
	Ul=G*(1.0-X)/rhol/(1.0-Alf);
	Ug=G*X/rhog/Alf;
	fout<<Ug<<",";
	fout<<Ul<<",";
	Mi=myul/myug;
	//Mi=Mi/double(i+1);
	fout<<Mi<<",";
	r=rhol/rhog;
	//r=r/double(i+1);
	fout<<r<<",";
	slp=Ug/Ul;
	slp=r*X*(1.0-Alf)/Alf/(1.0-X);
	fout<<slp<<",";
	Ugla=G/rhol*(Mi*pow(Alf,3.0/2.0)*(1.0-X)+r*(1.0-Alf)*(1.0-pow(Alf,1.0/2.0))*X)/(Alf*(1.0-Alf)*(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0))));
	fout<<Ugla<<",";
	Ugls=G/rhol*((Mi*(1.0-X)/(1.0-Alf)-r*X/Alf)*pow(1.0-(1.0-Alf)/(WD),0.5)+r*X/Alf)/((Mi-1.0)*pow(1.0-(1.0-Alf)/(WD),0.5)+1.0);
	fout<<Ugls<<",";
	
	
//	rg=X-1.0/(slp*(1.0-X)*Alf/X/(1.0-Alf))*(1.0-X);
	rg=X-1.0/r*(1.0-X);
//	rl=X*(slp*(1.0-X)*Alf/X/(1.0-Alf)+1.0)-1.0;
	rl=X*(r+1.0)-1.0;
	fout<<rg<<",";
	fout<<rl<<",";
	Up=G/rhol*pow(pow(r,2.0)*X*X*X/Alf+pow(1.0-X,3.0)/(1.0-Alf),0.5);
	fout<<Up<<",";
	rhop=1.0/(X/rhog+(1.0-X)/rhol);
	fout<<rhop<<",";
	ug=Ug/Up;
	ug=r*pow(pow(X,2.0)/(Alf)*(1.0-Alf)/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5);
	ul=Ul/Up;
	ul=pow(pow(1.0-X,2.0)/(1.0-Alf)*Alf/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5);
	fout<<ug<<",";
	fout<<ul<<",";
	
//	ugla=pow((Mi*Mi*pow(Alf,3.0)*pow(slp*(1.0-X)*Alf/X/(1.0-Alf)*(1.0-Alf)*(1.0-pow(Alf,0.5))*X,2.0))/((pow(slp*(1.0-X)*Alf/X/(1.0-Alf)*(1.0-Alf)*X,2.0)*X+pow(Alf*(1.0-X),3.0)/Alf)*(Mi*Mi*Alf+pow((1.0-pow(Alf,0.5)),2.0)+2.0*Mi*pow(Alf,0.5)*(1.0-pow(Alf,0.5)))),0.5);
	ugla=Ugla/Up;
	ugla=pow((pow(Mi*pow(Alf,3.0/2.0)*(1.0-X)+r*(1.0-Alf)*(1.0-pow(Alf,1.0/2.0))*X,2.0))/(Alf*(1.0-Alf)*pow(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)),2.0)*(r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))),0.5);
		
	ugls=Ugls/Up;
	ugls=pow(pow(((Mi*(1.0-X)/(1.0-Alf)-r*X/Alf)*pow(1.0-(1.0-Alf)/(WD),0.5)+r*X/Alf),2.0)/pow(((Mi-1.0)*pow(1.0-(1.0-Alf)/(WD),0.5)+1.0),2.0)*Alf*(1.0-Alf)/((r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))),0.5);
		
	fout<<ugla<<",";
	fout<<ugls<<",";

	Sgla=PI*D*pow(Alf,0.5);
	fout<<Sgla<<",";
	Sgls=D*(1.0+(PI*WD-1.0)*pow(1.0-(1.0-Alf)/WD,0.5));
	fout<<Sgla<<",";
	sgla=pow(Alf,0.5);
	sgls=((1.0+(PI*WD-1.0)*pow(1.0-(1.0-Alf)/WD,0.5)))/PI;
	fout<<sgla<<",";
	fout<<sgls<<",";

	We=rhop*D*Up*Up/st;
	We=0.5;
//	We=0.0000005*pow(5.0,double(i-1));
	fout<<We<<",";
	Wel=pow(G*(1.0-X),2.0)*D/rhol/st;
	fout<<Wel<<",";
	Weg=pow(G*X,2.0)*D/rhol/st;
	fout<<Weg<<",";
	Rel=(1.0-X)*m*4/PI/D/myul;
	fout<<Rel<<",";
	Reg=X*m*4/PI/D/myug;
	fout<<Reg<<endl;
	
for(i=1;i<=N;i++){///////////////
	//if (We<=0.1){
	//	We=We-0.0025;
	//fout<<We<<",";
	//}
	//else{
	//	We=0.9900005-0.1*double(i-1);
	//fout<<We<<",";
	//}
//			for(p=1;p<=Z-1;p++){

//				Alf=0.7+double(3*p)/double(10*Z);
//				fout<<Alf<<",";

//								}
//			fout<<0<<endl;

//			for(p=1;p<=Z-1;p++){

//				Alf=0.7+double(3*p)/double(10*Z);
//				eps=0.5*((X-1.0/r*(1.0-X))*Alf*pow(r*pow(pow(X,2.0)/(Alf)*(1.0-Alf)/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5),3.0)+(X*(r+1.0)-1.0)*(1.0-Alf)*pow(pow(pow(1.0-X,2.0)/(1.0-Alf)*Alf/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5),3.0))+4.0/We*(pow(Alf,0.5))*(pow((pow(Mi*pow(Alf,3.0/2.0)*(1.0-X)+r*(1.0-Alf)*(1.0-pow(Alf,1.0/2.0))*X,2.0))/(Alf*(1.0-Alf)*pow(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)),2.0)*(r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))),0.5));
					
//				fout<<eps<<",";

//								}
//			fout<<0<<endl;
	
			
			
			for(p=1;p<=Z-1;p++){

				WDD=double(p)/double(2*Z);
				fout<<WDD<<",";

								}
			fout<<0<<endl;

			for(j=1;j<=M-1;j++){
				Alf=0.7+double(3*j)/double(10*M);
				fout<<Alf<<",";
			for(p=1;p<=Z-1;p++){
				WDD=double(p)/double(2*Z);
				eps=0.5*((X-1.0/r*(1.0-X))*Alf*pow(r*pow(pow(X,2.0)/(Alf)*(1.0-Alf)/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5),3.0)+(X*(r+1.0)-1.0)*(1.0-Alf)*pow(pow(pow(1.0-X,2.0)/(1.0-Alf)*Alf/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5),3.0))+4.0/We*(((1.0+(PI*WDD-1.0)*pow(1.0-(1.0-Alf)/WDD,0.5)))/PI)*(pow(pow(((Mi*(1.0-X)/(1.0-Alf)-r*X/Alf)*pow(1.0-(1.0-Alf)/(WDD),0.5)+r*X/Alf),2.0)/pow(((Mi-1.0)*pow(1.0-(1.0-Alf)/(WDD),0.5)+1.0),2.0)*Alf*(1.0-Alf)/((r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))),0.5));
				fout<<eps<<",";

								}
			fout<<0<<endl;

								}


//			WD=0.5;
//	Alf=0.65; ///Initial value


				CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
	mnm.setup(2000,1+1e-12,1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

		p=0;

//			mnm.setValue( p , Alf);
//		p++;

			mnm.setValue( p , Alff);
		p++;

//			mnm.setValue( p , WD);
//		p++;

//			mnm.setValue( p , X);
//		p++;

		mnm.setAcc(0.01);				//加速度勾配の入力　※なくても可

	mnm.initial();					//計算を開始する前に必ず初期化してください
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// おまじない
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// おまじない

		p=0;
			//値の設定
//			Alf = mnm.getValue(p);
//		p++;

			Alff = mnm.getValue(p);
		p++;

//			WD = mnm.getValue(p);
//		p++;

//			X = mnm.getValue(p);
//		p++;

		p=0;

//			mnm.setError( p , 0.0, 0.5*((pow(r,3.0)*X-r*r*(1.0-X))*(pow((X*X*(1.0-Alf))/(Alf*(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0))),3.0/2.0)+3.0/2.0*Alf*(pow((X*X*(1.0-Alf))/(Alf*(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0))),1.0/2.0))*((-X*X*Alf*((r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)))-X*X*(1.0-Alf)*(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)+Alf*(-r*r*X*X*X+pow(1.0-X,3.0))))/pow(Alf*((r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0))),2.0)))+(X*(r+1.0)-1.0)*((pow((pow(1.0-X,2.0)*Alf)/((1.0-Alf)*(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0))),3.0/2.0))+(1.0-Alf)*3.0/2.0*(pow((pow(1.0-X,2.0)*Alf)/((1.0-Alf)*(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0))),1.0/2.0))*((pow(1.0-X,2.0)*((1.0-Alf)*(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)))-pow(1.0-X,2.0)*Alf*(-(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0))+(1.0-Alf)*(-r*r*pow(X,3.0)+pow(1.0-X,3.0))))/pow((1.0-Alf)*(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),2.0))))
//				+2.0/We*1.0/pow((pow(Mi*pow(Alf,3.0/2.0)*(1.0-X)+r*(1.0-Alf)*(1.0-pow(Alf,1.0/2.0))*X,2.0))/((1.0-Alf)*pow(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)),2.0)*(r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))),0.5)*((2.0*(Mi*pow(Alf,3.0/2.0)*(1.0-X)+r*(1.0-Alf)*(1.0-pow(Alf,1.0/2.0))*X)*(Mi*(1.0-X)*3.0/2.0*pow(Alf,0.5)+r*X*(-1.0+3.0/2.0*pow(Alf,0.5)-0.5*pow(Alf,-0.5)))*(((1.0-Alf)*pow(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)),2.0)*(r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))))-((pow(Mi*pow(Alf,3.0/2.0)*(1.0-X)+r*(1.0-Alf)*(1.0-pow(Alf,1.0/2.0))*X,2.0)))*(2.0*(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)))*(0.5*Mi*pow(Alf,-1.0/2.0)-0.5*pow(Alf,-1.0/2.0))*(1.0-Alf)*((r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0)))+(pow(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)),2.0))*(-X*X*X*r*r*2.0*(1.0-Alf)+pow(1.0-X,3.0)-2.0*Alf*pow(1.0-X,3.0))))/pow(((1.0-Alf)*pow(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)),2.0)*(r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))),2.0)));
//		p++;

//		N  (pow(Mi*pow(Alf,3.0/2.0)*(1.0-X)+r*(1.0-Alf)*(1.0-pow(Alf,1.0/2.0))*X,2.0))
//		D1 pow(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)),2.0)
//		D2 (r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))

		mnm.setError( p , 0.0, 0.5*((pow(r,3.0)*X-r*r*(1.0-X))*(pow((X*X*(1.0-Alff))/(Alff*(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0))),3.0/2.0)+3.0/2.0*Alff*(pow((X*X*(1.0-Alff))/(Alff*(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0))),1.0/2.0))*((-X*X*Alff*((r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0)))-X*X*(1.0-Alff)*(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0)+Alff*(-r*r*X*X*X+pow(1.0-X,3.0))))/pow(Alff*((r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0))),2.0)))+(X*(r+1.0)-1.0)*((pow((pow(1.0-X,2.0)*Alff)/((1.0-Alff)*(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0))),3.0/2.0))+(1.0-Alff)*3.0/2.0*(pow((pow(1.0-X,2.0)*Alff)/((1.0-Alff)*(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0))),1.0/2.0))*((pow(1.0-X,2.0)*((1.0-Alff)*(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0)))-pow(1.0-X,2.0)*Alff*(-(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0))+(1.0-Alff)*(-r*r*pow(X,3.0)+pow(1.0-X,3.0))))/pow((1.0-Alff)*(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0)),2.0))))
				+4.0/We*(((PI*WD-1.0)/(2.0*PI*WD))*1.0/pow(1.0-(1.0-Alff)/WD,0.5)*((pow(pow(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff),2.0)/pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0)*Alff*(1.0-Alff)/((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0))),0.5)))
				+0.5*((((1.0+(PI*WD-1.0)*pow(1.0-(1.0-Alff)/WD,0.5)))/PI))*1.0/((pow(pow(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff),2.0)/pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0)*Alff*(1.0-Alff)/((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0))),0.5)))
				*(((((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0)))*pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0))*((2.0*(Alff-Alff*Alff)*(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff))*((Mi*(1.0-X)/pow(1.0-Alff,2.0)+r*X/Alff/Alff)*pow(1.0-(1.0-Alff)/WD,0.5)+(Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*0.5/WD/pow(1.0-(1.0-Alff)/WD,0.5)-r*X/Alff/Alff))
				+(1.0-2.0*Alff)*(pow(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff),2.0)))
				-(pow(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff),2.0)*Alff*(1.0-Alff))
				*(2.0*(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0))*((Mi-1.0)/WD*0.5*pow(1.0-(1.0-Alff)/(WD),-0.5))*(((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0))))+(pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0))*(-r*r*X*X*X+pow(1.0-X,3.0))))/pow(((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0)))*pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0),2.0))));
		p++;
		/////change to ALFFFFFFFFFFFFFFFFFFFF

//		mnm.setError( p , 0.0, +4.0/We*((pow(1.0-(1.0-Alff)/WD,0.5)+(PI*WD-1.0)/(2.0*PI*WD*WD)*(1.0-Alff)*pow(1.0-(1.0-Alff)/WD,-0.5))*(pow(pow(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff),2.0)/pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0)*Alff*(1.0-Alff)/((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0))),0.5))
//				+0.5*1.0/((pow(pow(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff),2.0)/pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0)*Alff*(1.0-Alff)/((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0))),0.5)))
//				*((((1.0+(PI*WD-1.0)*pow(1.0-(1.0-Alff)/WD,0.5)))/PI))*((2.0*(((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0)))*pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0))*Alff*(1.0-Alff)*(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff))*(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff))/WD/WD*0.5*(1.0-Alff)/pow(1.0-(1.0-Alff)/WD,0.5))-2.0*(pow(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff),2.0)*Alff*(1.0-Alff))*(((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0))))*(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0))*(Mi-1.0)/WD/WD*0.5*(1.0-Alff)/pow(1.0-(1.0-Alff)/WD,0.5))
//				/pow(((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0)))*pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0),2.0))));
//		p++;

//		mnm.setError( p , 0.5*((X-1.0/r*(1.0-X))*Alf*pow(r*pow(pow(X,2.0)/(Alf)*(1.0-Alf)/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5),3.0)+(X*(r+1.0)-1.0)*(1.0-Alf)*pow(pow(pow(1.0-X,2.0)/(1.0-Alf)*Alf/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5),3.0))+4.0/We*(pow(Alf,0.5))*(pow((pow(Mi*pow(Alf,3.0/2.0)*(1.0-X)+r*(1.0-Alf)*(1.0-pow(Alf,1.0/2.0))*X,2.0))/(Alf*(1.0-Alf)*pow(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)),2.0)*(r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))),0.5)),
//				0.5*((X-1.0/r*(1.0-X))*Alff*pow(r*pow(pow(X,2.0)/(Alff)*(1.0-Alff)/(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0)),0.5),3.0)+(X*(r+1.0)-1.0)*(1.0-Alff)*pow(pow(pow(1.0-X,2.0)/(1.0-Alff)*Alff/(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0)),0.5),3.0))
//				+4.0/We*(((1.0+(PI*WD-1.0)*pow(1.0-(1.0-Alff)/WD,0.5)))/PI)*(pow(pow(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff),2.0)/pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0)*Alff*(1.0-Alff)/((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0))),0.5)));
//		p++;

//		R pow(1.0-(1.0-Alf)/WD,0.5)
//		F1 (((1.0+(PI*WD-1.0)*pow(1.0-(1.0-Alf)/WD,0.5)))/PI)
//			F2 (pow(pow(((Mi*(1.0-X)/(1.0-Alf)-r*X/Alf)*pow(1.0-(1.0-Alf)/(WD),0.5)+r*X/Alf),2.0)/pow(((Mi-1.0)*pow(1.0-(1.0-Alf)/(WD),0.5)+1.0),2.0)*Alf*(1.0-Alf)/((r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))),0.5))
//		N pow(((Mi*(1.0-X)/(1.0-Alf)-r*X/Alf)*pow(1.0-(1.0-Alf)/(WD),0.5)+r*X/Alf),2.0)*Alf*(1.0-Alf)
//			D1,2 pow(((Mi-1.0)*pow(1.0-(1.0-Alf)/(WD),0.5)+1.0),2.0) 
//			D2 ((r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0)))
//			D ((r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0)))*pow(((Mi-1.0)*pow(1.0-(1.0-Alf)/(WD),0.5)+1.0),2.0)
			// mnm.prt();				//エラー表示
		mnm.prt_sum();			//エラーの合計を表示

	}
	}

				fout<<Alf<<",";
				fout<<Alff<<",";
				fout<<WD<<",";
				fout<<X<<",";


				//ug=G*X/rhog/Alf;
				//fout<<ug<<",";
				//ul=G*(1.0-X)/rhol/(1.0-Alf);
				//fout<<ul<<",";
//				fout<<CA0<<",";
				
				
	Mi=myul/myug;
	fout<<Mi<<",";
	r=rhol/rhog;
	fout<<r<<",";
	rg=X-1.0/r*(1.0-X);
	rl=X*(r+1.0)-1.0;
	fout<<rg<<",";
	fout<<rl<<",";
	ug=r*pow(pow(X,2.0)/(Alf)*(1.0-Alf)/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5);
	ul=pow(pow(1.0-X,2.0)/(1.0-Alf)*Alf/(r*r*pow(X,3.0)*(1.0-Alf)+Alf*pow(1.0-X,3.0)),0.5);
	fout<<ug<<",";
	fout<<ul<<",";
	ug=r*pow(pow(X,2.0)/(Alff)*(1.0-Alff)/(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0)),0.5);
	ul=pow(pow(1.0-X,2.0)/(1.0-Alff)*Alff/(r*r*pow(X,3.0)*(1.0-Alff)+Alff*pow(1.0-X,3.0)),0.5);
	fout<<ug<<",";
	fout<<ul<<",";
	slp=r*X*(1.0-Alf)/Alf/(1.0-X);
	fout<<slp<<",";
	slp=r*X*(1.0-Alff)/Alff/(1.0-X);
	fout<<slp<<",";
	ugla=pow((pow(Mi*pow(Alf,3.0/2.0)*(1.0-X)+r*(1.0-Alf)*(1.0-pow(Alf,1.0/2.0))*X,2.0))/(Alf*(1.0-Alf)*pow(Mi*pow(Alf,1.0/2.0)+(1.0-pow(Alf,1.0/2.0)),2.0)*(r*r*(1.0-Alf)*X*X*X+Alf*pow(1.0-X,3.0))),0.5);
	ugls=pow(pow(((Mi*(1.0-X)/(1.0-Alff)-r*X/Alff)*pow(1.0-(1.0-Alff)/(WD),0.5)+r*X/Alff),2.0)/pow(((Mi-1.0)*pow(1.0-(1.0-Alff)/(WD),0.5)+1.0),2.0)*Alff*(1.0-Alff)/((r*r*(1.0-Alff)*X*X*X+Alff*pow(1.0-X,3.0))),0.5);
	fout<<ugla<<",";
	fout<<ugls<<",";
	sgla=pow(Alf,0.5);
	sgls=((1.0+(PI*WD-1.0)*pow(1.0-(1.0-Alff)/WD,0.5)))/PI;
	fout<<sgla<<",";
	fout<<sgls<<",";
	fout<<0<<endl;
	
			
//		Alff=0.05;///Initial value
//		Dh=4.0*(D*D/2.0*(1.0-pow(Alff,0.5))-pow(D/2.0*(1.0-pow(Alff,0.5)),2.0))/D;
//		Rel=G*(1.0-X)/(1.0-Alff)*4.0*(D*D/2.0*(1.0-pow(Alff,0.5))-pow(D/2.0*(1.0-pow(Alff,0.5)),2.0))/D/myul;
//			A=pow(2.457*log(1.0/(pow(7.0/(G*(1.0-X)/(1.0-Alff)*4.0*(D*D/2.0*(1.0-pow(Alff,0.5))-pow(D/2.0*(1.0-pow(Alff,0.5)),2.0))/D/rhol),0.9)+0.27*eps/D)),16.0);
//			B=pow(37530.0/(G*(1.0-X)/(1.0-Alff)*4.0*(D*D/2.0*(1.0-pow(Alff,0.5))-pow(D/2.0*(1.0-pow(Alff,0.5)),2.0))/D/myul),16.0);
//			f=2.0*pow(pow(8.0/(G*(1.0-X)/(1.0-Alff)*4.0*(D*D/2.0*(1.0-pow(Alff,0.5))-pow(D/2.0*(1.0-pow(Alff,0.5)),2.0))/D/myul),12.0)+pow(pow(2.457*log(1.0/(pow(7.0/(G*(1.0-X)/(1.0-Alff)*4.0*(D*D/2.0*(1.0-pow(Alff,0.5))-pow(D/2.0*(1.0-pow(Alff,0.5)),2.0))/D/myul),0.9)+0.27*eps/D)),16.0)+pow(37530.0/(G*(1.0-X)/(1.0-Alff)*4.0*(D*D/2.0*(1.0-pow(Alff,0.5))-pow(D/2.0*(1.0-pow(Alff,0.5)),2.0))/D/rhol),16.0),-1.5),1.0/12.0);
//		f=16.0*myul*G*(1.0-X)/D/rhol/pow(1.0-Alff,2.0);

//				CNewtonRaphsonMethod mnm2;		// CNewtonRaphsonMethodクラスの宣言
//	mnm2.setup(2000,1+1e-12,1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

//		p=0;

//			mnm2.setValue( p , Alff);
//		p++;
//			mnm2.setValue( p , X0);
//		p++;

//		mnm2.setAcc(0.03);				//加速度勾配の入力　※なくても可

//	mnm2.initial();					//計算を開始する前に必ず初期化してください
//	for(mnm2.main_loop_init();mnm2.main_loop_check();mnm2.main_loop_reinit()){		// おまじない
//		for(mnm2.sub_loop_init();mnm2.sub_loop_check();mnm2.sub_loop_reinit()){	// おまじない

//		p=0;
			//値の設定
//			Alff = mnm2.getValue(p);
//		p++;
//			X0 = mnm2.getValue(p);
//		p++;

//		p=0;

//			mnm2.setError( p , 0.0,  G*G*G*((1.0+64.0*L*myul*G*(1.0-X)/(D*D*pow(1.0-Alff,3.0)))*pow(1.0-X,3.0)/(rhol*rhol*pow(1.0-Alff,3.0))+pow(1.0-X,3.0)/(2.0*rhol*rhol*pow(1.0-Alff,2.0))*192.0*L*myul*G*(1.0-X)/(D*D*rhol*pow(1.0-Alff,4.0))-X*X*X/rhog/rhog/pow(Alff,3.0))-PI*D*st/2.0/pow(Alff,0.5));
//		p++;
//			mnm2.setError( p , pow(h0/R0,3.0), 2.0*pow(R0,1.0)*f/L0+2.0*(X0/2.0-R0*sin(CA0)/L0)*pow((1-cos(CA0)), 3.0) );
//		p++;

////////////////Continuity equation///////////////////////////////////////////////////////////////////////////////////////////////////////////


		// mnm2.prt();				//エラー表示
//		mnm2.prt_sum();			//エラーの合計を表示

//	}
//	}

//			fout<<Alff<<",";
//			ugf=G*X/rhog/Alff;
//			fout<<ugf<<",";
//			ulf=G*(1.0-X)/rhol/(1.0-Alff);
//			fout<<ulf<<",";


	

		

	}




}*/
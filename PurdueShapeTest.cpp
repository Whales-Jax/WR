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
PropertyWater wat2;

int main()
{
	int p,i,j,N,M,O,k,Z,E;
	N=2;	/////Stream-wise position mesh
	M=40;	/////Thickness's mesh
	Z=40;	/////Rivulet's transversal mesh
	O=1;	/////Mass flow rate
	E=30;	/////Energy step
	double	st,Dst,dR,dz,dX,dedX,gamma,gamma0,uav,uav0,uaveu,teta,Z0,L,G,L00,D,L0,rho,myu,beta,hmin,hminM,hminp,Re,h0,g,CA0,TI,TII,Ti,Xi,PI,X0,X0b,X0M,Xave,R0,R0b,R00,f,m,gl,n,ma,Reb,EM0,EMu,dp,P,h0p;
	double	uave[450],EM[450],EMM[100],ZZ[450],h[450],t[450],CA[450],X[450],XX[100],R[450],u[100][100],RR[100];

	ofstream fout("VWWRL.csv");
	if(!fout){
	cout << "ファイルをオープンできませんでした"<< endl;
    return 1;
	}
	else
		cout<< "ファイルをオープンしました" << endl;

	PI=(6*asin(0.5));//Pﾎ

////////////////////Experiment condition////////////////////

	//Wall inclination
	beta=PI/2.0;
	//Vapour/Gas pressure
	P=101.325;//kPa
		fout<<P<<",";
	//Inlet conditions
	Xi=00.0/100.0;
//		fout<<Xi<<",";
	ma=0.005;//[kg/s]
		fout<<ma<<",";
	Ti=25.0;
		fout<<Ti<<",";
	//Test section geometry
	L=0.2;	//surface width[m]
		fout<<L<<",";
	D=0.4;	//surface length[m]
		fout<<D<<",";

	
	for(n=1;n<=O;n++){///////////////Mass flow rate loop


////////////////////Fluid Properties/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Surface tension
	st=2.0*pow(10.0, -12.0)*pow(Ti, 4.0)-9.0*pow(10.0, -10.0)*pow(Ti, 3.0)-9.0*pow(10.0, -8.0)*pow(Ti, 2.0)-0.0002*Ti+0.0758;
	//Density
	rho=wat2.sc_roul(P,Ti);
	//Viscosity
	myu=wat2.sc_myul(P,Ti);	
	gl=9.81;//[m/s2]
	//Contact angle
	CA0=56.80*PI/180.0;
		fout<<CA0<<",";
	//Specific flowrate
	m=ma;
		fout<<m/L<<",";
	//Reynolds number
	Re=2.0*m/L/myu;
		fout<<Re<<endl;

//	k=libr.sc_thc_XT(Xi,Ti);
//	Reb=4.0*8.5027*pow(dp,0.0505);//Break-up Reynolds from Maron1982
//	fout<<Reb<<endl;


///////////////Steady wetting////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//			hmin=pow(3.0*Reb/4.0*pow(myu/rho,2.0)/9.81,1.0/3.0);//minimum thickness [mm] from Maron1982
//				fout<<hmin<<",";
			
	//Film thickness at the present flowrate
	h0=pow(3.0*m/L*myu/(rho*rho*gl*sin(beta)),1.0/3.0);
		fout<<h0<<",";

	//Minimum film rhickness from Mikielewicz (1976)
	hminp=0.1593*log(CA0)+0.6699;
	hminM=pow(((rho*rho*rho)*gl*gl*sin(beta)*sin(beta))/(15.0*myu*myu*st),-0.2)*hminp;
		fout<<hminM<<",";

	//Integral parameters/contact angle functions
	f=-1.0/4.0*pow(cos(CA0),3.0)*sin(CA0) -13.0/8.0*cos(CA0)*sin(CA0)-3.0/2.0*CA0*pow(sin(CA0),2.0) +15.0/8.0*CA0;//f(θo)
	g=(CA0*(5.0/16.0+15.0/4.0*pow(cos(CA0),2.0) +5.0/2.0*pow(cos(CA0),4.0))-sin(CA0)*(113.0/48.0*cos(CA0) +97.0/24.0*pow(cos(CA0),3.0) +1.0/6.0*pow(cos(CA0),5.0)));//Ψ(θo))
	
	//Wetting ability at the breaking point from Mikielewicz (1976)
	X0M=pow((hminM),3.0)*sin(CA0)/f*pow(2.0/45.0*(pow(rho,3.0)*pow((9.81*sin(beta)),2.0)/(st*pow(myu,2.0)))*g/sin(CA0)*pow((CA0/sin(CA0)-cos(CA0)),-1),3.0/5.0);
		fout<<X0M<<",";
	X0=pow((h0),3.0)*sin(CA0)/f*pow(2.0/45.0*(pow(rho,3.0)*pow((9.81*sin(beta)),2.0)/(st*pow(myu,2.0)))*g/sin(CA0)*pow((CA0/sin(CA0)-cos(CA0)),-1),3.0/5.0);
		fout<<X0<<",";

	//Rivulet radius at the breaking point from Mikielewicz (1976)
	R00=hminM*pow(sin(CA0)/X0M/f,1.0/3.0);
		fout<<R00<<",";

	//Rivulet spacing at the breaking point from Mikielewicz (1976)
	L00=2.0*R00*sin(CA0)/X0M;
		fout<<L00<<",";

	//Test section width [m] (reference rivulet spacing for this model)
	L0=L;
		fout<<L0<<",";

	//Average film velocity
	uaveu=0;

	for(j=1;j<=M;j++){

	u[0][j]=(rho*gl*h0*h0*sin(beta)/myu)*(double(j)/double(M)-1.0/2.0*pow(double(j)/double(M), 2.0));
		uaveu +=u[0][j]/double(M);

						}
		fout<<uaveu<<",";

	//Rivulet energy from Mikielewicz (1976)
	EM[0]=rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu/15.0*pow(X0M, -2.0/3.0)*pow(sin(CA0)/f,5.0/3.0)*pow(hminM, 5.0)*g/sin(CA0)+X0M*st*(CA0/sin(CA0)-cos(CA0))+st*cos(CA0);
		fout<<EM[0]<<",";
	EM[0]=rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu/15.0*pow(X0, -2.0/3.0)*pow(sin(CA0)/f,5.0/3.0)*pow(h0, 5.0)*g/sin(CA0)+X0*st*(CA0/sin(CA0)-cos(CA0))+st*cos(CA0);
		fout<<EM[0]<<",";
		
	//Film energy 
	EMu=1.0/15.0*rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(h0, 5.0)+st;
		fout<<EMu<<endl;
	
	//Forced wetting at the distributor inlet
	X[0]=1.0;
		fout<<X[0]<<",";
	X[1]=1.0;
		fout<<X[1]<<",";

	//Guessed radius			
	R0=R00/5;

	//Guessed hmin
	hmin=hminM;

	//Calculation of Rivuet Radius R at the breaking condition and minimum film thickness hmin from this model/////////////////////////////////////////
				CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
	mnm.setup(2000,1+1e-12,1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

		p=0;

			mnm.setValue( p , R0);
		p++;

			mnm.setValue( p , hmin );
		p++;

		mnm.setAcc(0.03);				//加速度勾配の入力　※なくても可

	mnm.initial();					//計算を開始する前に必ず初期化してください
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// おまじない
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// おまじない

		p=0;
	//値の設定
			R0 = mnm.getValue(p);
		p++;

			hmin = mnm.getValue(p);
		p++;

		p=0;

	//minimal rivulet energy
			mnm.setError( p , 0.0, rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*(+4.0/5.0*pow(R0, 5.0)/L0*g+2.0/15.0*pow(hmin, 3.0)*R0*pow(1.0-cos(CA0), 2.0)-4.0/5.0*pow(R0, 5.0)/L0*f*pow(1.0-cos(CA0), 2.0))+st*((-3.0*pow(hmin, 3.0)/(pow(R0, 4.0)*pow(1.0-cos(CA0), 3.0))-2.0*f/pow(1.0-cos(CA0), 3.0)/L0+2.0*sin(CA0)/L0)*(1.0-cos(CA0))+2.0/L0*(CA0-sin(CA0))) );
		p++;

	//energy balance
			mnm.setError( p , 1.0/15.0*rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(hmin, 5.0)+st, rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(R0, 5.0)*(+2.0/15.0*R0/L0*g+pow(1.0-cos(CA0), 2.0)/15.0*(pow(hmin/R0, 3.0)-2.0*R0*f/L0))+st*((pow(hmin/R0, 3.0)/pow(1.0-cos(CA0), 3.0)-2.0*R0*f/pow(1.0-cos(CA0), 3.0)/L0+2.0*R0*sin(CA0)/L0)*(1.0-cos(CA0))+cos(CA0)+2.0/L0*R0*(CA0-sin(CA0))) );
		p++;



		// mnm.prt();				//エラー表示
		mnm.prt_sum();			//エラーの合計を表示

	}
	}
	
		fout<<hmin<<",";
	
	R0b=R0;
	
	//Wetting ratio at the breaking condition from this model
	X0=pow(hmin/R0, 3.0)*1.0/pow(1-cos(CA0), 3.0)-2.0*R0*f/L0/pow(1-cos(CA0), 3.0)+2.0*R0*sin(CA0)/L0;
		fout<<X0<<",";
		fout<<R0<<",";

	X0b=X0;

	//Flat portion of the rivulet
	Z0=X0*L0/2.0-R0*sin(CA0);
		fout<<Z0<<",";

	//Dimensionless group for comparison with experiment from Maron 1982
	dp=(rho*pow(st,3.0))/(pow(myu,4.0)*gl);
		fout<<dp<<",";

	//Dimensionless minimum wetting rate for comparison with experiment from Maron 1982
	G=pow(hmin, 3.0)/pow(6.0*myu/(rho*rho*gl*sin(beta)),1.0);
		fout<<G/myu<<",";

	//Dimensionless film minimum thickness
	h0p=hmin*pow(rho*rho*rho*gl*gl*sin(beta)*sin(beta)/(15.0*myu*myu*st),1.0/5.0);
		fout<<h0p<<",";

	//Rivulet energy at film breaking from this model
	EM0=rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(R0, 5.0)*(+2.0/15.0*R0/L0*g+pow(1-cos(CA0), 2.0)/15.0*(pow(hmin/R0, 3.0)-2.0*R0*f/L0))+st*((pow(hmin/R0, 3.0)/pow(1-cos(CA0), 3.0)-2.0*R0*f/pow(1-cos(CA0), 3.0)/L0+2.0*R0*sin(CA0)/L0)*(1.0-cos(CA0))+cos(CA0)+2.0/L0*R0*(CA0-sin(CA0))) ;
		fout<<EM0<<endl;

		fout<<h0<<",";
		fout<<hmin<<",";

	R0=R0/10;
	//Calculation of Rivuet Radius R at the stable condition from this model/////////////////////////////////////////
				CNewtonRaphsonMethod mnm1;		// CNewtonRaphsonMethodクラスの宣言
	mnm1.setup(2000,1+1e-12,1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

		p=0;

			mnm1.setValue( p , R0);
		p++;

//			mnm1.setValue( p , hmin );
//		p++;

		mnm1.setAcc(0.03);				//加速度勾配の入力　※なくても可

	mnm1.initial();					//計算を開始する前に必ず初期化してください
	for(mnm1.main_loop_init();mnm1.main_loop_check();mnm1.main_loop_reinit()){		// おまじない
		for(mnm1.sub_loop_init();mnm1.sub_loop_check();mnm1.sub_loop_reinit()){	// おまじない

		p=0;
	//値の設定
			R0 = mnm1.getValue(p);
		p++;

//			hmin = mnm1.getValue(p);
//		p++;

		p=0;

	//minimal rivulet energy
			mnm1.setError( p , 0.0, rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*(+4.0/5.0*pow(R0, 5.0)/L0*g+2.0/15.0*pow(h0, 3.0)*R0*pow(1.0-cos(CA0), 2.0)-4.0/5.0*pow(R0, 5.0)/L0*f*pow(1.0-cos(CA0), 2.0))+st*((-3.0*pow(h0, 3.0)/(pow(R0, 4.0)*pow(1.0-cos(CA0), 3.0))-2.0*f/pow(1.0-cos(CA0), 3.0)/L0+2.0*sin(CA0)/L0)*(1.0-cos(CA0))+2.0/L0*(CA0-sin(CA0))) );
		p++;

	//energy balance
//			mnm1.setError( p , 1.0/15.0*rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(hmin, 5.0)+st, rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(R0, 5.0)*(+2.0/15.0*R0/L0*g+pow(1.0-cos(CA0), 2.0)/15.0*(pow(hmin/R0, 3.0)-2.0*R0*f/L0))+st*((pow(hmin/R0, 3.0)/pow(1.0-cos(CA0), 3.0)-2.0*R0*f/pow(1.0-cos(CA0), 3.0)/L0+2.0*R0*sin(CA0)/L0)*(1.0-cos(CA0))+cos(CA0)+2.0/L0*R0*(CA0-sin(CA0))) );
//		p++;



		// mnm.prt();				//エラー表示
		mnm1.prt_sum();			//エラーの合計を表示

	}
	}
	
	//Wetting ratio at the stable condition from this model
	X0=pow(h0/R0, 3.0)*1.0/pow(1-cos(CA0), 3.0)-2.0*R0*f/L0/pow(1-cos(CA0), 3.0)+2.0*R0*sin(CA0)/L0;
		fout<<X0<<",";
		fout<<R0<<",";

	//Flat portion of the rivulet
	Z0=X0*L0/2.0-R0*sin(CA0);
		fout<<Z0<<",";

	//Rivulet energy at the syable condition from this model
	EM0=rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(R0, 5.0)*(+2.0/15.0*R0/L0*g+pow(1-cos(CA0), 2.0)/15.0*(pow(h0/R0, 3.0)-2.0*R0*f/L0))+st*((pow(h0/R0, 3.0)/pow(1-cos(CA0), 3.0)-2.0*R0*f/pow(1-cos(CA0), 3.0)/L0+2.0*R0*sin(CA0)/L0)*(1.0-cos(CA0))+cos(CA0)+2.0/L0*R0*(CA0-sin(CA0))) ;
		fout<<EM0<<endl;

	if(h0>hmin){
		
		for(i=1;i<=N;i++){
		
	X[i]=1.0;
		fout<<X[i]<<endl;

//		CA[i]=CA0;
//		fout<<CA[i]<<",";				
//		fout<<h0<<",";

						}
				}

	else{

//Calculation of the rivulet local wetting evolution (approximate solution of the Lagrange equation)/////////////////////////////////////////
//		for i=1...

		R0=0.5*R0;

				CNewtonRaphsonMethod mnm2;		// CNewtonRaphsonMethodクラスの宣言
	mnm2.setup(2000,1+1e-12,1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

		p=0;
			mnm2.setValue( p , R0);
		p++;

//			mnm2.setValue( p , X0);
//		p++;

		mnm2.setAcc(0.03);				//加速度勾配の入力　※なくても可

	mnm2.initial();					//計算を開始する前に必ず初期化してください
	for(mnm2.main_loop_init();mnm2.main_loop_check();mnm2.main_loop_reinit()){		// おまじない
		for(mnm2.sub_loop_init();mnm2.sub_loop_check();mnm2.sub_loop_reinit()){	// おまじない

		p=0;
			R0 = mnm2.getValue(p);
		p++;

//			X0 = mnm2.getValue(p);
//		p++;

		p=0;
//			mnm2.setError( p , 0.0, rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*(+4.0/5.0*pow(R0, 5.0)/L0*g+2.0/15.0*pow(h0, 3.0)*R0*pow(1.0-cos(CA0), 2.0)-4.0/5.0*pow(R0, 5.0)/L0*f*pow(1.0-cos(CA0), 2.0))+st*((-3.0*pow(h0, 3.0)/(pow(R0, 4.0)*pow(1.0-cos(CA0), 3.0))-2.0*f/pow(1.0-cos(CA0), 3.0)/L0+2.0*sin(CA0)/L0)*(1.0-cos(CA0))+2.0/L0*(CA0-sin(CA0))) );
//		p++;
		mnm2.setError( p , pow(h0/R0,3.0), 2.0*pow(R0,1.0)*f/L0+2.0*(X[1]/2.0-R0*sin(CA0)/L0)*pow((1-cos(CA0)), 3.0) );
		p++;

		// mnm2.prt();				//エラー表示
		mnm2.prt_sum();			//エラーの合計を表示

	}
	}

		fout<<R0<<",";
	
	//Inlet flat part of the rivulet
	Z0=X[1]*L0/2.0-R0*sin(CA0);
		fout<<Z0<<",";

//	X0=pow(h0/R0,3.0)/pow((1-cos(CA0)), 3.0)- 2.0*pow(R0,1.0)*f/L0/pow((1-cos(CA0)), 3.0)+2.0*R0*sin(CA0)/L0;
//		fout<<X0<<",";

	//Inlet energy of the rivulet
	EM[1]=rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(R0, 5.0)*(+2.0/15.0*R0/L0*g+pow(1.0-cos(CA0), 2.0)/15.0*(pow(h0/R0, 3.0)-2.0*R0*f/L0))+st*((pow(h0/R0, 3.0)/pow(1.0-cos(CA0), 3.0)-2.0*R0*f/pow(1.0-cos(CA0), 3.0)/L0+2.0*R0*sin(CA0)/L0)*(1.0-cos(CA0))+cos(CA0)+2.0/L0*R0*(CA0-sin(CA0)));
		fout<<EM[1]<<",";

		fout<<0<<",";

	//Define dX and dz
	dX=-0.0001;
	dz=D/N;

	//Calculate dR
	dR=dX/(2.0*sin(CA0)/L0-3.0*pow(h0/(1.0-cos(CA0)), 3.0)/pow(R0, 4.0)-2.0*f/(L0*pow(1-cos(CA0), 3.0)));
		fout<<dR<<",";

	//Calculate eriv(R-dR)
	EMM[1]=rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow((R0+dR), 5.0)*(+2.0/15.0*(R0+dR)/L0*g+pow(1.0-cos(CA0), 2.0)/15.0*(pow(h0/(R0+dR), 3.0)-2.0*(R0+dR)*f/L0))+st*((pow(h0/(R0+dR), 3.0)/pow(1.0-cos(CA0), 3.0)-2.0*(R0+dR)*f/pow(1.0-cos(CA0), 3.0)/L0+2.0*(R0+dR)*sin(CA0)/L0)*(1.0-cos(CA0))+cos(CA0)+2.0/L0*(R0+dR)*(CA0-sin(CA0)));
		fout<<EMM[1]<<",";

	//Calculate dedX
	dedX=(EM[1]-EMM[1])/dX;
		fout<<dedX<<",";

	//Calculate gamma
	gamma=rho*2.0*(R0*R0/4.0*(CA0-sin(CA0))+Z0*R0*(1-cos(CA0)))/L0;
		fout<<gamma<<",";

	//Calculate uav
	uav=ma/(2.0*(rho*(R0*R0/4.0*(CA0-sin(CA0))+Z0*R0*(1-cos(CA0)))));
		fout<<uav<<",";

	uav0=uav;

	//Approximate Lagrange equation; calculate X(z+dz)
//	X[2]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[1]-X[0];
//	X[2]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[1]-X[0]-L0*myu/(gamma*uav*uav)*(X[1]-X[0]);
//	X[2]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[1]-X[0]-(X[1])*myu/(gamma*uav)*(X[1]-X[0]);
//	X[2]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[1]-X[0]-myu/((X[1])*gamma*uav)*(X[1]-X[0]);
	
//		X[2]=(dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[1]-X[0]-L0*myu/((X[1]-X0)*gamma*uav*uav)*(-X[1]))/(1+L0*myu/((X[1]-X0)*gamma*uav*uav));
//	X[2]=-X[0]+pow(X[0]*X[0]+2.0*(-0.5*(X[1]*X[0]-X[0]*X[0])-(X[1]-X[0])*dedX*dz*dz/(L0*L0*gamma*uav*uav)),0.5);	
//	X[2]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[1]-X[0]-(X[1]-X[0]);
	X[2]=X[1]-dz*pow(2.0*(EM[1]-EM0)/(gamma*uav*uav*L0*L0),0.5);
		fout<<X[2]<<endl;

	//variables initialization for plotting (eriv) vs (X) ??? 
		for(k=1;k<=E;k++){

	RR[k]=2.5*(R0/E*(E-21.0)+R0/E*42.0/double(E)*k);
		fout<<RR[k]<<",";

							}
	fout<<0<<endl;
			
		for(k=1;k<=E;k++){

	XX[k]=pow(h0/RR[k],3.0)/pow((1-cos(CA0)), 3.0)- 2.0*pow(RR[k],1.0)*f/L0/pow((1-cos(CA0)), 3.0)+2.0*RR[k]*sin(CA0)/L0;
		fout<<XX[k]<<",";

							}
	fout<<0<<endl;

	R[0]=R0;
			
	for(i=1;i<=N;i++){

	R[i]=R0;
				
						}

/////////////////////////Position loop///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for(i=2;i<=N;i++){

		fout<<D/double(N)*(i-1)<<",";
//				uave=pow(m/L/X[i],2.0/3.0)*pow(gl*sin(beta)/(rho*3.0*myu),1.0/3.0);
	
////////////Contact angle/////////////////////////////////////////////////////////
////////////////////NewtonRaphson Method//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			CNewtonRaphsonMethod mnm1;		// CNewtonRaphsonMethodクラスの宣言
	mnm1.setup(2000,1+1e-12,1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

		p=0;

			mnm1.setValue( p , R[i] );
		p++;
//			mnm1.setValue( p , CA[i] );
//		p++;
	
		mnm1.setAcc(0.5);				//加速度勾配の入力　※なくても可

	mnm1.initial();					//計算を開始する前に必ず初期化してください
	for(mnm1.main_loop_init();mnm1.main_loop_check();mnm1.main_loop_reinit()){		// おまじない
		for(mnm1.sub_loop_init();mnm1.sub_loop_check();mnm1.sub_loop_reinit()){	// おまじない

		p=0;
			//値の設定
			
			R[i] = mnm1.getValue(p);
		p++;
//			CA[i] = mnm1.getValue(p);
//		p++;

		p=0;

			mnm1.setError( p , pow(h0/R[i],3.0), (2.0*pow(R[i],1.0)*(-1.0/4.0*pow(cos(CA0),3.0)*sin(CA0) -13.0/8.0*cos(CA0)*sin(CA0)-3.0/2.0*CA0*pow(sin(CA0),2.0) +15.0/8.0*CA0)/L0+2.0*(X[i]/2.0-R[i]*sin(CA0)/L0)*pow((1-cos(CA0)), 3.0)) );
		p++;
//			mnm1.setError( p , 0.0, rho*rho*rho*gl*gl/myu/myu*(+4.0/5.0*pow(R0, 5.0)/L0*(((CA[i]*(5.0/16.0+15.0/4.0*pow(cos(CA[i]),2.0) +5.0/2.0*pow(cos(CA[i]),4.0))-sin(CA[i])*(113.0/48.0*cos(CA[i]) +97.0/24.0*pow(cos(CA[i]),3.0) +1.0/6.0*pow(cos(CA[i]),5.0)))))+2.0/15.0*pow(h0, 3.0)*R0*pow(1.0-cos(CA0), 2.0)-4.0/5.0*pow(R0, 5.0)/L0*((-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])-3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i]))*pow(1.0-cos(CA0), 2.0))+st*((-3.0*pow(h0, 3.0)/(pow(R0, 4.0)*pow(1.0-cos(CA0), 3.0))-2.0*((-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])-3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i]))/pow(1.0-cos(CA0), 3.0)/L0+2.0*sin(CA0)/L0)*(1.0-cos(CA0))+2.0/L0*(CA0-sin(CA0))));
//		p++;
//			mnm1.setError( p , 0.0, rho*rho*rho*gl*gl/myu/myu*pow(R[i], 4.0)*(-4.0/5.0*R[i]*((CA[i]*(5.0/16.0+15.0/4.0*pow(cos(CA[i]),2.0)+5.0/2.0*pow(cos(CA[i]),4.0))-sin(CA[i])*(113.0/48.0*cos(CA[i])+97.0/24.0*pow(cos(CA[i]),3.0)+1.0/6.0*pow(cos(CA[i]),5.0))))+1.0/3.0*X[i]*L0*pow(1.0-cos(CA[i]), 5.0)-4.0/5.0*R[i]*sin(CA[i])*pow(1.0-cos(CA[i]), 5.0))+2.0*st*(CA[i]-sin(CA[i])) );
//		p++;

////////////////Continuity equation///////////////////////////////////////////////////////////////////////////////////////////////////////////


		// mnm3.prt();				//エラー表示
		mnm1.prt_sum();			//エラーの合計を表示

	}
	}

		fout<<X[i]<<",";
		fout<<R[i]<<",";
				//fout<<uave<<endl;

	//calculate the flat part of the rivulet
	ZZ[i]=X[i]*L0/2.0-R[i]*sin(CA0);
		fout<<ZZ[i]<<",";

////////////////Rivulet average thickness///////////////////////////////////////////////////////////////////////////////////////////////////////////

//	h[i]=pow(3.0*m/L*myu/(X[i]*rho*rho*gl*sin(beta)),1.0/3.0);

////////////////Rivulet average stream-wise velocity///////////////////////////////////////////////////////////////////////////////////////////////////////////

			uave[i]=0;
			h[0]=R[i]*(1.0-cos(CA0));
			h[0]=R0b*(1.0-cos(CA0));
			t[i]=0;	//???

			for(k=1;k<=Z;k++){

				if(k*(X0b*L0/2.0)/double(Z)<=R0b*sin(CA0)) {

					teta=asin((k*X0b*L0/2.0)/R0b/double(Z));
				
					h[k]=(h[k-1]+R0b*(cos(teta)-cos(CA0)))/2.0;

					}
			
				else{

					h[k]=R0b*(1.0-cos(CA0));

				}

				for(j=1;j<=M;j++){

					u[k][j]=(rho*gl*h[k]*h[k]*sin(beta)/myu)*(double(j)/double(M)-1.0/2.0*pow(double(j)/double(M), 2.0));
					uave[i] +=u[k][j]/double(M)/double(Z);

					}

					t[i] +=h[k]/double(Z);
					fout<<h[k]<<",";
				}

		fout<<uave[i]<<",";
		
////////////////Rivulet Energy/Calculated quantities//////////////////////////////////////////////////////////////////////////////////////////////////////////

	EM[i]=rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(R[i], 5.0)*(+2.0/15.0*R[i]/L0*(((CA0*(5.0/16.0+15.0/4.0*pow(cos(CA0),2.0) +5.0/2.0*pow(cos(CA0),4.0))-sin(CA0)*(113.0/48.0*cos(CA0) +97.0/24.0*pow(cos(CA0),3.0) +1.0/6.0*pow(cos(CA0),5.0)))))+pow(1.0-cos(CA0), 2.0)/15.0*(pow(h0/R[i], 3.0)-2.0*R[i]*((-1.0/4.0*pow(cos(CA0),3.0)*sin(CA0) -13.0/8.0*cos(CA0)*sin(CA0)-3.0/2.0*CA0*pow(sin(CA0),2.0) +15.0/8.0*CA0))/L0))+st*((pow(h0/R[i], 3.0)/pow(1.0-cos(CA0), 3.0)-2.0*R[i]*((-1.0/4.0*pow(cos(CA0),3.0)*sin(CA0) -13.0/8.0*cos(CA0)*sin(CA0)-3.0/2.0*CA0*pow(sin(CA0),2.0) +15.0/8.0*CA0))/pow(1.0-cos(CA0), 3.0)/L0+2.0*R[i]*sin(CA0)/L0)*(1.0-cos(CA0))+cos(CA0)+2.0/L0*R[i]*(CA0-sin(CA0)));
		fout<<EM[i]<<",";

	//Calculate dR
	dR=dX/(2.0*sin(CA0)/L0-3.0*pow(h0/(1.0-cos(CA0)), 3.0)/pow(R[i], 4.0)-2.0*f/(L0*pow(1-cos(CA0), 3.0)));
		fout<<dR<<",";

	//Calculate eriv(R-dR)
	EMM[i]=rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow((R[i]+dR), 5.0)*(+2.0/15.0*(R[i]+dR)/L0*g+pow(1.0-cos(CA0), 2.0)/15.0*(pow(h0/(R[i]+dR), 3.0)-2.0*(R[i]+dR)*f/L0))+st*((pow(h0/(R[i]+dR), 3.0)/pow(1.0-cos(CA0), 3.0)-2.0*(R[i]+dR)*f/pow(1.0-cos(CA0), 3.0)/L0+2.0*(R[i]+dR)*sin(CA0)/L0)*(1.0-cos(CA0))+cos(CA0)+2.0/L0*(R[i]+dR)*(CA0-sin(CA0)));
		fout<<EMM[i]<<",";

	//Calculate dedX
	dedX=(EM[i]-EMM[i])/dX;
		fout<<dedX<<",";

	gamma0=gamma;

	//Calculate gamma
	gamma=rho*2.0*(R[i]*R[i]/4.0*(CA0-sin(CA0))+ZZ[i]*R[i]*(1-cos(CA0)))/L0;
		fout<<gamma<<",";

	uav0=uav;

	//Calculate uav
	uav=ma/(2.0*(rho*(R[i]*R[i]/4.0*(CA0-sin(CA0))+ZZ[i]*R[i]*(1-cos(CA0)))));
		fout<<uav<<",";

	//Ti
	TI=0.5*L0*L0*gamma*pow(uav*(X[i]-X[i-1])/dz,2.0);

	//Ti-1
	TII=0.5*L0*L0*gamma0*pow(uav0*(X[i-1]-X[i-2])/dz,2.0);

	//Approximate Lagrange equation; calculate X(z+dz)
//	X[i+1]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[i]-X[i-1];
//	X[i+1]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[i]-X[i-1]-L0*myu/(gamma*uav*uav)*(X[i]-X[i-1]);	
//	X[i+1]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[i]-X[i-1]-(X[i])*myu/(gamma*uav)*(X[i]-X[i-1]);	
	X[i+1]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[i]-X[i-1]-myu/((X[i])*gamma*uav)*(X[i]-X[i-1]);	

//	X[i+1]=dedX*dz*dz/(L0*L0*gamma*uav*uav)-(TI-TII)/(X[i]-X[i-1])*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[i]-X[i-1]-myu/((X[i])*gamma*uav)*(X[i]-X[i-1]);	
//		X[i+1]=(dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[i]-X[i-1]-L0*myu/((X[i]-X0)*gamma*uav*uav)*(-X[i]))/(1+L0*myu/((X[i]-X0)*gamma*uav*uav));	
//	X[i+1]=-X[i-1]+pow(X[i-1]*X[i-1]+2.0*(-0.5*(X[i]*X[i-1]-X[i-1]*X[i-1])-(X[i]-X[i-1])*dedX*dz*dz/(L0*L0*gamma*uav*uav)),0.5);	
//	X[i+1]=dedX*dz*dz/(L0*L0*gamma*uav*uav)+2.0*X[i]-X[i-1]-(X[i]-X[i-1]);
//	X[i+1]=X[i]-dz*pow(2.0*(EM[i]-EM0)/(gamma*uav*uav*L0*L0),0.5);
//		fout<<X[i+1]<<",";



	fout<<0<<",";

		for(k=1;k<=E;k++){	

	EMM[k]=rho*rho*rho*gl*gl*sin(beta)*sin(beta)/myu/myu*pow(RR[k], 5.0)*(+2.0/15.0*RR[k]/L0*(((CA0*(5.0/16.0+15.0/4.0*pow(cos(CA0),2.0) +5.0/2.0*pow(cos(CA0),4.0))-sin(CA0)*(113.0/48.0*cos(CA0) +97.0/24.0*pow(cos(CA0),3.0) +1.0/6.0*pow(cos(CA0),5.0)))))+pow(1.0-cos(CA0), 2.0)/15.0*(pow(h0/RR[k], 3.0)-2.0*RR[k]*((-1.0/4.0*pow(cos(CA0),3.0)*sin(CA0) -13.0/8.0*cos(CA0)*sin(CA0)-3.0/2.0*CA0*pow(sin(CA0),2.0) +15.0/8.0*CA0))/L0))+st*((pow(h0/RR[k], 3.0)/pow(1.0-cos(CA0), 3.0)-2.0*RR[k]*((-1.0/4.0*pow(cos(CA0),3.0)*sin(CA0) -13.0/8.0*cos(CA0)*sin(CA0)-3.0/2.0*CA0*pow(sin(CA0),2.0) +15.0/8.0*CA0))/pow(1.0-cos(CA0), 3.0)/L0+2.0*RR[k]*sin(CA0)/L0)*(1.0-cos(CA0))+cos(CA0)+2.0/L0*RR[k]*(CA0-sin(CA0)));
		fout<<EMM[k]<<",";

							}
	fout<<0<<endl;

////////////////Local wetting ratio///////////////////////////////////////////////////////////////////////////////////////////////////////////	
		
//			X[i+1]=X[i]-pow(2.0*t[i]*(EM[i]-EM[0])/(L0*L0*m/L*uave[i]*uave[i]),0.5)*D/double(N);	
//				fout<<X[i+1]<<",";

	}

	

	} //else

//	fout<<0<<endl;

	Xave=0;

	for(i=1;i<=N;i++){

				Xave+=X[i]/double(N);
				
				}
	fout<<Xave<<",";
	fout<<0<<endl;

	

	}




}*/
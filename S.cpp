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
	int G,l,i,j,N,M,O,n,p,t,o;
	N=60;/////Angle mesh
	M=30;/////Thickness mesh
	O=1;/////Parameter Analysis
	G=1;/////Mass flow rate
	double u[100][100],v[100][100],Ti[100][100],Xi[100][100],Cw[100][100],Et[100][100],Ef[100][100],Ec[100][100],Ed[100][100];
	double g,st, L,rho,myu,beta,hmin,hminp,Re,PI,f,gl,k,Reb,qf,mf,dp,Tw,d,hv,ro,cp,Tsat,m,P,_Tw,Gw,Twi,ri,aave,bave,Xave,Tave,uave,a,htc,mtc,Tin,Xin,B,S,C,Mw,hw0,sw0,Cwi;
	double Gv[100],A[100],b[100],X[100],CA[100],h[100],dr[100],r[100],R[100],dtx[100],dty[100],dxx[100],dxy[100],dux[100],duy[100],dvx[100],dvy[100];

	ofstream fout("S.csv");
	if(!fout){
	cout << "ファイルをオープンできませんでした"<< endl;
    return 1;
	}
	else
		cout<< "ファイルをオープンしました" << endl;

/////////////////////////flowrate loop////////////////////////////////////////////////////////////////////////////////////////////////////////////
for(l=1;l<=G;l++){

	m=0.09;///[kg/s]
	fout<<m<<",";

////////////////////Experiment conditions////////////////////

	P=12.5;//[kPa]
	t=10;//tube number
    L=0.878;//tube length[m]
	//ri=8.0/1000.0;
	PI=(6*asin(0.5));//π
	beta=PI/N;//Angle of tube//////for 180 degree

////////////////////grid generation in the radial direction///

	for(j=1;j<=M;j++){
	r[j]=1.0/2.0*(1.0-cos((j-1)/(double(M)-1)*PI));
						}
	for(j=1;j<=M-1;j++){
	dr[j]=r[j+1]-r[j];
						}
	for(j=1;j<=M-2;j++){
	R[j]=dr[j+1]/dr[j];
						}
////////////////////first conditions/////////////////////////

	//Twi=32.0;//water inlet temperature
	Tin=95.0;//Solution inlet temperature [^C]
	Xin=0.60;//Solution inlet concentration
	Tave=Tin;
	Xave=Xin;
	//Gw=1.0;///kg/s
/////////////////////properties calculation//////////////////

	st=libr.sc_st_XT(Xave,Tave);//////[N/m]
	rho=libr.sc_rho_XT(Xave,Tave);////[kg/m3]
	myu=libr.sc_visc_XT(Xave,Tave);///[Pas]	
	gl=9.81;//[m/s2]
	k=libr.sc_thc_XT(Xave,Tave);////[W/mK]
	d=libr.sc_d_XT(Xave,Tave);//////[m2/s]
	cp=libr.sc_cp_XT(Xave,Tave)*1000.0;////[J/kgK]
	a=k/rho/cp;
	Tsat=wat.p_t( P );
	Mw=0.018015;/////Molar weight [kg/mol]
	hw0=-241830;/////standard water molar enthalpy [J/mol]
	sw0=188.84;//////standard water molar entropy [J/molK]
	Cwi=(rho)/Mw*(1.0-Xin);

///////////////////variables initialization///////////////////

	for(i=1;i<=N-1;i++){
		Gv[i]=0.001;//interface mass flux per unit surface initialization [kg/m2s]
		for(j=1;j<=M;j++){
    		Xi[i][j]=Xin;//initialising concentration
			Ti[i][j]=Tin;//initialising temperature
							}
						}

//////////////////////Parametric Analysis loop//////////////////////////////////////////////////////////////////////////////////////////////
	for (o=1;o<=O;o++){
		
	Tw=85.0;//Tube Wall temperature[^C]
	fout<<Tw<<",";
	ro=9.0/1000.0*o;//outer tube diameter [m]
	fout<<ro<<",";

//////////////////////Partial wetting model Parameters////////

	Re=4.0*m/L/myu;
		fout<<Re<<",";
	dp=(rho*pow(st,3.0))/(pow(myu,4.0)*gl);//dimensionless parameter group for minimum stable thickness calculation
	Reb=4.0*8.5027*pow(dp,0.0505);//Break-up Reynolds from Maron1982
	hmin=pow(3.0*Reb/4.0*pow(myu/rho,2.0)/9.81,1.0/3.0);//minimum thickness [m] from Maron1982

//////////////////////Average tube values/////////////////////
		
		aave=0;
		bave=0;
		qf=0;
		mf=0;

////////////////////////Angle loop//////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(i=1;i<=N-1;i++){	

	h[i]=pow(3.0*m/L*myu/(rho*rho*gl*sin(beta*i)),1.0/3.0);//uniform film thickness

////////////////////////Wetting model///////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////first half of tube surface////////////
	if(i<N/2){

		hminp=pow(((rho*rho*rho)*gl*gl*sin(beta*i)*sin(beta*i))/(15.0*myu*myu*st),0.2)*hmin;// dimensionless critical thickness
		CA[i]= 1665.4*pow(hminp,6.0) - 3446.8*pow(hminp,5.0) + 2399.3*pow(hminp,4.0) - 459.73*pow(hminp,3.0) + 49.538*pow(hminp,2.0) + 2.9741*hminp + 0.8616;//Local contact Angle
		CA[i]=CA[i]/180*PI;//contact angle from critical condition and minimum stable thickness
		
		if(h[i]>=hmin){//complete wetting/////////////

			X[i]=1.0;
			b[i]=h[i];
						}
		else{////////////Partial wetting//////////////

				f=-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])
					 -3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i];//f(θo)
				g=(CA[i]*(5.0/16.0+15.0/4.0*pow(cos(CA[i]),2.0) +5.0/2.0*pow(cos(CA[i]),4.0))
					 -sin(CA[i])*(113.0/48.0*cos(CA[i]) +97.0/24.0*pow(cos(CA[i]),3.0) +1.0/6.0*pow(cos(CA[i]),5.0)));//Ψ(θo))

				X[i]=pow((h[i]),3.0)*sin(CA[i])/f*pow(2.0/45.0*(pow(rho,3.0)
					 *pow((9.81*sin(beta*i)),2.0)/(st*pow(myu,2.0)))*g/sin(CA[i])*pow((CA[i]/sin(CA[i])-cos(CA[i])),-1),3.0/5.0);//wetting ratio calculation for 10th tube
				
				if(X[i-1]==1){///Linear transition zone between uniform and rivulets configurations///////
						B=-15.0/PI;
						C=15.0/double(N)*(i-1)+1;
						S=B*PI*double(i)/double(N)+C;
						}
						else{
							S=B*PI*double(i)/double(N)+C;
							if(S<0){
							S=0;
							}
							else{}
						}

				X[i]=S*(1-X[i])+X[i];
				A[i]=(log(X[i]+(1-X[i])*pow(Re/Reb,1.0))-log(X[i]))/9.0;//exponential coefficient describing wetting ratio increase for previous tubes
				X[i]=(X[i]+(1-X[i])*pow(Re/Reb,1.0))*exp(-A[i]*(t-1));//wetting ratio calculation for tube t
				b[i]=h[i]/pow(X[i],1.0/3.0);//average rivulet thickness
						
			}
				}////////////first half of tube surface////////////

////////////////////second half of tube surface////////////////////
	else{

		CA[i]= CA[i-1];

		if(X[N/2-1]==1){///Constant contact angle criteria break-up///////
			hminp=0.1593*log(CA[i])+0.6699;
			hmin=pow(((rho*rho*rho)*gl*gl*sin(beta*i)*sin(beta*i))/(15.0*myu*myu*st),-0.2)*hminp;

			if(h[i]>=hmin){///////complete wetting////////////////////////
			X[i]=1.0;
			b[i]=h[i];
							}
			else{/////////////////partial wetting/////////////////////////
	
			f=-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])
				 -3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i];//f(θo)
			g=(CA[i]*(5.0/16.0+15.0/4.0*pow(cos(CA[i]),2.0) +5.0/2.0*pow(cos(CA[i]),4.0))
				 -sin(CA[i])*(113.0/48.0*cos(CA[i]) +97.0/24.0*pow(cos(CA[i]),3.0) +1.0/6.0*pow(cos(CA[i]),5.0)));//Ψ(θo))

			X[i]=pow((h[i]),3.0)*sin(CA[i])/f*pow(2.0/45.0*(pow(rho,3.0)
				 *pow((9.81*sin(beta*i)),2.0)/(st*pow(myu,2.0)))*g/sin(CA[i])*pow((CA[i]/sin(CA[i])-cos(CA[i])),-1),3.0/5.0);
			
				if(X[i-1]==1){///Linear transition zone between uniform and rivulets configurations///////
					B=-15.0/PI;
					C=15.0/double(N)*(i-1)+1;
					S=B*PI*double(i)/double(N)+C;		
								}					
				else{
					S=B*PI*double(i)/double(N)+C;
					if(S<0){
						S=0;
							}
					else{}		
						}
						
			X[i]=S*(1-X[i])+X[i];
			A[i]=(log(X[i]+(1-X[i])*pow(Re/Reb,1.0))-log(X[i]))/9.0;//exponential coefficient describing wetting ratio increase for previous tubes
			X[i]=(X[i]+(1-X[i])*pow(Re/Reb,1.0))*exp(-A[i]*(t-1));//wetting ratio calculation for tube t
			b[i]=h[i]/pow(X[i],1.0/3.0);//average rivulet thickness
									
				}
			}

///////////////////////dynamic wetting ratio calculation (without minimum stable thickness comparison)////////////////////////////////////////////
		else{
	
			f=-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])
				 -3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i];//f(θo)
			g=(CA[i]*(5.0/16.0+15.0/4.0*pow(cos(CA[i]),2.0) +5.0/2.0*pow(cos(CA[i]),4.0))
				 -sin(CA[i])*(113.0/48.0*cos(CA[i]) +97.0/24.0*pow(cos(CA[i]),3.0) +1.0/6.0*pow(cos(CA[i]),5.0)));//Ψ(θo))

			X[i]=pow((h[i]),3.0)*sin(CA[i])/f*pow(2.0/45.0*(pow(rho,3.0)
				 *pow((9.81*sin(beta*i)),2.0)/(st*pow(myu,2.0)))*g/sin(CA[i])*pow((CA[i]/sin(CA[i])-cos(CA[i])),-1),3.0/5.0);

				if(X[i-1]==1){
					B=-15.0/PI;
					C=15.0/double(N)*(i-1)+1;
					S=B*PI*double(i)/double(N)+C;
								}
				else{
					S=B*PI*double(i)/double(N)+C;
					if(S<0){
						S=0;
							}
					else{}
						}

			X[i]=S*(1-X[i])+X[i];
			A[i]=(log(X[i]+(1-X[i])*pow(Re/Reb,1.0))-log(X[i]))/9.0;//exponential coefficient describing wetting ratio increase for previous tubes
			X[i]=(X[i]+(1-X[i])*pow(Re/Reb,1.0))*exp(-A[i]*(t-1));//wetting ratio calculation for tube t
			b[i]=h[i]/pow(X[i],1.0/3.0);//average rivulet thickness
						
			}	
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
///////////////////////numerical soluyion of energy and species transport equations///////////////////////////////////////////////////////////////

	///////////////////Boundary conditions/////////////////

	u[i][1]=0;
	v[i][1]=0;
	Ti[i][1]=Tw;//wall temperature boundary condition

///////////////////////Velocity field///////////////////////
	for(j=2;j<=M;j++){
		u[i][j]=(rho*gl*h[i]*h[i]*sin(PI*i/double(N))/myu)*(r[j]-1.0/2.0*r[j]*r[j]);
		v[i][j]=-(rho*gl*h[i]*h[i]*r[j]*r[j])/(2.0*myu*ro)*(1.0/PI*(-pow(myu*m/L*PI*PI*PI/(9.0*rho*rho*gl),1.0/3.0)*1.0/(pow(sin(PI*double(i)/double(N)),1.0/3.0)*tan(PI*double(i)/double(N))))*sin(PI*double(i)/double(N))+h[i]*(1.0-r[j]/3.0)*cos(PI*double(i)/double(N)));
	}
//////////////////////heat and mass transfer coefficients at the previous angle///////////////////////////////////////////////////////////////////////
if (i>=2){

	htc= k*(-Ti[i-1][3]-R[1]*(2.0+R[1])*Ti[i-1][1]+(1.0+R[1]*R[1]+2.0*R[1])*Ti[i-1][2])/((1.0+R[1])*dr[2])/h[i-1]/abs(Tave-Tw) ;/////////////////h or b
	//fout<<htc<<",";
	if(X[i-1]==1){
		htc= k*(-Ti[i-1][3]-R[1]*(2.0+R[1])*Ti[i-1][1]+(1.0+R[1]*R[1]+2.0*R[1])*Ti[i-1][2])/((1.0+R[1])*dr[2])/h[i-1]/abs(Tave-Tw) ;/////////////////h or b
					}
	else{
		htc= X[i-1]*k*(-Ti[i-1][3]-R[1]*(2.0+R[1])*Ti[i-1][1]+(1.0+R[1]*R[1]+2.0*R[1])*Ti[i-1][2])/((1.0+R[1])*dr[2])/b[i-1]/abs(Tave-Tw) ;/////////////////h or b
			}
		//fout<<htc<<",";
	mtc= -d/(Xi[i-1][M])*((1+2.0*R[M-2])*Xi[i-1][M]+R[M-2]*R[M-2]*Xi[i-1][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*Xi[i-1][M-1])/((1.0+R[M-2])*dr[M-1])/h[i-1]/abs(Xi[i-1][1]-Xi[i-1][M]) ;/////////////////h or b
	//fout<<mtc*3600<<",";
	if(X[i-1]==1){
		mtc= -d/(Xi[i-1][M])*((1+2.0*R[M-2])*Xi[i-1][M]+R[M-2]*R[M-2]*Xi[i-1][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*Xi[i-1][M-1])/((1.0+R[M-2])*dr[M-1])/h[i-1]/abs(Xi[i-1][1]-Xi[i-1][M]) ;/////////////////h or b
					}
	else{
		mtc= -X[i-1]*d/(Xi[i-1][M])*((1+2.0*R[M-2])*Xi[i-1][M]+R[M-2]*R[M-2]*Xi[i-1][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*Xi[i-1][M-1])/((1.0+R[M-2])*dr[M-1])/b[i-1]/abs(Xi[i-1][1]-Xi[i-1][M]) ;/////////////////h or b
			}
		//fout<<mtc*3600<<",";

/////////////////////Tube average calculation////////////////////////

		mf +=Gv[i]*PI*ro/(double(N)-2.0);
			if (i>=3){
			aave +=htc/(double(N)-3.0);
			bave +=mtc/(double(N)-3.0);
			qf +=PI*ro/(double(N)-3.0)*htc*(Tave-Tw);
			}

////////////////////NewtonRaphson Method//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
	mnm.setup(2000,1+1e-12,1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

	p=0;

	for(j=2;j<=M;j++){
	mnm.setValue( p , Ti[i][j],3.0 );
	p++;
	mnm.setValue( p , Xi[i][j],0.1 );
	p++;
	}
	mnm.setValue( p ,Gv[i],0.05  );
	p++;

	mnm.setAcc(0.3);				//加速度勾配の入力　※なくても可

	mnm.initial();					//計算を開始する前に必ず初期化してください
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// おまじない
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// おまじない
			p=0;
			//値の設定
	for(j=2;j<=M;j++){
			Ti[i][j] = mnm.getValue(p);
			p++;
			Xi[i][j] = mnm.getValue(p);
			p++;
						}
	Gv[i] = mnm.getValue(p);
	p++;

	Xi[i][1]=((1.0+R[1]*R[1]+2.0*R[1])*Xi[i][2]-Xi[i][3])/(R[1]*(2.0+R[1]));/// concentration boundary condition

/////////////////////Absorption heat///////////////////////////////////

	hv = wat.sat_hv( Tsat )*1000.0;//-libr.sc_h_XT(Xi[i][M],Ti[i][M])*1000.0;/////[J/kg]?????check calculation formula

	p=0;

for(j=2;j<=M-1;j++){

			mnm.setError( p , (Ti[i][j]-Ti[i-1][j])/(1.0/double(N))*1000.0 , 1000.0*((PI*ro*a)/(u[i][j]*h[i]*h[i])*((2.0*R[j-1]*Ti[i][j+1]+2.0*R[j-1]*R[j-1]*Ti[i][j-1]-2.0*R[j-1]*(1.0+R[j-1])*Ti[i][j])/(dr[j]*dr[j]*(1.0+R[j-1])))+(r[j]/h[i]*(-pow(myu*m/L*PI*PI*PI/(9.0*rho*rho*gl),1.0/3.0)*1.0/(pow(sin(PI*double(i)/double(N)),1.0/3.0)*tan(PI*double(i)/double(N))))-PI*ro*v[i][j]/(h[i]*u[i][j]))*((Ti[i][j+1]-R[j-1]*R[j-1]*Ti[i][j-1]-(1.0-R[j-1]*R[j-1])*Ti[i][j])/((1.0+R[j-1])*dr[j]))) );
			p++;
////////////////Energy transport equation///////////////////////////////////////////////////////////////////////////////////////////////////////////		
		
			mnm.setError( p ,1000.0*(Xi[i][j]-Xi[i-1][j])/(1.0/double(N)) ,1000.0*((PI*ro*d)/(u[i][j]*h[i]*h[i])*((2.0*R[j-1]*Xi[i][j+1]+2.0*R[j-1]*R[j-1]*Xi[i][j-1]-2.0*R[j-1]*(1.0+R[j-1])*Xi[i][j])/(dr[j]*dr[j]*(1.0+R[j-1])))+(r[j]/h[i]*(-pow(myu*m/L*PI*PI*PI/(9.0*rho*rho*gl),1.0/3.0)*1.0/(pow(sin(PI*double(i)/double(N)),1.0/3.0)*tan(PI*double(i)/double(N))))-PI*ro*v[i][j]/(h[i]*u[i][j]))*((Xi[i][j+1]-R[j-1]*R[j-1]*Xi[i][j-1]-(1.0-R[j-1]*R[j-1])*Xi[i][j])/((1.0+R[j-1])*dr[j]))));
			p++;
////////////////Species transport equation//////////////////////////////////////////////////////////////////////////////////////////////////////////		

	////////////Interface////////////////////////////////////////////////////////
	if(j==M-1){

	////////////Phases Equilibrium///////////////////////////////////////////////
			mnm.setError( p , Ti[i][M] , libr.sc_T_XTsat(Xi[i][M],wat.p_t( P )) );
			p++;
		
	////////////Fick's diffusion law/////////////////////////////////////////////
			mnm.setError( p , Gv[i] , -rho*d/(Xi[i][M])*((1+2.0*R[j-1])*Xi[i][j+1]+R[j-1]*R[j-1]*Xi[i][j-1]-(1.0+R[j-1]*R[j-1]+2.0*R[j-1])*Xi[i][j])/((1.0+R[j-1])*dr[j])/h[i]) ;//////h or b
			p++;
		
	///////////Absorption heat Thermal Diffusion/////////////////////////////////		
			mnm.setError( p , Gv[i]*hv , k*((1+2.0*R[j-1])*Ti[i][j+1]+R[j-1]*R[j-1]*Ti[i][j-1]-(1.0+R[j-1]*R[j-1]+2.0*R[j-1])*Ti[i][j])/((1.0+R[j-1])*dr[j])/h[i] );/////////////////h or b
			p++;

				}	
			}

		// mnm.prt();				//エラー表示
		mnm.prt_sum();			//エラーの合計を表示

	}
	}

////////////////Temperature and concentration bulk values/////////////

		uave=0;
		for(n=2;n<=M;n++){
			uave+=u[i][n]*dr[n-1];
							}
		Tave=0;
		Xave=0;
		for(n=2;n<=M;n++){
			Tave+=u[i][n]*Ti[i][n]*dr[n-1]/uave;//average temperature for previous angle step
			Xave+=u[i][n]*Xi[i][n]*dr[n-1]/uave;//average concentration for previous angle step
		}
		//fout<<Xave<<",";
		//fout<<Tave<<endl;

///////////////Molar concentration///////////////////////
		
		for(j=1;j<=M;j++){
		Cw[i][j]=(rho)/Mw*(1.0-Xi[i][j])-Cwi;
							}
/////////////////////Velocity Temperature and concentration Gradients calculation////////////////////////////////////////////////////////////////////////
/////////////////////wall derivatives////////////////////
	dtx[1]=0;
	dty[1]=(-Ti[i][3]-R[1]*(2.0+R[1])*Ti[i][1]+(1.0+R[1]*R[1]+2.0*R[1])*Ti[i][2])/((1.0+R[1])*dr[2])/h[i];
	dux[1]=0;
	dvx[1]=0;
	duy[1]=(-u[i][3]-R[1]*(2.0+R[1])*u[i][1]+(1.0+R[1]*R[1]+2.0*R[1])*u[i][2])/((1.0+R[1])*dr[2])/h[i];
	dvy[1]=(-v[i][3]-R[1]*(2.0+R[1])*v[i][1]+(1.0+R[1]*R[1]+2.0*R[1])*v[i][2])/((1.0+R[1])*dr[2])/h[i];
	dxy[1]=0;
	dxx[1]=1.0/(PI*ro)*(Cw[i][1]-Cw[i-1][1])/(1.0/double(N));
/////////////////////Inner nodes///////////////////////////////////////
	for(j=2;j<=M-1;j++){
	dtx[j]=1.0/(PI*ro)*(Ti[i][j]-Ti[i-1][j])/(1.0/double(N));
	dty[j]=1.0/h[i]*(Ti[i][j+1]-R[j-1]*R[j-1]*Ti[i][j-1]-(1.0-R[j-1]*R[j-1])*Ti[i][j])/((1.0+R[j-1])*dr[j]);
	dxx[j]=1.0/(PI*ro)*(Cw[i][j]-Cw[i-1][j])/(1.0/double(N));
	dxy[j]=1.0/h[i]*(Cw[i][j+1]-R[j-1]*R[j-1]*Cw[i][j-1]-(1.0-R[j-1]*R[j-1])*Cw[i][j])/((1.0+R[j-1])*dr[j]);
	dux[j]=1.0/(PI*ro)*(u[i][j]-u[i-1][j])/(1.0/double(N));
	duy[j]=1.0/h[i]*(u[i][j+1]-R[j-1]*R[j-1]*u[i][j-1]-(1.0-R[j-1]*R[j-1])*u[i][j])/((1.0+R[j-1])*dr[j]);
	dvx[j]=1.0/(PI*ro)*(v[i][j]-v[i-1][j])/(1.0/double(N));
	dvy[j]=1.0/h[i]*(v[i][j+1]-R[j-1]*R[j-1]*v[i][j-1]-(1.0-R[j-1]*R[j-1])*v[i][j])/((1.0+R[j-1])*dr[j]);
	}
////////////////////Interface//////////////////////////////////////////
	dtx[M]=1.0/(PI*ro)*(Ti[i][M]-Ti[i-1][M])/(1.0/double(N));
	dty[M]=1.0/h[i]*((1+2.0*R[M-2])*Ti[i][M]+R[M-2]*R[M-2]*Ti[i][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*Ti[i][M-1])/((1.0+R[M-2])*dr[M-1]);
	dxx[M]=1.0/(PI*ro)*(Cw[i][M]-Cw[i-1][M])/(1.0/double(N));
	dxy[M]=1.0/h[i]*((1+2.0*R[M-2])*Cw[i][M]+R[M-2]*R[M-2]*Cw[i][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*Cw[i][M-1])/((1.0+R[M-2])*dr[M-1]);
	dux[M]=1.0/(PI*ro)*(u[i][M]-u[i-1][M])/(1.0/double(N));
	duy[M]=1.0/h[i]*((1+2.0*R[M-2])*u[i][M]+R[M-2]*R[M-2]*u[i][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*u[i][M-1])/((1.0+R[M-2])*dr[M-1]);
	dvx[M]=1.0/(PI*ro)*(v[i][M]-v[i-1][M])/(1.0/double(N));
	dvy[M]=1.0/h[i]*((1+2.0*R[M-2])*v[i][M]+R[M-2]*R[M-2]*v[i][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*v[i][M-1])/((1.0+R[M-2])*dr[M-1]);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////Entropy generation groups/////////////////////////////////////////////////////////////////////////////////////////////////////
	for(j=1;j<=M;j++){

		Et[i][j]=(libr.sc_thc_XT(Xi[i][j],Ti[i][j]))/(Ti[i][j]*Ti[i][j])*(dtx[j]*dtx[j]+dty[j]*dty[j]);
		Ef[i][j]=(libr.sc_visc_XT(Xi[i][j],Ti[i][j]))/Ti[i][j]*(2.0*(dux[j]*dux[j]+dvy[j]*dvy[j])+(duy[j]+dvx[j])*(duy[j]+dvx[j]));
		Ec[i][j]=(wat.sat_cpv(Ti[i][j])*1000*Mw*(Ti[i][j]-298.15)+hw0+Ti[i][j]*(wat.sat_cpv(Ti[i][j])*1000*Mw*log(Ti[i][j]/298.15)+sw0))*(v[i][j]*dty[j]+u[i][j]*dtx[j])*(Cw[i][j])/(Ti[i][j]*Ti[i][j]);
		Ed[i][j]=-(wat.sat_cpv(Ti[i][j])*1000*Mw*(Ti[i][j]-298.15)+hw0+Ti[i][j]*(wat.sat_cpv(Ti[i][j])*1000*Mw*log(Ti[i][j]/298.15)+sw0))*(dxy[j]*dty[j]+dxx[j]*dtx[j])*(libr.sc_d_XT(Xi[i][j],Ti[i][j]))/(Ti[i][j]*Ti[i][j]);

						}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

else {

//////////////Inlet position parameters////////////////////////////////////
	for(j=1;j<=M;j++){

		Cw[i][j]=(rho)/Mw*(1.0-Xi[i][j])-Cwi;
		Et[i][j]=0;
		Ef[i][j]=0;
		Ec[i][j]=0;
		Ed[i][j]=0;

						}

		}
	
	}////////////////////////Angle loop//////////////////////////////////////////////////////////////////////////////////////////////////////////

	fout<<qf<<",";
	fout<<aave<<",";
	fout<<bave*3600<<",";	
	fout<<mf<<endl;

///////////////////////////////////////cooling water

//	hwi=wat.sc_hl(100.15,Twi)*1000.0;
//	hwo=hwi+Qsum/Gw;
//	Two=wat.T_Ph(100.15,hwo/1000.0);
//	Rew=Gw/(wat.sc_myul(100.15,Twi)*PI*ri);
//	aw=0.023*pow(wat.sc_Prl(100.15,Two),0.4)*pow(Rew,0.8)*wat.sc_laml(100.15,Twi)/(ri*2.0);
//	K=1/(ri*2.0)/(1/aw/(ri*2.0)+1.0/2.0/400.0*log(ro/ri));	
//	_Tw=(Twi*(exp(K/Qsum*2.0*PI*ri*L*(Twi-Two)))-Two)/((exp(K/Qsum*2.0*PI*ri*L*(Twi-Two)))-1);

//////////////////////////Results printing out//////////////////////////////////////////////////////////////////////////////////////////////

	for(j=1;j<=M;j++){
	for(i=1;i<=N-1;i++){
		fout<<Et[i][j]<<",";
						}
	fout<<0<<endl;
						}

		for(j=1;j<=M;j++){
	for(i=1;i<=N-1;i++){
		fout<<Ef[i][j]<<",";
						}
	fout<<0<<endl;
						}

			for(j=1;j<=M;j++){
	for(i=1;i<=N-1;i++){
		fout<<Ec[i][j]<<",";
						}
	fout<<0<<endl;
						}

				for(j=1;j<=M;j++){
	for(i=1;i<=N-1;i++){
		fout<<Ed[i][j]<<",";
						}
	fout<<0<<endl;
						}

					for(j=1;j<=M;j++){
	for(i=1;i<=N-1;i++){
		fout<<Et[i][j]+Ef[i][j]+Ec[i][j]+Ed[i][j]<<",";
						}
	fout<<0<<endl;
						}

	//for(i=1;i<=N-1;i++){
		//fout<<Gv[i]<<endl;
		//fout<<X[i]<<",";
	//if(X[i]==1){
		//	fout<<h[i]<<endl;
	//}
	//else {
		//	fout<<b[i]<<endl;
	//}
	//}
	
//for(j=1;j<=M;j++){
	//for(i=1;i<=N-1;i++){
		//fout<<Xi[i][j]<<",";
//	}
	//fout<<0<<endl;
//}

//for(j=1;j<=M;j++){
	//for(i=1;i<=N-1;i++){
		//fout<<Ti[i][j]<<",";
	//}
	//fout<<0<<endl;
//}

//for(j=1;j<=M;j++){
	//for(i=1;i<=N-1;i++){
		//fout<<u[i][j]<<",";
	//}
	//fout<<0<<endl;
//}

//for(j=1;j<=M;j++){
	//for(i=1;i<=N-1;i++){
		//fout<<v[i][j]<<",";
	//}
	//fout<<0<<endl;
//}
/////////////////////////////Results Printing out///////////////////////////////////////////////////////////////////////////////////////////
	
}////////////////////////////Parametric Analysis loop///////////////////////////////////////////////////////////////////////////////////////
		
}////////////////////////////flowrate loop//////////////////////////////////////////////////////////////////////////////////////////////////

}*/
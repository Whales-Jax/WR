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
	int i,j,N,M,O,n,p,t,o;
	N=30;/////Angle's mesh
	M=15;/////Thickness's mesh
	O=1;/////Mass flow rate
	double st, L,rho,myu,beta,hmin,hminp,Re,h[100],g,CA[100],Ti[100][100],Xi[100][100],PI,X[100],f,gl,k,A[100],Reb,dp,b[100],Tw,d,hv,ro,cp,Tsat,m,P,Gv[100],_Tw,Gw,Twi,ri,Xave,Tave,uave,a,htc,mtc,Tin,Xin,B,S,C;
	double dr[100],r[100],R[100],u[100][100],v[100][100];

	ofstream fout("GridWRCP.csv");
	if(!fout){
	cout << "ファイルをオープンできませんでした"<< endl;
    return 1;
	}
	else
		cout<< "ファイルをオープンしました" << endl;


////////////////////Experiment condition////////////////////
	P=12.5;//[kPa]
	Tw=83.54;//
	Tin=100.0;//
	Xin=0.60;//
///////////////first condition
	Twi=32.0;//water inlet temperature????

	for (o=1;o<=O;o++){//////////////////////flowrate loop
		
	//////flowrate initialization
	m=0.02408*o;///[kg/s]
	
	
	fout<<m<<endl;

///////////////geometric conditions

	t=10;
    L=0.878;//tube length[m]
	ro=9.0/1000.0;
	ri=8.0/1000.0;
	PI=(6*asin(0.5));//π
	beta=PI/N;//Angle of tube//////for 180 degree

	Tave=Tin;
	Xave=Xin;
	
	//Gw=1.0;///kg/s

	for(i=1;i<=N-1;i++){
		for(j=1;j<=M;j++){
    		Xi[i][j]=Xin;//initialising concentration
			Ti[i][j]=Tin;//initialising temperature
		}
	}
	for(i=1;i<=N-1;i++){
		Gv[i]=0.001;//[kg/m2s]
	}
/////////////////////////////////grid generation in the radial direction 

	for(j=1;j<=M;j++){
	r[j]=1.0/2.0*(1.0-cos((j-1)/(double(M)-1)*PI));
						}
	for(j=1;j<=M-1;j++){
	dr[j]=r[j+1]-r[j];
						}
	for(j=1;j<=M-2;j++){
	R[j]=dr[j+1]/dr[j];
						}
//////////////////properties calculation	
	st=libr.sc_st_XT(Xave,Tave);//////[N/m]
	rho=libr.sc_rho_XT(Xave,Tave);////[kg/m3]
	myu=libr.sc_visc_XT(Xave,Tave);///[Pas]	
	gl=9.81;//[m/s2]
	k=libr.sc_thc_XT(Xave,Tave);////[W/mK]
	d=libr.sc_d_XT(Xave,Tave);//////[m2/s]
	cp=libr.sc_cp_XT(Xave,Tave)*1000.0;////[J/kgK]
	a=k/rho/cp;
	Tsat=wat.p_t( P );
	

	dp=(rho*pow(st,3.0))/(pow(myu,4.0)*gl);//dimensionless parameter group for minimum stable thickness calculation
	Reb=4.0*8.5027*pow(dp,0.0505);//Break-up Reynolds from Maron1982
	hmin=pow(3.0*Reb/4.0*pow(myu/rho,2.0)/9.81,1.0/3.0);//minimum thickness [m] from Maron1982
//////////////////////////////////////////////////////Angle loop WR
	for(i=1;i<=N-1;i++){
	
    Re=4.0*m/L/myu;
	h[i]=pow(3.0*m/L*myu/(rho*rho*gl*sin(beta*i)),1.0/3.0);//uniform film thickness

		if(i<N/2){////////////first half of tube surface
			hminp=pow(((rho*rho*rho)*gl*gl*sin(beta*i)*sin(beta*i))/(15.0*myu*myu*st),0.2)*hmin;//partial wetting

				CA[i]= 1665.4*pow(hminp,6.0) - 3446.8*pow(hminp,5.0) + 2399.3*pow(hminp,4.0) - 459.73*pow(hminp,3.0) + 49.538*pow(hminp,2.0) + 2.9741*hminp + 0.8616;//////////////////////////////////////////////////////
				CA[i]=CA[i]/180*PI;//contact angle from critical condition and minimum stable thickness
		if(h[i]>=hmin){
			X[i]=1.0;//complete wetting
			b[i]=h[i];
		}
		else{
			
					f=-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])
					 -3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i];//f(θo)
					g=(CA[i]*(5.0/16.0+15.0/4.0*pow(cos(CA[i]),2.0) +5.0/2.0*pow(cos(CA[i]),4.0))
					 -sin(CA[i])*(113.0/48.0*cos(CA[i]) +97.0/24.0*pow(cos(CA[i]),3.0) +1.0/6.0*pow(cos(CA[i]),5.0)));//Ψ(θo))

						X[i]=pow((h[i]),3.0)*sin(CA[i])/f*pow(2.0/45.0*(pow(rho,3.0)
						 *pow((9.81*sin(beta*i)),2.0)/(st*pow(myu,2.0)))*g/sin(CA[i])*pow((CA[i]/sin(CA[i])-cos(CA[i])),-1),3.0/5.0);//wetting ratio calculation for 10th tube
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
	
	else{////////////////second half of tube surface
		CA[i]= CA[i-1];


		if(X[N/2-1]==1){
			hminp=0.1593*log(CA[i])+0.6699;
		hmin=pow(((rho*rho*rho)*gl*gl*sin(beta*i)*sin(beta*i))/(15.0*myu*myu*st),-0.2)*hminp;

	if(h[i]>=hmin){
		X[i]=1.0;
		b[i]=h[i];
		}
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
		else{//dynamic wetting ratio calculation (without minimum stable thickness comparison)

		
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
	}
/////////
	
/////////
for(i=1;i<=N-1;i++){
	u[i][1]=0;
	v[i][1]=0;
	Ti[i][1]=Tw;//wall temperature boundary condition
}
for(i=1;i<=N-1;i++){
for(j=2;j<=M;j++){
	
	
		u[i][j]=(rho*gl*h[i]*h[i]*sin(PI*i/double(N))/myu)*(r[j]-1.0/2.0*r[j]*r[j]);
	

	
		v[i][j]=-(rho*gl*h[i]*h[i]*r[j]*r[j])/(2.0*myu*ro)*(1.0/PI*(-pow(myu*m/L*PI*PI*PI/(9.0*rho*rho*gl),1.0/3.0)*1.0/(pow(sin(PI*double(i)/double(N)),1.0/3.0)*tan(PI*double(i)/double(N))))*sin(PI*double(i)/double(N))+h[i]*(1.0-r[j]/3.0)*cos(PI*double(i)/double(N)));
	
}}
//	v[i][M]=0;

//////////////////////////////////////////////////////////////move!!!


////////////////////////////////////////////////////////////


	CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
	mnm.setup(2000,1+1e-12,1e-6);					//セットアップ　※引数はニュートン法の変数以上にする

	p=0;

	// 初期値の代入
	for(i=2;i<=N-1;i++){
	for(j=2;j<=M;j++){
	mnm.setValue( p , Ti[i][j],3.0 );
	p++;
	mnm.setValue( p , Xi[i][j],0.1 );
	p++;
	}}
	for(i=2;i<=N-1;i++){
	mnm.setValue( p ,Gv[i],0.05  );
	p++;
	}

	mnm.setAcc(0.3);				//加速度勾配の入力　※なくても可

	mnm.initial();					//計算を開始する前に必ず初期化してください
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// おまじない
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// おまじない
			p=0;
			//値の設定
			for(i=2;i<=N-1;i++){
	for(j=2;j<=M;j++){
			Ti[i][j] = mnm.getValue(p);
			p++;
			Xi[i][j] = mnm.getValue(p);
			p++;
	}}
			for(i=2;i<=N-1;i++){
	Gv[i] = mnm.getValue(p);
	p++;
			}
			for(i=1;i<=N-1;i++){
	Xi[i][1]=((1.0+R[1]*R[1]+2.0*R[1])*Xi[i][2]-Xi[i][3])/(R[1]*(2.0+R[1]));/// concentration boundary condition
			}
	hv = wat.sat_hv( Tsat )*1000.0;//-libr.sc_h_XT(Xi[i][M],Ti[i][M])*1000.0;/////[J/kg]?????check calculation formula

	p=0;

	for(i=2;i<=N-2;i++){
for(j=2;j<=M-1;j++){
	
		
			mnm.setError( p , (Ti[i+1][j]-Ti[i-1][j])/(2.0/double(N))*1000.0 , 1000.0*((PI*ro*a)/(u[i][j]*h[i]*h[i])*((2.0*R[j-1]*Ti[i][j+1]+2.0*R[j-1]*R[j-1]*Ti[i][j-1]-2.0*R[j-1]*(1.0+R[j-1])*Ti[i][j])/(dr[j]*dr[j]*(1.0+R[j-1])))+(r[j]/h[i]*(-pow(myu*m/L*PI*PI*PI/(9.0*rho*rho*gl),1.0/3.0)*1.0/(pow(sin(PI*double(i)/double(N)),1.0/3.0)*tan(PI*double(i)/double(N))))-PI*ro*v[i][j]/(h[i]*u[i][j]))*((Ti[i][j+1]-R[j-1]*R[j-1]*Ti[i][j-1]-(1.0-R[j-1]*R[j-1])*Ti[i][j])/((1.0+R[j-1])*dr[j]))) );
			p++;
		
		
			mnm.setError( p ,1000.0*(Xi[i+1][j]-Xi[i-1][j])/(2.0/double(N)) ,1000.0*((PI*ro*d)/(u[i][j]*h[i]*h[i])*((2.0*R[j-1]*Xi[i][j+1]+2.0*R[j-1]*R[j-1]*Xi[i][j-1]-2.0*R[j-1]*(1.0+R[j-1])*Xi[i][j])/(dr[j]*dr[j]*(1.0+R[j-1])))+(r[j]/h[i]*(-pow(myu*m/L*PI*PI*PI/(9.0*rho*rho*gl),1.0/3.0)*1.0/(pow(sin(PI*double(i)/double(N)),1.0/3.0)*tan(PI*double(i)/double(N))))-PI*ro*v[i][j]/(h[i]*u[i][j]))*((Xi[i][j+1]-R[j-1]*R[j-1]*Xi[i][j-1]-(1.0-R[j-1]*R[j-1])*Xi[i][j])/((1.0+R[j-1])*dr[j]))));
			p++;
		
if(j==M-1){

	mnm.setError( p , Ti[i][M] , libr.sc_T_XTsat(Xi[i][M],wat.p_t( P )) );
    p++;
		
		
			mnm.setError( p , Gv[i] , -rho*d/(Xi[i][M])*((1+2.0*R[j-1])*Xi[i][j+1]+R[j-1]*R[j-1]*Xi[i][j-1]-(1.0+R[j-1]*R[j-1]+2.0*R[j-1])*Xi[i][j])/((1.0+R[j-1])*dr[j])/h[i]) ;//////h or b
			p++;
	
			


				}
	
}
if(i==N-2){
	for(j=2;j<=M-1;j++){
	
		
			mnm.setError( p , (3.0*Ti[i+1][j]-4.0*Ti[i][j]+Ti[i-1][j])/(2.0/double(N))*1000.0 , 1000.0*((PI*ro*a)/(u[i+1][j]*h[i+1]*h[i+1])*((2.0*R[j-1]*Ti[i+1][j+1]+2.0*R[j-1]*R[j-1]*Ti[i+1][j-1]-2.0*R[j-1]*(1.0+R[j-1])*Ti[i+1][j])/(dr[j]*dr[j]*(1.0+R[j-1])))+(r[j]/h[i+1]*(-pow(myu*m/L*PI*PI*PI/(9.0*rho*rho*gl),1.0/3.0)*1.0/(pow(sin(PI*double(i+1)/double(N)),1.0/3.0)*tan(PI*double(i+1)/double(N))))-PI*ro*v[i+1][j]/(h[i+1]*u[i+1][j]))*((Ti[i+1][j+1]-R[j-1]*R[j-1]*Ti[i+1][j-1]-(1.0-R[j-1]*R[j-1])*Ti[i+1][j])/((1.0+R[j-1])*dr[j]))) );
			p++;
		
		
		
			mnm.setError( p ,1000.0*(3.0*Xi[i+1][j]-4.0*Xi[i][j]+Xi[i-1][j])/(2.0/double(N)) ,1000.0*((PI*ro*d)/(u[i+1][j]*h[i+1]*h[i+1])*((2.0*R[j-1]*Xi[i+1][j+1]+2.0*R[j-1]*R[j-1]*Xi[i+1][j-1]-2.0*R[j-1]*(1.0+R[j-1])*Xi[i+1][j])/(dr[j]*dr[j]*(1.0+R[j-1])))+(r[j]/h[i+1]*(-pow(myu*m/L*PI*PI*PI/(9.0*rho*rho*gl),1.0/3.0)*1.0/(pow(sin(PI*double(i+1)/double(N)),1.0/3.0)*tan(PI*double(i+1)/double(N))))-PI*ro*v[i+1][j]/(h[i+1]*u[i+1][j]))*((Xi[i+1][j+1]-R[j-1]*R[j-1]*Xi[i+1][j-1]-(1.0-R[j-1]*R[j-1])*Xi[i+1][j])/((1.0+R[j-1])*dr[j]))));
			p++;
		
if(j==M-1){

	mnm.setError( p , Ti[i+1][M] , libr.sc_T_XTsat(Xi[i+1][M],wat.p_t( P )) );
    p++;
		
		
			mnm.setError( p , Gv[i] , -rho*d/(Xi[i+1][M])*((1+2.0*R[j-1])*Xi[i+1][j+1]+R[j-1]*R[j-1]*Xi[i+1][j-1]-(1.0+R[j-1]*R[j-1]+2.0*R[j-1])*Xi[i+1][j])/((1.0+R[j-1])*dr[j])/h[i+1]) ;//////h or b
			p++;
		
			


				}
	
	}}}
for(i=2;i<=N-1;i++){
		
			mnm.setError( p , Gv[i]*hv , k*((1+2.0*R[M-2])*Ti[i][M]+R[M-2]*R[M-2]*Ti[i][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*Ti[i][M-1])/((1.0+R[M-2])*dr[M-1])/h[i] );/////////////////h or b
			p++;
			
}
		// mnm.prt();				//エラー表示
		mnm.prt_sum();			//エラーの合計を表示

	}
	}


for(i=1;i<=N-1;i++){
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
		fout<<Xave<<",";
		fout<<Tave<<",";

		if(X[i]==1){
			htc= k*(-Ti[i][3]-R[1]*(2.0+R[1])*Ti[i][1]+(1.0+R[1]*R[1]+2.0*R[1])*Ti[i][2])/((1.0+R[1])*dr[2])/h[i]/(Tave-Tw) ;/////////////////h or b
			}
			else{
			htc= X[i]*k*(-Ti[i][3]-R[1]*(2.0+R[1])*Ti[i][1]+(1.0+R[1]*R[1]+2.0*R[1])*Ti[i][2])/((1.0+R[1])*dr[2])/b[i]/(Tave-Tw) ;/////////////////h or b
		
			}
			fout<<htc<<",";
			if(X[i]==1){
			mtc= -d/(Xi[i][M])*((1+2.0*R[M-2])*Xi[i][M]+R[M-2]*R[M-2]*Xi[i][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*Xi[i][M-1])/((1.0+R[M-2])*dr[M-1])/h[i]/(Xi[i][1]-Xi[i][M]) ;/////////////////h or b
			}
			else{
			mtc= -X[i]*d/(Xi[i][M])*((1+2.0*R[M-2])*Xi[i][M]+R[M-2]*R[M-2]*Xi[i][M-2]-(1.0+R[M-2]*R[M-2]+2.0*R[M-2])*Xi[i][M-1])/((1.0+R[M-2])*dr[M-1])/b[i]/(Xi[i][1]-Xi[i][M]) ;/////////////////h or b
		
			}
			fout<<mtc<<endl;

}

		for(i=1;i<=N-1;i++){
		fout<<Gv[i]<<",";
	fout<<X[i]<<",";
	if(X[i]==1){
		fout<<h[i]<<endl;
	}
	else {
		fout<<b[i]<<endl;
	}
	}

for(j=1;j<=M;j++){
	for(i=1;i<=N-1;i++){
		fout<<Xi[i][j]<<",";
	}
	fout<<0<<endl;
}

for(j=1;j<=M;j++){
	for(i=1;i<=N-1;i++){
		fout<<Ti[i][j]<<",";
	}
	fout<<0<<endl;
}

for(j=1;j<=M;j++){
	for(i=1;i<=N-1;i++){
		fout<<u[i][j]<<",";
	}
	fout<<0<<endl;
}

for(j=1;j<=M;j++){
	for(i=1;i<=N-1;i++){
		fout<<v[i][j]<<",";
	}
	fout<<0<<endl;
}

						
	}
///////////////////////////////////////cooling water

//	hwi=wat.sc_hl(100.15,Twi)*1000.0;
//	hwo=hwi+Qsum/Gw;
//	Two=wat.T_Ph(100.15,hwo/1000.0);
//	Rew=Gw/(wat.sc_myul(100.15,Twi)*PI*ri);
//	aw=0.023*pow(wat.sc_Prl(100.15,Two),0.4)*pow(Rew,0.8)*wat.sc_laml(100.15,Twi)/(ri*2.0);
//	K=1/(ri*2.0)/(1/aw/(ri*2.0)+1.0/2.0/400.0*log(ro/ri));	
//	_Tw=(Twi*(exp(K/Qsum*2.0*PI*ri*L*(Twi-Two)))-Two)/((exp(K/Qsum*2.0*PI*ri*L*(Twi-Two)))-1);



	
}*/
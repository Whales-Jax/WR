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
	int i,j,N,M,O,n,p,t;
	N=10;/////Angle's mesh
	M=10;/////Thickness's mesh
	O=10;/////Mass flow rate
	double st, L,rho,myu,beta,hmin,hminp,Re,h[100],g,CA[100],Ti[100][100],Xi[100][100],PI,X[100],f,gl,k,A[100],Reb,dp,b[100],Tw,d,hv,ro,cp,Tsat,m[100],P,Gv,_Tw,Gw,Twi,ri,Xave,Tave,a;
	double dr[100],r[100],R[100],dhde,dtde,dt2de,dxde,dx2de,u,v;

	ofstream fout("GridInv.csv");
	if(!fout){
	cout << "ファイルをオープンできませんでした"<< endl;
    return 1;
	}
	else
		cout<< "ファイルをオープンしました" << endl;


////////////////////Experiment condition////////////////////

///////////////first condition
	Twi=32.0;//water inlet temperature????
	for(i=1;i<=N-1;i++){
		m[i]=0.05;///[kg/s]
						}
	P=1.2;//[kPa]
	Gv=0.1;//???????
///////////////
	t=10;
    L=0.878;//tube length[m]
	ro=9.0/1000.0;
	ri=8.0/1000.0;
	PI=(6*asin(0.5));//π
	beta=PI/double(N);//Angle of tube//////for 180 degree

	Tw=32.0;
	_Tw=30;
	Gw=1.0;///kg/s

	for(i=1;i<=N-1;i++){
		for(j=1;j<=M;j++){
    		Xi[i][j]=0.55;//initialising concentration
			Ti[i][j]=38.0;//initialising temperature
		}
	}
/////////////////////////////////grid generation in the radial direction (?outside loop?)

	for(j=1;j<=M;j++){
	r[j]=1.0/2.0*(1.0-cos((j-1)/(double(M)-1)*PI));
						}
	for(j=1;j<=M-1;j++){
	dr[j]=r[j+1]-r[j];
						}
	for(j=1;j<=M-2;j++){
	R[j]=dr[j+1]/dr[j];
						}
//////////////////////////////////////////////////////Angle loop
	for(i=1;i<=N-1;i++){
		Tave=0;
		Xave=0;
		for(n=1;n<=M;n++){
			Tave+=Ti[i][n]/M;//average temperature for previous angle step
			Xave+=Xi[i][n]/M;//average concentration for previous angle step
							}
//////////////////properties calculation
	st=libr.sc_st_XT(Xave,Tave);//////[N/m]
	rho=libr.sc_rho_XT(Xave,Tave);////[kg/m3]
	myu=libr.sc_visc_XT(Xave,Tave);///[Pas]	
	gl=9.81;//[m/s2]
	k=libr.sc_thc_XT(Xave,Tave);////[W/mK]
	d=libr.sc_d_XT(Xave,Tave);//////[m2/s]
	cp=libr.sc_cp_XT(Xave,Tave)*1000.0;////[J/kgK]

	Tsat=wat.p_t( P );

	dp=(rho*pow(st,3.0))/(pow(myu,4.0)*gl);//dimensionless parameter group for minimum stable thickness calculation
	Reb=2.0*(0.2314*pow(dp,0.1761));//Break-up Reynolds from Maron1982
	hmin=pow(3.0*(0.2314*pow(dp,0.1761))*pow(myu/rho,2.0)/9.81,1.0/3.0);//minimum thickness [m] from Maron1982
	
    Re=2.0*m[i]/L/myu;
	h[i]=pow(3.0*m[i]*myu/(rho*rho*gl*sin(beta*i)),1.0/3.0);//uniform film thickness

		if(i<N/2){////////////first half of tube surface

		if(h[i]>=hmin){
			X[i]=1.0;//complete wetting
						}
		else{
			hminp=pow(((rho*rho*rho)*gl*gl*sin(beta*i)*sin(beta*i))/(15.0*myu*myu*st),0.2)*hmin;//partial wetting

				CA[i]= 1665.4*pow(hminp,6.0) - 3446.8*pow(hminp,5.0) + 2399.3*pow(hminp,4.0) - 459.73*pow(hminp,3.0) + 49.538*pow(hminp,2.0) + 2.9741*hminp + 0.8616;//////////////////////////////////////////////////////
				CA[i]=CA[i]/180*PI;//contact angle from critical condition and minimum stable thickness
					f=-1.0/4.0*pow(cos(CA[i]),3.0)*sin(CA[i]) -13.0/8.0*cos(CA[i])*sin(CA[i])
					-3.0/2.0*CA[i]*pow(sin(CA[i]),2.0) +15.0/8.0*CA[i];//f(θo)
					g=(CA[i]*(5.0/16.0+15.0/4.0*pow(cos(CA[i]),2.0) +5.0/2.0*pow(cos(CA[i]),4.0))
					-sin(CA[i])*(113.0/48.0*cos(CA[i]) +97.0/24.0*pow(cos(CA[i]),3.0) +1.0/6.0*pow(cos(CA[i]),5.0)));//Ψ(θo))

						X[i]=pow((h[i]),3.0)*sin(CA[i])/f*pow(2.0/45.0*(pow(rho,3.0)
						*pow((9.81*sin(beta*i)),2.0)/(st*pow(myu,2.0)))*g/sin(CA[i])*pow((CA[i]/sin(CA[i])-cos(CA[i])),-1),3.0/5.0);//wetting ratio calculation for 10th tube

							A[i]=(log(X[i]+(1-X[i])*pow(Re/Reb,2.0))-log(X[i]))/9.0;//exponential coefficient describing wetting ratio increase for previous tubes
							X[i]=(X[i]+(1-X[i])*pow(Re/Reb,2.0))*exp(-A[i]*(t));//wetting ratio calculation for tube t
							b[i]=(h[i]*(CA[i]/2.0-1.0/2.0*cos(CA[i])*sin(CA[i])))/(pow(X[i]*f,1.0/3.0)*pow(sin(CA[i]),2.0/3.0));//average rivulet thickness
			}
				}
	
	else{////////////////second half of tube surface

		if(h[N/2-1]>=hmin){
			X[i]=1.0;//complete wetting
							}
		else{//dynamic wetting ratio calculation (without minimum stable thickness comparison)
			hminp=pow(((rho*rho*rho)*gl*gl*sin(beta*i)*sin(beta*i))/(15.0*myu*myu*st),0.2)*hmin;

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
									}
						else{
							A[i]=(log(X[i]+(1-X[i])*pow(Re/Reb,2.0))-log(X[i]))/9.0;//exponential coefficient describing wetting ratio increasing for previous tubes
							X[i]=(X[i]+(1-X[i])*pow(Re/Reb,2.0))*exp(-A[i]*(j));
							b[i]=(h[i]*(CA[i]/2.0-1.0/2.0*cos(CA[i])*sin(CA[i])))/(pow(X[i]*f,1.0/3.0)*pow(sin(CA[i]),2.0/3.0));
							}
			}	
		}
/////////
	//m[i+1]=0.051;/////[kg/s]??????
	//Gv=0.0001;//[kg/m2s] ??outside loop??
/////////

	Ti[i][1]=Tw;//wall temperature boundary condition
if (i>=2){
	CNewtonRaphsonMethod mnm;		// CNewtonRaphsonMethodクラスの宣言
	mnm.setup(2000);					//セットアップ　※引数はニュートン法の変数以上にする

	p=0;

	// 初期値の代入
	mnm.setValue( p , Xi[i][M] );
	p++;
	mnm.setValue( p , Ti[i][M] );
	p++;
	for(j=2;j<M;j++){
		mnm.setValue( p , Xi[i][j] );
	p++;
	mnm.setValue( p , Ti[i][j] );
	p++;
	
	}
	mnm.setValue( p , Gv  );
	p++;
	mnm.setValue( p , m[i+1]  );
	p++;

	mnm.setAcc(0.5);				//加速度勾配の入力　※なくても可

	mnm.initial();					//計算を開始する前に必ず初期化してください
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){		// おまじない
		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){	// おまじない
			p=0;
			//値の設定
			Xi[i][M] = mnm.getValue(p);
			p++;
			Ti[i][M] = mnm.getValue(p);
			p++;
	for(j=2;j<M;j++){
			Xi[i][j] = mnm.getValue(p);
			p++;
			Ti[i][j] = mnm.getValue(p);
			p++;
			
						}
	Gv = mnm.getValue(p);
	p++;
	m[i+1] = mnm.getValue(p);
	p++;

	Xi[i][1]=Xi[i][2];/// concentration boundary condition

	hv = wat.sat_hv( Tsat )*1000.0;//-libr.sc_h_XT(Xi[i][M],Ti[i][M])*1000.0;/////[J/kg]?????check calculation formula
	
	p=0;
	if(X[i]==1){
			mnm.setError( p , Gv , -rho*d/(Xi[i][M])*(Xi[i][M]-Xi[i][M-1])/(dr[M-1])/h[i] );//////h or b
			p++;
			}
			else{
			mnm.setError( p , Gv , -rho*d/(Xi[i][M])*(Xi[i][M]-Xi[i][M-1])/(dr[M-1])/b[i] );//////h or b
			p++;			
			}
			mnm.setError( p , Ti[i][M] , libr.sc_T_XTsat(Xi[i][M],wat.p_t( P )) );
            p++;
for(j=2;j<=M-1;j++){
	st=libr.sc_st_XT(Xi[i][j],Ti[i][j]);//////[N/m]
	rho=libr.sc_rho_XT(Xi[i][j],Ti[i][j]);////[kg/m3]
	myu=libr.sc_visc_XT(Xi[i][j],Ti[i][j]);///[Pas]	
	gl=9.81;//[m/s2]
	k=libr.sc_thc_XT(Xi[i][j],Ti[i][j]);////[W/mK]
	d=libr.sc_d_XT(Xi[i][j],Ti[i][j]);//////[m2/s]
	cp=libr.sc_cp_XT(Xi[i][j],Ti[i][j])*1000.0;////[J/kgK]
	a=k/rho/cp;
	u=(rho*gl*h[i]*h[i]*sin(PI*i/double(N))/myu)*(r[j]-1.0/2.0*r[j]*r[j]);
	dhde=-pow(myu*m[i]/L*PI*PI*PI/(9.0*rho*rho*gl),1.0/3.0)*1/(pow(sin(PI*i/double(N)),1.0/3.0)*tan(PI*i/double(N)));
	v=-(rho*gl*h[i]*h[i]*r[j]*r[j])/(2.0*myu*ro)*(1.0/PI*dhde*sin(PI*i/double(N))+h[i]*(1-r[j]/3.0)*cos(PI*i/double(N)));
	dtde=(Ti[i][j+1]-R[j-1]*R[j-1]*Ti[i][j-1]-(1-R[j-1]*R[j-1])*Ti[i][j])/((1+R[j-1])*dr[j]);
	dxde=(Xi[i][j+1]-R[j-1]*R[j-1]*Xi[i][j-1]-(1-R[j-1]*R[j-1])*Xi[i][j])/((1+R[j-1])*dr[j]);
	dt2de=(2.0*R[j-1]*Ti[i][j+1]-2.0*R[j-1]*R[j-1]*Ti[i][j-1]-2.0*R[j-1]*(1+R[j-1])*Ti[i][j])/(dr[j]*dr[j]*(1+R[j-1]));
	dx2de=(2.0*R[j-1]*Xi[i][j+1]-2.0*R[j-1]*R[j-1]*Xi[i][j-1]-2.0*R[j-1]*(1+R[j-1])*Xi[i][j])/(dr[j]*dr[j]*(1+R[j-1]));

	mnm.setError( p , (Xi[i][j]-Xi[i-1][j])/(1/double(N)) , (PI*ro*d)/(u*h[i]*h[i])*dx2de+(r[j]/h[i]*dhde-PI*ro*v/(h[i]*u))*dxde );
	p++;
	mnm.setError( p , (Ti[i][j]-Ti[i-1][j])/(1/double(N)) , (PI*ro*a)/(u*h[i]*h[i])*dt2de+(r[j]/h[i]*dhde-PI*ro*v/(h[i]*u))*dtde );
	p++;
	
	if(j==M-1){
						
			if(X[i]==1){
			mnm.setError( p , Gv*hv , k*(Ti[i][M]-Ti[i][M-1])/(dr[M-1])/h[i] );/////////////////h or b
			p++;
			}
			else{
			mnm.setError( p , Gv*hv , k*(Ti[i][M]-Ti[i][M-1])/(dr[M-1])/b[i] );/////////////////h or b
			p++;
			}
			mnm.setError( p , m[i+1], m[i]+Gv*L*(ro+h[i])*beta );//?????????????? use just 1 eq. ????
			p++;

				}

}


		mnm.prt();				//エラー表示
		mnm.prt_sum();			//エラーの合計を表示

	}
	}
}
else {}
	}
///////////////////////////////////////cooling water

	hwi=wat.sc_hl(100.15,Twi)*1000.0;
	hwo=hwi+Qsum/Gw;
	Two=wat.T_Ph(100.15,hwo/1000.0);
	Rew=Gw/(wat.sc_myul(100.15,Twi)*PI*ri);
	aw=0.023*pow(wat.sc_Prl(100.15,Two),0.4)*pow(Rew,0.8)*wat.sc_laml(100.15,Twi)/(ri*2.0);
	K=1/(ri*2.0)/(1/aw/(ri*2.0)+1.0/2.0/400.0*log(ro/ri));	
	_Tw=(Twi*(exp(K/Qsum*2.0*PI*ri*L*(Twi-Two)))-Two)/((exp(K/Qsum*2.0*PI*ri*L*(Twi-Two)))-1);



}*/
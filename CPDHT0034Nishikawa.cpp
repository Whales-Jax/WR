#include "CPDHT0034Nishikawa.h"


/**********************************************************
核沸騰熱伝達整理式

Ref 日本機械学会論文集42巻361号(昭51-9) 各沸騰熱伝達の整理式における圧力補正項について
calc2:低熱流速
calc:高熱流速

作成 : 中西	(2012/11)
**********************************************************/

CPDHT0034Nishikawa::CPDHT0034Nishikawa(){
}
CPDHT0034Nishikawa::~CPDHT0034Nishikawa(){
}

void CPDHT0034Nishikawa::calc(){
////////つかわないよ
	p=1.699;///kcal/h
	p=p/3600.0/4.184*1000.0;//[J/s]
	M=900.0;//[m-1]
	Pa=101.3;//[kPa]

	// Latent heat [J/kg]
	h_fg = ( fp->hV - fp->h ) * 1000.0;

	temp_alpha = 1000.0;
	while(1){
		// Heat flux [W/m2]
		fp->q = temp_alpha * ( fp->Tw - fp->T );
		f=pow(fp->P/Pa,0.7);////*(1+3*pow(P/P,3.0))/(1+3*pow(Pa/P,3.0));

	X =sqrt(fp->cp*1000.0*fp->D*fp->D/(M*M*p*fp->thc*fp->st*h_fg*fp->DV))*pow(sqrt(sp->D2),3.0)*fp->q;
	// Nusselt number [-]
	fp->Nu1 =0.66*pow(sp->D2,-2.0/5.0)*pow((f*X),4.0/5.0);


	// Heat transfer Coefficient [W/m2K]	
	fp->alpha = fp->Nu1 * fp->thc / sp->D2;
		if( fabs( fp->alpha - temp_alpha ) < 1.0e-8 ){ break; }		
		temp_alpha = fp->alpha;
	}
	// Conversion[W/m2K]->[kW/m2K]
	fp->alpha /= 1000.0;

}


void CPDHT0034Nishikawa::calc2(){

	p=1.699;///kcal/h
	p=p/3600.0/4.184*1000.0;//[J/s]
	M=900.0;//[m-1]
	Pa=101.3;//[kPa]

	// Latent heat [J/kg]
	h_fg = ( fp->hV - fp->h ) * 1000.0;

	temp_alpha = 1000.0;
	while(1){
		// Heat flux [W/m2]
		fp->q = temp_alpha * ( fp->Tw - fp->T );
		f=pow(fp->P/Pa,0.7);////*(1+3*pow(P/P,3.0))/(1+3*pow(Pa/P,3.0));

	X =sqrt(fp->cp*1000.0*fp->D*fp->D/(M*M*p*fp->thc*fp->st*h_fg*fp->DV))*pow(sqrt(sp->D2),3.0)*fp->q;
	// Nusselt number [-]
	fp->Nu1 =6.24*pow((f*X),2.0/3.0);


	// Heat transfer Coefficient [W/m2K]	
	fp->alpha = fp->Nu1 * fp->thc / sp->D2;
		if( fabs( fp->alpha - temp_alpha ) < 1.0e-8 ){ break; }		
		temp_alpha = fp->alpha;
	}
	// Conversion[W/m2K]->[kW/m2K]
	fp->alpha /= 1000.0;

}

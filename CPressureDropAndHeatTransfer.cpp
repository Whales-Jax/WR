#include "CPressureDropAndHeatTransfer.h"

void sannjikannsuu( double x1 , double y1 , double y11 , double x2 , double y2 , double y22 , double &a , double &b , double &c , double &d ){

	/*
	x1  : x1の値
	y1  : x1の時のyの値
	y11 : x1の時のyの一回微分
	x2  : x2の値
	y2  : x2の時のyの値
	y22 : x2の時のyの一回微分
	*/

	a = (-2*y1)/pow(x1 - x2,3) + (x1*y11)/pow(x1 - x2,3) - (x2*y11)/pow(x1 - x2,3) + 
   (2*y2)/pow(x1 - x2,3) + (x1*y22)/pow(x1 - x2,3) - (x2*y22)/pow(x1 - x2,3);

	b = (3*x1*y1)/pow(x1 - x2,3) + (3*x2*y1)/pow(x1 - x2,3) - 
   (pow(x1,2)*y11)/pow(x1 - x2,3) - (x1*x2*y11)/pow(x1 - x2,3) + 
   (2*pow(x2,2)*y11)/pow(x1 - x2,3) - (3*x1*y2)/pow(x1 - x2,3) - 
   (3*x2*y2)/pow(x1 - x2,3) - (2*pow(x1,2)*y22)/pow(x1 - x2,3) + 
   (x1*x2*y22)/pow(x1 - x2,3) + (pow(x2,2)*y22)/pow(x1 - x2,3);

	c = (-6*x1*x2*y1)/pow(x1 - x2,3) + (3*pow(x1,2)*x2*y11)/pow(x1 - x2,3) - 
   (x2*y11)/(x1 - x2) - (3*x1*pow(x2,2)*y11)/pow(x1 - x2,3) + 
   (6*x1*x2*y2)/pow(x1 - x2,3) + (x1*y22)/(x1 - x2) + 
   (3*pow(x1,2)*x2*y22)/pow(x1 - x2,3) - (3*x1*pow(x2,2)*y22)/pow(x1 - x2,3);

	d = (3*x1*pow(x2,2)*y1)/pow(x1 - x2,3) - (pow(x2,3)*y1)/pow(x1 - x2,3) - 
   (pow(x1,2)*pow(x2,2)*y11)/pow(x1 - x2,3) + 
   (x1*pow(x2,3)*y11)/pow(x1 - x2,3) + (pow(x1,3)*y2)/pow(x1 - x2,3) - 
   (3*pow(x1,2)*x2*y2)/pow(x1 - x2,3) - (pow(x1,3)*x2*y22)/pow(x1 - x2,3) + 
   (pow(x1,2)*pow(x2,2)*y22)/pow(x1 - x2,3);

	return;
}


void CPressureDropAndHeatTransfer::init(){
	dbo.fp = &fp;
	dbo.sp = &sp; 

	blas.fp = &fp;
	blas.sp = &sp; 

	fji.fp = &fp;
	fji.sp = &sp;

	gwi.fp = &fp;
	gwi.sp = &sp;

	fjt.fp = &fp;
	fjt.sp = &sp;

	nsl.fp = &fp;
	nsl.sp = &sp;

	msb.fp = &fp;
	msb.sp = &sp;

	isl.fp = &fp;
	isl.sp = &sp;

	grm.fp = &fp;
	grm.sp = &sp;

	kms.fp = &fp;
	kms.sp = &sp;

	sab.fp = &fp;
	sab.sp = &sp;

	chn.fp = &fp;
	chn.sp = &sp;

	knd.fp = &fp;
	knd.sp = &sp;

	sha.fp = &fp;
	sha.sp = &sp;

	jon.fp = &fp;
	jon.sp = &sp;

	liu.fp = &fp;
	liu.sp = &sp;

	got.fp = &fp;
	got.sp = &sp;

	Grv.fp = &fp;
	Grv.sp = &sp;

	Thr.fp = &fp;
	Thr.sp = &sp;

	Nis.fp = &fp;
	Nis.sp = &sp;


}
void CPressureDropAndHeatTransfer::setup(){

	dbo.ref.CopyTable( ref );

	fji.ref.CopyTable( ref );

	noz.ref.CopyTable( ref );

	gwi.ref.CopyTable( ref );

	blas.ref.CopyTable( ref );

	ito.ref.CopyTable( ref );

	ysh.ref.CopyTable( ref );

	chi.ref.CopyTable( ref );

	car.ref.CopyTable( ref );

	Yon.ref.CopyTable( ref );

	Yoh.ref.CopyTable( ref );

	Kub.ref.CopyTable( ref );

	Mor.ref.CopyTable( ref );

	got.ref.CopyTable( ref );

	Grv.ref.CopyTable( ref );
	Grv.setup();

	IST.ref.CopyTable( ref );
	IST.setup();

	gni.ref.CopyTable( ref );


}

void CPressureDropAndHeatTransfer::setup_water(double P, double h ){

	fp.P	 = P;
	fp.h	 = h;
	fp.satT  = wat2.p_t( fp.P );
	fp.hL	 = wat2.sat_hl( fp.satT );
	fp.hV	 = wat2.sat_hv( fp.satT );

	fp.DL	 = wat2.sat_roul( fp.satT );
	fp.DV	 = wat2.sat_rouv( fp.satT );
	fp.cpL	 = wat2.sat_cpl ( fp.satT );
	fp.cpV	 = wat2.sat_cpv ( fp.satT );
	fp.viscL = wat2.sat_myul( fp.satT );
	fp.viscV = wat2.sat_myuv( fp.satT );
	fp.thcL	 = wat2.sat_laml( fp.satT );
	fp.thcV	 = wat2.sat_lamv( fp.satT );
	fp.PrL	 =  fp.viscL * fp.cpL * 1000.0 / fp.thcL;
	fp.PrV	 =  fp.viscV * fp.cpV * 1000.0 / fp.thcV;

	fp.x	 = ( fp.h - fp.hL ) / ( fp.hV - fp.hL );

	if( h <= fp.hL ){
		fp.T	 = wat2.sc_tl  ( fp.P, fp.h );
		fp.v	 = wat2.sc_vl  ( fp.P, fp.T );
		fp.D     = 1.0 / fp.v;
		fp.cp	 = wat2.sc_cpl ( fp.P, fp.T );
		fp.visc  = wat2.sc_myul( fp.P, fp.T );
		fp.thc   = wat2.sc_laml( fp.P, fp.T );
		fp.Pr1	 = fp.visc * fp.cp *1000.0 / fp.thc;
	}else if( h < fp.hV ){
		fp.T	 = fp.satT;
		fp.v	 = ( 1.0 - fp.x ) * ( 1.0 / fp.DL )  + fp.x * ( 1.0 / fp.DV );
		fp.D     = 1.0 / fp.v;
	}
	else{
		fp.T	 = wat2.sh_tv  ( fp.P, fp.h );
		fp.v	 = wat2.sh_vv  ( fp.P, fp.T );
		fp.D     = 1.0 / fp.v;
		fp.cp	 = wat2.sh_cpv ( fp.P, fp.T );
		fp.visc  = wat2.sh_myuv( fp.P, fp.T );
		fp.thc   = wat2.sh_lamv( fp.P, fp.T );
		fp.Pr1	 = fp.visc * fp.cp * 1000.0 / fp.thc;
	}
}

//作動流体が臭化リチウムの場合
void CPressureDropAndHeatTransfer::setup_libr(double P, double h, double X ){

	double Le = 100.0;

	fp.P	 = P;
	fp.h	 = h;
	fp.X	 = X;
	fp.T     = libr.sc_T_Xh( fp.X, fp.h );
	fp.D     = libr.sc_rho_XT( fp.X, fp.T );
	fp.visc  = libr.sc_visc_XT( fp.X, fp.T );//[Pa*s]
	fp.thc   = libr.sc_thc_XT( fp.X, fp.T );
	fp.cp    = libr.sc_cp_XT( fp.X, fp.T );
	fp.st    = libr.sc_st_XT( fp.X, fp.T );//
	fp.d     = libr.sc_d_XT( fp.X, fp.T );//
	fp.nu    = fp.visc / fp.D;
	fp.Pr1	 = fp.visc * fp.cp * 1000.0 / fp.thc;
	
	fp.satT  = wat2.p_t( fp.P );
	fp.hV    = wat2.sh_hv( fp.P, fp.T );
	fp.DV    = wat2.sh_rouv( fp.P, fp.T );
}


void CPressureDropAndHeatTransfer::calc_ht2(){
	
	if( strcmp( fp.HtransferCode, "Nusselt" ) == 0 )
		nsl.calc2();
	else if( strcmp( fp.HtransferCode, "Fujita" ) == 0 )
		fjt.calc2();
	else if( strcmp( fp.HtransferCode, "Dittus Boelter" ) == 0 )
		dbo.calc2();
	else if( strcmp( fp.HtransferCode, "Gungor & Winterton" ) == 0 )
		gwi.calc2();
	else if( strcmp( fp.HtransferCode, "Islam" ) == 0 )
		isl.calc2();
	else if( strcmp( fp.HtransferCode, "Garimella" ) == 0 )
		grm.calc2();
	else if( strcmp( fp.HtransferCode, "Garimella_form" ) == 0 )
		grm.calc_form();
	else if( strcmp( fp.HtransferCode, "Kamoshida" ) == 0 )
		kms.calc2();
	else if( strcmp( fp.HtransferCode, "Kamoshida2003" ) == 0 )
		kms.calc();
	else if( strcmp( fp.HtransferCode, "Kamoshida2004" ) == 0 )
		kms.calc3();
	else if( strcmp( fp.HtransferCode, "Stephan & Abdelsalam_LiBr" ) == 0 )
		sab.calc2();
	else if( strcmp( fp.HtransferCode, "Stephan & Abdelsalam_Water" ) == 0 )
		sab.calc();
	else if( strcmp( fp.HtransferCode, "Chen" ) == 0 )
		chn.calc2();
	else if( strcmp( fp.HtransferCode, "Kandlikar" ) == 0 )
		knd.calc2();
	else if( strcmp( fp.HtransferCode, "Shah" ) == 0 )
		sha.calc2();
	else if( strcmp( fp.HtransferCode, "JongTaekOh" ) == 0 )
		sha.calc2();
	else if( strcmp( fp.HtransferCode, "LiuWinterton" ) == 0 )
		sha.calc2();
	else if( strcmp( fp.HtransferCode, "Thrope" ) == 0 )
		Thr.calc();
	else if( strcmp( fp.HtransferCode, "Nishikawa" ) == 0 )
		Nis.calc();
	else if( strcmp( fp.HtransferCode, "Nishikawa2" ) == 0 )
		Nis.calc2();
	else if( strcmp( fp.HtransferCode, "Meisenburg" ) == 0 )
		msb.calc2();

	else
		return;
}

/*
//熱伝達率の計算
void CPressureDropAndHeatTransfer::calc_ht(){

	if( G < 0.000001 ){
		G = 0.000001;
	}

	ref.state_ph2( P , h );
	x = ref.Rc.x;
	if( ref.Rc.T > T_wi ){		//もしも冷媒の温度が壁温よりも高ければ，
		kanetu = 0;				//hantei == 0　流体を加熱した場合(?)
	}else{						//冷媒の温度が壁温よりも高くなければ，
		kanetu = 1;				//hantei == 1 流体を冷却した場合(?)
	}




	if( 0.02 < x && x < 0.98 ){				//二相流 乾き度が0.02～0.98の場合，式をそのまま使用する．

		if( kanetu == 0 ){					//もしも冷媒の温度が壁温よりも高ければ，凝縮器として働く．管内凝縮熱伝達率の計算に藤井らの式を用いる．
			noz.G = G;
			noz.P = P;
			noz.h = h;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){			//もしも冷媒の温度が壁温よりも高くなければ，壁面からの熱で蒸発させる．蒸発器として働く．二相蒸発の場合は吉田らの式を用いる．
			ysh.G = G;
			ysh.P = P;
			ysh.h = h;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

	}else if( -0.02 < x && x <= 0.02 ){		//乾き度が-0.02～0.02の区間では三次関数による補間を行う．

		double x1 = -0.02;
		double x2 =  0.02;
		double a,b,c,d;
		double y1,y11,y2,y22;
		double h1;
		double h2;
		double delta = 0.0001;

		ref.sat_p2( P );

		h1 = ref.Rl.h - ( ref.Rv.h - ref.Rl.h ) * 0.02;
		h2 = ref.Rl.h + ( ref.Rv.h - ref.Rl.h ) * 0.02;

		//h1過熱判定
		ref.state_ph2( P , h1 );
		if( ref.Rc.T > T_wi ){
			kanetu = 0;
		}else{
			kanetu = 1;
		}

		dbo.G = G;
		dbo.P = P;
		dbo.h = h1;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;

		y1 = alpha;

		dbo.G = G;
		dbo.P = P;
		dbo.h = h1-h1*delta;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;

		y11 = ( y1-alpha ) / ( h1*delta );

		//h2過熱判定
		ref.state_ph2( P , h2 );
		if( ref.Rc.T > T_wi ){
			kanetu = 0;
		}else{
			kanetu = 1;
		}
		if( kanetu == 0 ){				//凝縮
			noz.G = G;
			noz.P = P;
			noz.h = h2;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){		//蒸発
			ysh.G = G;
			ysh.P = P;
			ysh.h = h2;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

		y2 = alpha;



		if( kanetu == 0 ){//凝縮
			noz.G = G;
			noz.P = P;
			noz.h = h2+h2*delta;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){//蒸発
			ysh.G = G;
			ysh.P = P;
			ysh.h = h2+h2*delta;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

		y22 = ( alpha-y2 ) / ( h2*delta );


		sannjikannsuu( x1 , y1 , y11 , x2 , y2 , y22 , a , b , c , d );

		alpha = a * pow( x , 3 ) + b * pow( x , 2 ) + c * x + d;

	}else if( 0.98 <= x && x < 1.02 ){			//乾き度が0.98～1.02の区間では三次関数による補間を行う．


		double x1 = 0.98;
		double x2 = 1.02;
		double a,b,c,d;
		double y1,y11,y2,y22;
		double h1;
		double h2;
		double delta = 0.0001;

		ref.sat_p2( P );

		h1 = ref.Rv.h - ( ref.Rv.h - ref.Rl.h ) * 0.02;
		h2 = ref.Rv.h + ( ref.Rv.h - ref.Rl.h ) * 0.02 ;



		//h1過熱判定
		ref.state_ph2( P , h1 );
		if( ref.Rc.T > T_wi ){
			kanetu = 0;
		}else{
			kanetu = 1;
		}
		if( kanetu == 0 ){//凝縮
			noz.G = G;
			noz.P = P;
			noz.h = h1;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){//蒸発
			ysh.G = G;
			ysh.P = P;
			ysh.h = h1;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

		y1 = alpha;



		if( kanetu == 0 ){//凝縮
			noz.G = G;
			noz.P = P;
			noz.h = h1-h1*delta;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){//蒸発
			ysh.G = G;
			ysh.P = P;
			ysh.h = h1-h1*delta;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

		y11 = ( alpha-y1 ) / ( h1*delta );


		//h2過熱判定
		ref.state_ph2( P , h2 );
		if( ref.Rc.T > T_wi ){
			kanetu = 0;
		}else{
			kanetu = 1;
		}
		dbo.G = G;
		dbo.P = P;
		dbo.h = h2;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;

		y2 = alpha;


		dbo.G = G;
		dbo.P = P;
		dbo.h = h2+h2*delta;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;

		y22 = ( y2-alpha ) / ( h2*delta );



		sannjikannsuu( x1 , y1 , y11 , x2 , y2 , y22 , a , b , c , d );

		alpha = a * pow( x , 3 ) + b * pow( x , 2 ) + c * x + d;

		cout << "";




	}else{										//単相流の場合，管内熱伝達率の計算にDittus-Boelterの式を用いる．
		dbo.G = G;
		dbo.P = P;
		dbo.h = h;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;
	}



}

*/


/*
void CPressureDropAndHeatTransfer::calc_pd(){	//冷媒側の圧力損失の計算

	if( G < 0.000001 ){
		G = 0.000001;
	}

	ref.state_ph2( P , h );
	x = ref.Rc.x;

	ref.sat_p2( P );

	if( ref.Rc.T > T_wi ){
		kanetu = 0;
	}else{
		kanetu = 1;
	}

	if( 0.0 < x && x < 1.0 ){		//二相流の場合は，圧力損失の計算にはChisolmの式と均質流の式を用い，それぞれのdpdzで小さいほうをdpdzとして採用する．
		chi.G = G;
		chi.P = P;
		chi.h = h;
		chi.D = D;
		chi.kanetu = kanetu;
		chi.calc();


		blas.G = G;
		blas.P = P;
		blas.h = h;
		blas.D = D;
		blas.kanetu = kanetu;
		blas.calc();


//		if( chi.dpdz < blas.dpdz ){
//			dpdz = blas.dpdz;
//		}else{
			dpdz = chi.dpdz;
//		}

	}else{							//単相流の場合，均質流の式で計算する．
		blas.G = G;
		blas.P = P;
		blas.h = h;
		blas.D = D;
		blas.kanetu = kanetu;
		blas.calc();
		dpdz = blas.dpdz;
	}


}*/

//溝付き管の熱伝達率の計算
/*
void Grooved_tube::calc_ht(){

	this->noz.ref.data_table   = ref.data_table;
	this->ysh.ref.data_table   = ref.data_table;
	this->dbo.ref.data_table   = ref.data_table;


	if( G < 0.000001 ){
		G = 0.000001;
	}

	ref.state_ph2( P , h );
	x = ref.Rc.x;
	if( ref.Rc.T > T_wi ){		//もしも冷媒の温度が壁温よりも高ければ，
		kanetu = 0;				//hantei == 0　流体を加熱した場合(?)
	}else{						//冷媒の温度が壁温よりも高くなければ，
		kanetu = 1;				//hantei == 1 流体を冷却した場合(?)
	}



	if( 0.02 < x && x < 0.98 ){				//二相流 乾き度が0.02～0.98の場合，式をそのまま使用する．

		if( kanetu == 0 ){					//もしも冷媒の温度が壁温よりも高ければ，凝縮器として働く．管内凝縮熱伝達率の計算に藤井らの式を用いる．
			noz.G = G;
			noz.P = P;
			noz.h = h;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){			//もしも冷媒の温度が壁温よりも高くなければ，壁面からの熱で蒸発させる．蒸発器として働く．二相蒸発の場合は吉田らの式を用いる．
			ysh.G = G;
			ysh.P = P;
			ysh.h = h;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

	}else if( -0.02 < x && x <= 0.02 ){		//乾き度が-0.02～0.02の区間では三次関数による補間を行う．

		double x1 = -0.02;
		double x2 =  0.02;
		double a,b,c,d;
		double y1,y11,y2,y22;
		double h1;
		double h2;
		double delta = 0.0001;

		ref.sat_p2( P );

		h1 = ref.Rl.h - ( ref.Rv.h - ref.Rl.h ) * 0.02;
		h2 = ref.Rl.h + ( ref.Rv.h - ref.Rl.h ) * 0.02;

		//h1過熱判定
		ref.state_ph2( P , h1 );
		if( ref.Rc.T > T_wi ){
			kanetu = 0;
		}else{
			kanetu = 1;
		}

		dbo.G = G;
		dbo.P = P;
		dbo.h = h1;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;

		y1 = alpha;

		dbo.G = G;
		dbo.P = P;
		dbo.h = h1-h1*delta;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;

		y11 = ( y1-alpha ) / ( h1*delta );

		//h2過熱判定
		ref.state_ph2( P , h2 );
		if( ref.Rc.T > T_wi ){
			kanetu = 0;
		}else{
			kanetu = 1;
		}
		if( kanetu == 0 ){				//凝縮
			noz.G = G;
			noz.P = P;
			noz.h = h2;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){		//蒸発
			ysh.G = G;
			ysh.P = P;
			ysh.h = h2;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

		y2 = alpha;



		if( kanetu == 0 ){//凝縮
			noz.G = G;
			noz.P = P;
			noz.h = h2+h2*delta;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){//蒸発
			ysh.G = G;
			ysh.P = P;
			ysh.h = h2+h2*delta;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

		y22 = ( alpha-y2 ) / ( h2*delta );


		sannjikannsuu( x1 , y1 , y11 , x2 , y2 , y22 , a , b , c , d );

		alpha = a * pow( x , 3 ) + b * pow( x , 2 ) + c * x + d;

	}else if( 0.98 <= x && x < 1.02 ){			//乾き度が0.98～1.02の区間では三次関数による補間を行う．


		double x1 = 0.98;
		double x2 = 1.02;
		double a,b,c,d;
		double y1,y11,y2,y22;
		double h1;
		double h2;
		double delta = 0.0001;

		ref.sat_p2( P );

		h1 = ref.Rv.h - ( ref.Rv.h - ref.Rl.h ) * 0.02;
		h2 = ref.Rv.h + ( ref.Rv.h - ref.Rl.h ) * 0.02 ;



		//h1過熱判定
		ref.state_ph2( P , h1 );
		if( ref.Rc.T > T_wi ){
			kanetu = 0;
		}else{
			kanetu = 1;
		}
		if( kanetu == 0 ){//凝縮
			noz.G = G;
			noz.P = P;
			noz.h = h1;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){//蒸発
			ysh.G = G;
			ysh.P = P;
			ysh.h = h1;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

		y1 = alpha;



		if( kanetu == 0 ){//凝縮
			noz.G = G;
			noz.P = P;
			noz.h = h1-h1*delta;
			noz.D = D;
			noz.T_wi = T_wi;
			noz.kanetu = kanetu;
			noz.calc();
			alpha = noz.alpha;
		}else if( kanetu == 1 ){//蒸発
			ysh.G = G;
			ysh.P = P;
			ysh.h = h1-h1*delta;
			ysh.D = D;
			ysh.T_wi = T_wi;
			ysh.kanetu = kanetu;
			ysh.calc();
			alpha = ysh.alpha;
		}

		y11 = ( alpha-y1 ) / ( h1*delta );


		//h2過熱判定
		ref.state_ph2( P , h2 );
		if( ref.Rc.T > T_wi ){
			kanetu = 0;
		}else{
			kanetu = 1;
		}
		dbo.G = G;
		dbo.P = P;
		dbo.h = h2;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;

		y2 = alpha;


		dbo.G = G;
		dbo.P = P;
		dbo.h = h2+h2*delta;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;

		y22 = ( y2-alpha ) / ( h2*delta );



		sannjikannsuu( x1 , y1 , y11 , x2 , y2 , y22 , a , b , c , d );

		alpha = a * pow( x , 3 ) + b * pow( x , 2 ) + c * x + d;

		cout << "";




	}else{										//単相流の場合，管内熱伝達率の計算にDittus-Boelterの式を用いる．
		dbo.G = G;
		dbo.P = P;
		dbo.h = h;
		dbo.D = D;
		dbo.kanetu = kanetu;
		dbo.calc();
		alpha = dbo.alpha;
	}



}
void Grooved_tube::calc_pd(){	//冷媒側の圧力損失の計算

	if( G < 0.000001 ){
		G = 0.000001;
	}

	ref.state_ph2( P , h );
	x = ref.Rc.x;

	ref.sat_p2( P );

	if( ref.Rc.T > T_wi ){
		kanetu = 0;
	}else{
		kanetu = 1;
	}

	if( 0.0 < x && x < 1.0 ){		//二相流の場合は，圧力損失の計算にはChisolmの式と均質流の式を用い，それぞれのdpdzで小さいほうをdpdzとして採用する．
		chi.G = G;
		chi.P = P;
		chi.h = h;
		chi.D = D;
		chi.kanetu = kanetu;
		chi.calc();

		dpdz = chi.dpdz;

	}else{							//単相流の場合，均質流の式で計算する．
		blas.G = G;
		blas.P = P;
		blas.h = h;
		blas.D = D;
		blas.kanetu = kanetu;
		blas.calc();
		dpdz = blas.dpdz;
	}


}
*/

/*
void yoshida::calc(){//吉田らの式
	ref.state_ph2( P , h );
	ref.sat_p2( P );

	lambda_l = ref.Rl.thc;

	T = ref.Rc.T;

	//流路断面積の計算S(m^2)
	S = D * D * atan ( 1.0 );

	//冷媒質量速度G2(kg/(m^2・s))
	G2 = G / S;

	rho_l = ref.Rl.rho;			//液相での密度の導出
	rho_v = ref.Rv.rho;			//気相での密度の導出

	v_l = 1.0 / rho_l;			//液相での比体積の計算
	v_v = 1.0 / rho_v;			//気相での比体積の計算
	
	//粘度viscosityの計算(単位の変換 viscの単位がはマイクロなので)
	mu_l = ref.Rl.visc / 1000000;
	mu_v = ref.Rv.visc / 1000000;

	//乾き度
	x = ref.Rc.x;
	Xtt = pow( mu_l / mu_v , 0.1 ) * pow( v_l / v_v , 0.5 ) * pow( (1-x)/x , 0.9 ); 

	hv = ( ref.Rv.h - ref.Rl.h ) * 1000.0;

	Cp_l = ref.Rl.cp * 1000.0;

	Pr_l = Cp_l * mu_l / lambda_l;
	Re_l = G2 * ( 1 - x ) * D / mu_l;

	alpha_l = 0.023 * pow( Re_l , 0.8 ) * pow( Pr_l , 0.4 ) * lambda_l / D;

	Bo = 0.002;//過程
	_Bo = 0.002;//過程

	mnm.setValue( 0, Bo );
	mnm.initial();
	mnm.acc = 0.1;
	mnm.setSolver(1);

	
	for(mnm.main_loop_init();mnm.main_loop_check();mnm.main_loop_reinit()){

		mnm.acc *= 1.1;
		if( mnm.acc > 1.0 ){
			mnm.acc = 1.0;
		}

		for(mnm.sub_loop_init();mnm.sub_loop_check();mnm.sub_loop_reinit()){
			Bo = mnm.getValue(0);
			E_new = 3.7 * pow( Bo * 10000 + 0.23 * pow( Bo * 10000 , 0.67 ) * pow( 1.0 / Xtt , 2.0 ) , 0.44 );
			alpha = alpha_l * E_new;
			_Bo = alpha * ( T_wi - T ) / G2 / (hv);

			mnm.setError( 0 , Bo , _Bo );

		}
	}


//	int i = 0;
//	while( i < 1000 ){

//		E_new = 3.7 * pow( Bo * 10000 + 0.23 * pow( Bo * 10000 , 0.67 ) * pow( 1.0 / Xtt , 2.0 ) , 0.44 );
//		alpha = alpha_l * E_new;
//		_Bo = alpha * ( T_wi - T ) / G2 / (hv);

//		if( fabs( _Bo - Bo ) < 1.0e-8 ){
//			break;
//		}
//		Bo = _Bo;
//		i++;
//	}

	alpha /= 1000.0;

}
*/
/*
void chisholm::calc(){		//Chisolmの式　圧力損失に用いられる．

	ref.state_ph2( P , h );
	ref.sat_p2( P );
	//断面積の計算S(m^2)
	S = D * D * atan ( 1.0 ) ;

	//冷媒質量速度G2の計算G2(kg/(m^2・s))
	G2 = G / S;

	rho_l = ref.Rl.rho;
	rho_v = ref.Rv.rho;

	v_l = 1.0 / rho_l;
	v_v = 1.0 / rho_v;
	
	//粘度
	mu_l = ref.Rl.visc / 1000000;
	mu_v = ref.Rv.visc / 1000000;

	//乾き度
	x = ref.Rc.x;
	Xtt = pow( mu_l / mu_v , 0.1 ) * pow( v_l / v_v , 0.5 ) * pow( (1-x)/x , 0.9 ); 

	phai2_l = 1.0 + 20.0 / Xtt + 1.0 / ( Xtt * Xtt);
	phai2_v = 1.0 + 20.0 * Xtt + 1.0 * ( Xtt * Xtt);

	S = D * D * atan ( 1.0 ) ;

	G2_l = (1-x) * G / S;
	G2_v = x * G / S;

	Re_l = G2_l * D / mu_l;
	Re_v = G2_v * D / mu_v;

//	if( Re_l > 2000 ){
		f_l = 0.079 * pow( Re_l , -0.25 );
//	}else{
//		f_l = 16.0 / Re_l;
//	}
//	if( Re_v > 2000 ){
		f_v = 0.079 * pow( Re_v , -0.25 );
//	}else{
//		f_v = 16.0 / Re_v;
//	}

	dpdz_l = 2.0 * f_l * pow( G2_l / rho_l , 2.0 ) / D * rho_l;
	dpdz_v = 2.0 * f_v * pow( G2_v / rho_v , 2.0 ) / D * rho_v;

	dpdz = dpdz_l * phai2_l;
//	dpdz = dpdz_v * phai2_v;

	dpdz /= 1000.0;

//	if( 0.95 < x && x <= 1.0 ){
//		dpdz *= -sin( ( x - 1 ) / 0.1 * 4.0 * atan ( 1.0 ) );
//	}
//	if( 0.0 < x && x <= 0.05 ){
//		dpdz *= sin( ( x ) / 0.1 * 4.0 * atan ( 1.0 ) );
//	}
}
*/
/*
void itou_pd::calc(){		//伊藤の式　曲がり配管に使われている． 圧力損失の式

	ref.state_ph2( P , h );
	ref.sat_p2( P );

	x = ref.Rc.x;
	rho = ref.Rc.rho;


	mu = ref.Rc.visc / 1000000;

	mu_v = ref.Rv.visc / 1000000;
	mu_l = ref.Rl.visc / 1000000;

	if( 0.0 < x && x < 1.0 ){
		mu = 1.0 / ( x / mu_v + (1-x) / mu_l );
//		mu = x * mu_v + (1-x) * mu_l;
	}

	double pi;
	pi = atan( 1.0 ) * 4.0;
	S = D * D * pi / 4.0;

	G2 = G / S;

	Re = G2 * D / mu;

	double temp1;

	temp1 = Re * pow( D/R , 2.0 );//摩擦損失の式(伊藤の式)Re(D/R)

	alpha = 1.0 + 5.06 * pow( R/D , -4.52 );

//	double theta_rad;
//	theta_rad = theta/360*pi;
//	theta_rad = theta;

	if( temp1 <=1.5){

		ramda = 0.00873 * alpha * (0.029+0.304*pow( Re*pow(D/2/R,2) ,-0.25))*pow(2*R/D,0.5)*theta;//伊藤の乱れた曲流れの管摩擦係数を用いた実験公式 管摩擦係数は "速研報告　第14巻146頁式(2)"より，摩擦損失の係数を求める式"速研報告　第15巻16頁 式(9)"より ただし曲管摩擦係数の適用範囲は300>Re*(a/R)^2>0.034

	}else if( (1.5 < temp1)&&(temp1 <= 364 )){

		ramda = 0.00515 * alpha * theta * pow( Re , -0.2 ) * pow( R/D , 0.9 );//伊藤の乱れた曲流れの管摩擦係数を用いた実験公式 管摩擦係数は "速研報告　第14巻146頁式(3)"より，摩擦損失の係数を求める式"速研報告　第15巻16頁 式(9)"より　ただし曲管摩擦係数の適用範囲はRe(a/R)^2>6で摩擦損失の係数の式の適用範囲はRe*(a/R)^2<91

	}else{
		ramda = 0.00431 * alpha * theta * pow( Re , -0.17 ) * pow( R/D , 0.84 );//伊藤の摩擦損失の係数を求める式"速研報告　第15巻16頁 式(10)"より ただし式の適用範囲はRe*(a/R)^2>91
	}

	dpdz = ramda / 2.0 * pow( G2 / rho, 2.0 ) * rho;
	dpdz /= 1000.0;

}*/

/*
void blassitsuryu::calc(){		//均質流の式 圧力損失の式

	ref.state_ph2( P , h );
	ref.sat_p2( P );

	x = ref.Rc.x;
	rho = ref.Rc.rho;


	mu = ref.Rc.visc / 1000000;//visc:粘度の計算(viscの単位がマイクロなので10^6で割っている)

	mu_v = ref.Rv.visc / 1000000;
	mu_l = ref.Rl.visc / 1000000;

	if( 0.0 < x && x < 1.0 ){					//湿り蒸気の場合の粘度の計算?　二相流と考え，Macadamusの式をもちいる．
		mu = 1.0 / ( x / mu_v + (1-x) / mu_l );	//macadamus:二相流の式の定義
//		mu = x * mu_v + (1-x) * mu_l;//cicchitti
	}

	//断面積の計算
	S = D * D * atan ( 1.0 ) ;
	//冷媒流速の計算
	G2 = G / S;
	//レイノルズ数の計算
	Re = G2 * D / mu;

//	if( Re > 2000 ){
		f = 0.079 * pow( Re , -0.25 );//ブラシウスの式(単相流？)
//	}else{
//		f = 16.0 / Re;
//	}

	dpdz = 2.0 * f * pow( G2 / rho, 2.0 ) / D * rho;//冷媒側圧力損失の式
	dpdz /= 1000.0;

}

*/

/**********************************************************
carnavosの式
内面溝付管の単相流の圧損の式

2012.06　渡辺寿朗作成

Aef : 実効流路断面積 [m2]
D   : 最大内径 [m]
Dn  : 公称内径 (4*PI*Aef)^0.5 [m]
Dh  : 水力相当直径 4*Aef/(PI*Dn) [m]
A   : 最大内径の時の流路断面積 (D/2)^2*PI [m2]
beta: 溝のリード角[deg]

**********************************************************/
/*
void blassitsuryu2::calc(){		//均質流の式 圧力損失の式

	ref.state_ph2( P , h );
	ref.sat_p2( P );

	x = ref.Rc.x;
	rho = ref.Rc.rho;

	mu = ref.Rc.visc / 1000000;//visc:粘度の計算(viscの単位がマイクロなので10^6で割っている)

	mu_v = ref.Rv.visc / 1000000;
	mu_l = ref.Rl.visc / 1000000;

	if( 0.0 < x && x < 1.0 ){					//湿り蒸気の場合の粘度の計算?　二相流と考え，Macadamusの式をもちいる．
		mu = 1.0 / ( x / mu_v + (1-x) / mu_l );	//macadamus:二相流の式の定義

	}

	//公称内径の計算
	Dn = pow ( 4 * 4 * atan (1.0) * Aef , 0.5 );

	//冷媒流速の計算
	G2 = G / Aef;

	//レイノルズ数の計算
	Re = G2 * Dn / mu;

	//水力相当直径の計算
	Dh = 4 * Aef / ( 4 * atan ( 1.0 ) * Dn );

	//fの計算
	f = 0.046 * pow( Re , -0.20 ) * ( Da / Dh ) * pow ( 1.0 / cos( beta ) , 0.75 ) ;//carnavosの式(単相流？)

	dpdz = 2.0 * f * pow( G2 / rho, 2.0 ) / Da * rho;//冷媒側圧力損失の式
	dpdz /= 1000.0;

}*/

/*
void gungor_winterton::calc(){//Gungor_Wintertonの式　alphaの計算

	ref.state_ph2( P , h );
	ref.sat_p2( P );

	lambda_l = ref.Rl.thc;

	T = ref.Rc.T;

	//流路断面積
	S = D * D * atan ( 1.0 );

	//冷媒質量速度
	G2 = G / S;

	rho_l = ref.Rl.rho;
	rho_v = ref.Rv.rho;

	//粘土
	mu_l = ref.Rl.visc / 1000000;
	mu_v = ref.Rv.visc / 1000000;

	//乾き度
	x = ref.Rc.x;

	hv = ( ref.Rv.h - ref.Rl.h ) * 1000.0;//蒸発潜熱(?)の計算 および,比エンタルピーの単位を(kJ/kg)から(J/kg)へ

	Cp_l = ref.Rl.cp * 1000.0;//液相での定積比熱(kJ/kg)を(J/kg)へ

	Pr_l = Cp_l * mu_l / lambda_l;//液相でのプラントル数の計算
	Re_l = G2 * ( 1 - x ) * D / mu_l;//液相でのレイノルズ数の計算

	alpha_l = 0.023 * pow( Re_l , 0.8 ) * pow( Pr_l , 0.4 ) * lambda_l / D;//単相流で管内が加熱される場合(蒸発?吉田らの式を適用?)のNu数からalphaの計算

	Bo = 0.0002;//過程
	_Bo = 0.0002;//過程

	while(1){

		E_new = 1.0 + 3000.0 * pow( Bo , 0.86 ) + 1.12 * pow( (x)/(1-x) , 0.75 ) * pow( rho_l / rho_v , 0.41 );
		alpha = alpha_l * E_new;
		_Bo = alpha * ( T_wi - T ) / G2 / (hv);

		if( fabs( _Bo - Bo ) < 1.0e-6 ){
			break;
		}
		Bo = _Bo;
	}

	alpha /= 1000.0;

}

/*
void fujii::calc(){

	ref.state_ph2( P , h );
	ref.sat_p2( P );

	lambda_l = ref.Rl.thc;

	//流路断面積
	S = D * D * atan ( 1.0 );

	//冷媒質量速度
	G2 = G / S;

	//比体積
	v_l = 1.0 / ref.Rl.rho;
	v_g = 1.0 / ref.Rv.rho;

	rho_l = ref.Rl.rho;
	rho_v = ref.Rv.rho;

	//粘度
	mu_l = ref.Rl.visc / 1000000;
	mu_v = ref.Rv.visc / 1000000;

	//乾き度
	x = ref.Rc.x;

	T_wi += 273.15;
	T_sat = ref.Rv.T + 273.15;
	Cp_l = ref.Rl.cp * 1000.0;


	//レイノルズ数
	Re_l = G2 * ( 1 - x ) * D / mu_l;
	Re = G2 * D / mu_l;

	Pr_l = Cp_l * mu_l / lambda_l;

	g = 9.8;

	Ga = g * pow( rho_l , 2.0 ) * pow( D , 3.0 ) / pow( mu_l , 2 );

	hv = ( ref.Rv.h - ref.Rl.h ) * 1000.0;

	H_l = Cp_l * ( T_sat - T_wi ) / hv; 

	double temp;
	temp = ( rho_l / rho_v + 0.4 * (1.0-x)/x ) / (1.0 + 0.4 * (1-x)/x );

	gzai = 1.0 + rho_v / rho_l * (1.0-x)/x * ( 0.4 + 0.6 * sqrt( temp ) );
	gzai = 1.0 / gzai;

	A = 10.0 * pow( 1.0 - gzai , 0.1 ) - 8.9;

//	H = gzai + A * sqrt( gzai ) * ( 1.0 - sqrt( gzai ) );

	H = gzai + ( 10.0 * ( pow( 1.0 - gzai , 0.1) -1.0 ) + 1.7 * 1e-4 * Re ) * sqrt( gzai ) * ( 1.0 - sqrt( gzai ) );

	Nu_b = 0.725 * H * pow( Ga * Pr_l / H_l , 0.25 );


	X_tt = pow( (1.0-x)/x  , 0.9) * pow( rho_v / rho_l , 0.5 ) * pow( mu_l / mu_v , 0.5 );

	fai_v = 1.0 + 0.5 * pow( G2 / sqrt( g * D * rho_v * ( rho_l - rho_v ) ) , 0.75 ) *
		pow( X_tt , 0.35 );
	
	taw_wv = 0.023 * pow( G2 , 2.0 ) * pow( x , 2.0 ) * rho_v / 
		pow( G2 * x * D / mu_v , 0.2 );

	taw_w = taw_wv * pow( fai_v , 2.0 );

	Re_lstar = rho_l * pow( taw_w / rho_l , 0.5 ) * D / mu_l;
	T_iplus = rho_l * Cp_l * pow( taw_w / rho_l , 0.5 ) * ( T_i - T_wi ) / q;

	Nu_f = Re_lstar * Pr_l / T_iplus;
	Nu_f = 0.0152 * ( 1.0 + 0.6 * pow( Pr_l , 0.8 ) ) * ( fai_v / X_tt ) * pow( Re_l , 0.77 );

	m = 2.0;

	Nu = pow ( pow( Nu_f , m ) + pow( Nu_b , m ) , 1.0 / m );

	alpha = Nu * ( lambda_l / D );
	alpha /= 1000.0;

}

*/


/*
void fujii2::calc(){

	ref.state_ph2( P , h );
	ref.sat_p2( P );

	lambda_l = ref.Rl.thc;
	lambda_v = ref.Rv.thc;
	superheat = fabs( T_wi - ref.Rc.T );

	Cp = ref.Rl.cp * 1000.0;
	x = ref.Rc.x;

	//流路断面積
	S = D * D * atan ( 1.0 );

	//冷媒質量流速
	G2 = G / S;

	//比体積
	v_l = 1.0 / ref.Rl.rho;
	v_g = 1.0 / ref.Rv.rho;

	rho_l = ref.Rl.rho;
	rho_v = ref.Rv.rho;

	//粘度
	mu_l = ref.Rl.visc / 1000000;
	mu_v = ref.Rv.visc / 1000000;

	//乾き度
	x = ref.Rc.x;

	T_wi += 273.15;
	T_sat = ref.Rv.T + 273.15;
	Cp_l = ref.Rl.cp * 1000.0;
	Cp_v = ref.Rv.cp * 1000.0;


	//レイノルズ数
	Re_l = G2 * ( 1 - x ) * D / mu_l;
	Re = G2 * D / mu_l;

	Pr_l = Cp_l * mu_l / lambda_l;
	Pr_v = Cp_v * mu_v / lambda_v;

	g = 9.8;

	Ga = g * pow( rho_l , 2.0 ) * pow( D , 3.0 ) / pow( mu_l , 2 );


	hv = ( ref.Rv.h - ref.Rl.h ) * 1000.0;

	H_l = Cp_l * ( T_sat - T_wi ) / hv; 

	C5 = G / D / mu_l;

	C4 = 20.0 * exp( - C5 / 3000.0 );

	Ja = Cp_l / hv * superheat;
	Ja_v = Cp_v / hv * superheat;
	C3 = 0.47 * sqrt( v_g / v_l ) * pow( Ja / Pr_l , 1.0/12.0 ) * pow( Re_l * x/(1.0-x) , 0.9 )
		/pow( Ga * Pr_l / Ja , 1.1/4.0 );

	C2 = pow( 1.0 + 1.6 * 1e11 * pow( Ja / Pr_l , 5.0 ) , 0.25 ) / sqrt( v_g / v_l )
		* pow( pow( Ga * Pr_l / Ja , 0.25 ) / ( G*(1.0-x) / D / mu_l * x/(1.0-x) ) , 1.8 );

	C1 = 0.071 * pow( Re_l , 0.1 ) * pow( v_g / v_l , 0.55 ) * pow(x/(1.0-x) ,  0.2 - 0.1 * x)
		* pow( Pr_l + 8000.0 / pow( Re_l , 1.5 ) , 0.333333 );

	Nu_f = 0.018 * pow( Re_l * sqrt( v_g / v_l ) , 0.9 ) * pow( x/(1.0-x) , 0.1 * x + 0.8 )
		* pow( Pr_l + ( 8000.0 / pow( Re_l , 1.5 ) ) , 0.3333 ) * ( 1.0 + C1 * Ja / Pr_l  - 0.2 * Ja_v / Pr_v );

	Nu_b = 0.725 * pow( Ga * Pr_l / Ja , 0.25 ) * pow( 1.0 + 0.003 * sqrt( Pr_l ) * pow( C3 , 3.1 - 0.5 / Pr_l ) , 0.3 )
		/ pow( 1.0 + C2 * C4 , 0.25 );

	if( Nu_b > Nu_f ){
		Nu = Nu_b;
	}else{
		Nu = Nu_f;
	}

	alpha = Nu * ( lambda_l / D );
	alpha /= 1000.0;

	T_wi -= 273.15;

	if( 0.95 < x && x <= 1.0 ){
		alpha *= -sin( ( x - 1 ) / 0.1 * 4.0 * atan ( 1.0 ) );
	}
}
*/

/*
void dittus_boelter::calc(){

	if( kanetu == 0 ){
		n = 0.3;
	}else if( kanetu == 1 ){
		n = 0.4;
	}
	ref.state_ph2( P , h );
	lambda = ref.Rc.thc;
	mu = ref.Rc.visc / 1000000;
	S = D * D * atan(1.0);
	rho = ref.Rc.rho;
	u = G / rho / S;
	Cp = ref.Rc.cp * 1000;
	Re = u * D * rho / mu;
	Pr = mu * Cp / lambda;
	alpha = 0.023 * pow( Re , 0.8 ) * pow( Pr , n ) * ( lambda / D ); 
	alpha /= 1000.0;
}
*/

/*
void seshimo::calc(){

	///////////////////////////////////////////////////////////////////////////
	//		”低レイノルズ数域のプレートフィンチューブ　第1～3”
	//		瀬下　裕，藤井　雅雄
	//　2010.07 kikuchi
	//  大野改造
	//    変数名を瀬下らの論文と同じにした
	//    いったん霜なしにした
	//    一つのフィンのみにした
	//	
	//	****************************************************************
	//	<input>
	//	AO			:空気側伝熱面積[m^2]
	//	Di			:管内径[m]
	//	Do			:管外径[m]
	//	Fp			:フィンピッチ[m]
	//	Ft			:フィン厚さ[m]
	//	L			:熱交換器長さ[m]
	//	L_x			:列ピッチ[m]
	//	L_z			:段ピッチ[m]	
	//	e_fin		:フィン効率[-]
	//	pNum_y		:y方向配管列数[-]
	//
	//	L_f			:霜厚さ[m]
	//	ramda_f		:霜熱伝達率[kW/mK]
	//
	//	rho_air		:空気密度[kg/m^3]入口温度と出口温度の平均から計算
	//	rho_airi	:空気入口密度[kg/m^3]
	//	visc_air	:空気動粘性係数[uPa*s]
	//	v_air		:入口空気流速[m/s]
	//	thc_air		:入口空気熱伝導率[W/mK]
	//
	//			
	//		       
	//   　＼-------- -- 
	//	L	 /	   /| 　L_y
	//		/	  /||
	// ＼  /   z /||| --
	//	 -------/||/
	//	  |　　|||/	←Air
	//	  | ○ ||/ 
	//	  |　　|/x
	//	y--------
	//	  |    |
	//		L_x (L2)
	//
	//
	//		[x][y][z]
	//
	//
	//	<output>
	//	alpha	:冷媒側熱伝達率[W/m^2K]	
	//	
	///////////////////////////////////////////////////////////////////////////

//	double pi = 3.141592653589793;
//	double tau = eps / ( 1.0 - pow(1.0 - eps, 0.5) );////Mezedureの屈曲度

	//フィンカラー径
	Dc = Do + 2.0 * Ft;

	//フィン伝熱面積
	Ae = 2.0 * ( L / Fp ) * ( L_x * L_y - B.PI/4.0 * Dc * Dc );

	//管伝熱面積
	Ap = Dc * B.PI * L - Dc * B.PI * L / Fp * Ft;

	//全伝熱面積
	Ao = Ae + Ap;

	//全内伝熱面積
	Ai = L * Di * B.PI;

	//前面面積
	Af = L_y * L;

	//平均自由通過断面積
	Ac = Af * ( L_x * L_y - B.PI * Dc * Dc / 4.0 ) / ( L_x * L_y ) * ( Fp - Ft ) / Fp;

	//代表風速
	Vac = Af * rho_airi * v_air / ( Ac * rho_air );

	//代表寸法
	Dec = 4.0 * Ac * L_x / Ao;

	//レイノルズ数
	Re_air = Vac * Dec / ( visc_air / rho_air );


	//代表寸法
//	De_min = 4.0 * Ac * L_x / Ao;

	//レイノルズ数
//	Re_air_star = Va_max * De_min / ( visc_air / rho_air );

	//代表風速
//	Va_max = Af * rho_airi * v_air / ( Ac * rho_air );




	if( 0.5 <= pNum_y && pNum_y < 1.5 ){

		//ヌッセルト数（一列）
		Nu = 2.1 * pow( Re_air * Pr_air * Dec / L_x , 0.38 );
		//流動損失係数
		f = ( 0.43 + 35.1 * pow( Re_air * Dec / L_x , -1.07 )) * Dec / L_x;

	}else if( 1.5 <= pNum_y ){

		//ヌッセルト数（二列）
		Nu = 2.1 * pow( Re_air * Pr_air * Dec / L_x , 0.47 );
		//流動損失係数
		f = ( 0.83 + 24.7 * pow( Re_air * Dec / L_x , -0.89 )) * Dec / L_x;

	}


	//熱伝達率
	alpha = Nu * (thc_air) / Dec;

	//圧力損失
	dpdz = 2.0 * rho_air * pow(Vac,2.0) * f / Dec;


	//フィン相当直径
	Df = pow( 4.0 * L_x * L_y / B.PI , 0.5 );
	//fin効率
	e_fin = 1.0 / ( 1.0 + alpha * pow( Df - Dc, 2.0 ) * pow( Df/Dc , 0.5) / ( 6.0 * ramda_fin * Ft ) );





}


*/
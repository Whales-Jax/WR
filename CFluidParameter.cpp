#include "CFluidParameter.h"
#include <ctime>

CFluidParameter::CFluidParameter(){
	
	srand((unsigned) time(NULL));

	G = 0.050;
	P = 1000.0;
	h = 300.0;
	x = 0.0;
	X = 0.0011;
	s = 0.0;
	D = 0.0;
	u = 0.0;
	M = 0.0;
	Q = 0.0;
	q = 0.0;
	T = 30.0;
	v = 0.0;
	V = 0.0;
	velocity = 0.0;
	visc = 0.0;


	CycleNUM = 0;


/*

	ConnectModule = 0;
	ConnectPort = -1;
	IOCode = -1;*/

}
CFluidParameter::~CFluidParameter(){
}

bool CFluidParameter::valueCheck( CFluidParameter a ){

	if( G == a.G &&
		P == a.P &&
		h == a.h &&
		x == a.x &&
		X == a.X &&
		s == a.s &&
		D == a.D &&
		u == a.u &&
		M == a.M &&
		Q == a.Q &&
		q == a.q &&
		T == a.T &&
		v == a.v &&
		V == a.V &&
		velocity == a.velocity &&
		visc == a.visc&&
		Go == a.Go &&
		Gi == a.Gi  ){
		return true;
	}else{
		return false;
	}
	return false;
}


void CFluidParameter::setOutletConnection( int module , int port ){

	ConnectModule = module;
	ConnectPort = port;

	return;
}

void CFluidParameter::copy( CFluidParameter a ){
	G = a.G;
	P = a.P;
	h = a.h;
	x = a.x;
	X = a.X;
	s = a.s;
	D = a.D;
	u = a.u;
	M = a.M;
	Q = a.Q;
	q = a.q;
	T = a.T;
	v = a.v;
	V = a.V;
	velocity = a.velocity;
	visc = a.visc;

	alpha = a.alpha;
	pd = a.pd;
	dp = a.dp;

	Gi = a.Gi;
	Pi = a.Pi;
	hi = a.hi;
	Xi = a.Xi;
	Gj = a.Gj;
	Pj = a.Pj;
	hj = a.hj;
	Xj = a.Xj;
	Go = a.Go;
	Po = a.Po;
	ho = a.ho;
	Xo = a.Xo;
	Gp = a.Gp;
	Pp = a.Pp;
	hp = a.hp;
	Xp = a.Xp;

}

void CFluidParameter::subset( DXFluid a ){

	P    = a.P;
	h    = a.h;
	D    = a.rho;
	T    = a.T;
	s    = a.s;
	x    = a.x;
	u    = a.h - a.P / a.rho;
	cp   = a.cp;
	thc  = a.thc;
	visc = a.visc;
	st   = a.st;
	cv   = a.cv;
	
}

void CFluidParameter::subset_i( DXFluid a ){

	Pi    = a.P;
	hi    = a.h;
	Di    = a.rho;
	Ti    = a.T;
	si    = a.s;
	
}

void CFluidParameter::subset_j( DXFluid a ){

	Pj    = a.P;
	hj    = a.h;
	Dj    = a.rho;
	Tj    = a.T;
	sj    = a.s;
	
}

void CFluidParameter::subset_o( DXFluid a ){

	Po    = a.P;
	ho    = a.h;
	Do    = a.rho;
	To    = a.T;
	so    = a.s;
	
}

void CFluidParameter::subset_p( DXFluid a ){

	Pp    = a.P;
	hp    = a.h;
	Dp    = a.rho;
	Tp    = a.T;
	sp    = a.s;
	
}

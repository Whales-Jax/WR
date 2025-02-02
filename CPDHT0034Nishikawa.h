#ifndef __CPDHT0034Nishikawa_h
#define __CPDHT0034Nishikawa_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
using namespace std;

#include "DXProperty_ver06.h"
#include "CFluidParameter.h"
#include "CStructureParameter.h"
#include "CNewtonRaphsonMethod.h"

class CPDHT0034Nishikawa{
public :
	CBasicParameter		B;
	CFluidParameter		*fp;
	CStructureParameter *sp;

	double X;
	double f;
	double M;       // 900 [m-1]
	double p;		// ���x�`����[m2/s]
	double Pa;		// ��C��[kPa]
	double P;		// �C��[kPa]
	double h_fg;	// �������M[kJ/kg]
	double temp_alpha;	// �M�`�B��(����l)W/m2////////


	CPDHT0034Nishikawa();
	~CPDHT0034Nishikawa();

	void calc();
	void calc2();
};


#endif
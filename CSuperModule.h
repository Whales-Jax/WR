#ifndef __CSuperModule_h_
#define __CSuperModule_h_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>

#include "CStructureParameter.h"
#include "CFileoutControlModule.h"
#include "DXProperty_ver06.h"


using namespace std;

class CValuePointerModule{
public:
	string name;
	int ID;//0:int 1:double
	int *IValue;
	double *DValue;
	char* SValue;
};
class CValuePointer{
public:
	CValuePointer(){
		Pnum = 0;
	}
	~CValuePointer(){}

	void setPointer( int _ID , string _name , int& _Value ){
		pushback();
		VPM[ Pnum-1 ].ID = _ID;
		VPM[ Pnum-1 ].name = _name;
		VPM[ Pnum-1 ].IValue = &_Value;
	}
	void setPointer( int _ID , string _name , double& _Value ){
		pushback();
		VPM[ Pnum-1 ].ID = _ID;
		VPM[ Pnum-1 ].name = _name;
		VPM[ Pnum-1 ].DValue = &_Value;
	}
	void setPointer( int _ID , string _name , char* _Value ){
		pushback();
		VPM[ Pnum-1 ].ID = _ID;
		VPM[ Pnum-1 ].name = _name;
		VPM[ Pnum-1 ].SValue = _Value;
	}
	void pushback(){
		Pnum++;
		CValuePointerModule temp;
		VPM.push_back( temp );
		return;
	}

	int Pnum;
	vector<CValuePointerModule> VPM;
};



class CSuperModule{
public:
	CSuperModule(){
		ErrorCode=1;

		er.resize(30);
		for( int i = 0 ; i < (int)er.size() ; i++ ){
			er[i].resize(2);
		}

	}
	~CSuperModule(){}

	CStructureParameter _P;
	CStructureParameter P;
	CStructureParameter P1;
	CStructureParameter P2;
	CStructureParameter P3;
	CBasicParameter B;

	string name;
	char name2[128];
	char code[128];
	int ErrorCode;


	vector <vector<double>> er;
	int en;
	void set_error( int i , double a, double b){//ÉÅÉ\ÉbÉhÇÃêÈåæ
		er[i][0] = a;
		er[i][1] = b;
	}

	//iniÇ©ÇÁéÛÇØìnÇ∑É|ÉCÉìÉ^Å[
	CValuePointer CVP;

	double t;
	double dt;
	double M;//è[ìUó 
	double MX;
	double *Mass;//è[ìUó 

	int DriveMode;

	int nDXP;
	DXProperty *DXP;
	
	int R134a;
	int R410A;

	CFileoutControlModule CFCM;


	int unsteady;
	int Htransfer;
	int Pdrop;

	double ht_keisu;
	double pd_keisu;



};

#endif
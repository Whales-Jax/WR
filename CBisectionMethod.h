#ifndef _CBisectionMethod_h_
#define _CBisectionMethod_h_

#include <cmath>
#include <ctime>


class CBisectionMethod {
public:
	int LoopCount;	//ループ回数
	int SubLoopCount;	//ループ回数
	double Error;
	double Value[3];
	double Answer[3];
	double ReturnValue;
	double StartValue[3];

	CBisectionMethod();
	~CBisectionMethod();
	void setupError( double _Error = 1e-6 );
	void init();
	bool check();
	void setValue( double LValue , double RValue );
	double getValue( void );
	void setError( double Error );
	void prt();
};


class CBisectionMethod02 {
public:
	int LoopCount;	//ループ回数
	int SubLoopCount;	//ループ回数
	int area;
	double Error;
	double Value[2][3];
	double Answer[2][3][3];
	double StartValue[2][3];
	double ReturnValue[2];

	CBisectionMethod02();
	~CBisectionMethod02();
	void setupError( double _Error = 1e-6 );
	void init();
	bool check();
	bool check2();
	bool check2_sub(double,double,double,double);
	void setValue( int i , double LValue , double RValue );
	double getValue( int i );
	void setError( int i , double Error );
	void prt();
};



#endif
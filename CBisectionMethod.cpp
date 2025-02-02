#include <iostream>
#include "CBisectionMethod.h"
/**************************************************
二分法　2009/08/07　大野
使い方の例
y = x^3 - 11　の根を求める場合
#include <iostream>
#include <cmath>
#include "CBisectionMethod.h"

void main(){
	double x;
	double y;
	CBisectionMethod cbm;				//宣言します
	cbm.setupError( 1e-8 );				//許容誤差をセット（なくても可，その場合はError=1e-6）
	cbm.setValue( 0 , 100 );			//二分法ではさむ上限下限を入力
	for(cbm.init();cbm.check();){		//おまじない
		x = cbm.getValue();				//二分法の修正値を受け取る関数
		y = pow( x , 3 ) - 11.0;		//このあたりで計算
		cbm.setError( y );				//エラーセット　このyの値が0になるように二分法します
		cbm.prt();						//エラー表示，コメント推奨
	}
	return;
}


二変数の二分法も作ろうと思ったけど，やめました．
そのためこのcppファイルには結構いらない部分もあります．
そのうち消します　大野

****************************************************/
CBisectionMethod::CBisectionMethod(){

	srand((unsigned) time(NULL));

	LoopCount = 0;
	SubLoopCount = 0;
	Error = 1e-6;

	Value[0] = 0.0;
	Value[1] = 0.0;
	Value[2] = 0.0;

	Answer[0] = 1000.0;
	Answer[1] = 1000.0;
	Answer[2] = 1000.0;

	ReturnValue = 0.0;

}

CBisectionMethod::~CBisectionMethod(){
}

void CBisectionMethod::setupError( double _Error ){
	Error = _Error;
}

void CBisectionMethod::setValue( double LValue, double RValue ){
	if( LValue > RValue ){
		double temp;
		temp = LValue;
		LValue = RValue;
		RValue = temp;
	}


	Value[0] = LValue;
	Value[2] = RValue;
	Value[1] = ( Value[0] + Value[2] ) / 2.0;

	StartValue[0] = Value[0];
	StartValue[1] = Value[1];
	StartValue[2] = Value[2];
}


double CBisectionMethod::getValue( void ){
	return ReturnValue;
}

void CBisectionMethod::setError( double Error ){

	if( LoopCount == 1 ){
		Answer[0] = Error; 
	}else if( LoopCount == 2 ){
		Answer[2] = Error;
	}else{
		Answer[1] = Error;
	}

}

void CBisectionMethod::init(){
	LoopCount = 0;
	SubLoopCount = 0;
	Answer[0] = -1000.0;
	Answer[1] = 1000.0;
	Answer[2] = 1000.0;
}

bool CBisectionMethod::check(){

	LoopCount++;

	if( LoopCount == 4 ){

		if( fabs( Answer[1] ) < Error ){
			return false;
		}

		if( Answer[0]*Answer[2] > 0.0 ){
			std::cout << "二分法初期値Error" << std::endl;
			exit(1);
		}
	}

	if( LoopCount > 3 ){
		if( Answer[0]*Answer[1] < 0.0 ){
			Answer[2] = Answer[1];
			Value[2] = Value[1];
			Value[1] = ( Value[0] + Value[2] ) / 2.0 ;
		}else{
			Answer[0] = Answer[1];
			Value[0] = Value[1];
			Value[1] = ( Value[0] + Value[2] ) / 2.0 ;
		}
	}

	if( LoopCount == 1 ){
		ReturnValue = Value[0];
	}else if( LoopCount == 2 ){
		ReturnValue = Value[2];
	}else{
		ReturnValue = Value[1];
	}


	if( fabs( Answer[1] ) < Error ){
		return false;
	}
	return true;
}

void CBisectionMethod::prt(){

	std::cout << " Loop = " << LoopCount;
	std::cout << " Error = " << Answer[1];
	std::cout << std::endl;

}



CBisectionMethod02::CBisectionMethod02(){

	srand((unsigned) time(NULL));

	LoopCount = 0;
	SubLoopCount = 0;
	Error = 1e-6;

	Value[0][0] = 0.0;
	Value[0][1] = 0.0;
	Value[0][2] = 0.0;
	Value[1][0] = 0.0;
	Value[1][1] = 0.0;
	Value[1][2] = 0.0;

	Answer[0][0][0] = 1000.0;
	Answer[0][0][1] = 1000.0;
	Answer[0][0][2] = 1000.0;
	Answer[0][1][0] = 1000.0;
	Answer[0][1][1] = 1000.0;
	Answer[0][1][2] = 1000.0;
	Answer[0][2][0] = 1000.0;
	Answer[0][2][1] = 1000.0;
	Answer[0][2][2] = 1000.0;
	Answer[1][0][0] = 1000.0;
	Answer[1][0][1] = 1000.0;
	Answer[1][0][2] = 1000.0;
	Answer[1][1][0] = 1000.0;
	Answer[1][1][1] = 1000.0;
	Answer[1][1][2] = 1000.0;
	Answer[1][2][0] = 1000.0;
	Answer[1][2][1] = 1000.0;
	Answer[1][2][2] = 1000.0;

	ReturnValue[0] = 0.0;
	ReturnValue[1] = 0.0;

}

CBisectionMethod02::~CBisectionMethod02(){
}

void CBisectionMethod02::setupError( double _Error ){
	Error = _Error;
}

void CBisectionMethod02::setValue( int i , double LValue, double RValue ){

	Value[i][0] = LValue;
	Value[i][2] = RValue;
	Value[i][1] = ( Value[i][0] + Value[i][2] ) / 2.0;

}


double CBisectionMethod02::getValue( int i ){
	return ReturnValue[i];
}

void CBisectionMethod02::setError( int i , double Error ){


	if( LoopCount == 1 ){
		if( SubLoopCount == 1 ){
			ReturnValue[0] = 
			Answer[i][0][0] = Error;
		}else if( SubLoopCount == 2 ){
			Answer[i][0][1] = Error;
		}else if( SubLoopCount == 3 ){
			Answer[i][0][2] = Error;
		}else if( SubLoopCount == 4 ){
			Answer[i][1][0] = Error;
		}else if( SubLoopCount == 5 ){
			Answer[i][1][1] = Error;
		}else if( SubLoopCount == 6 ){
			Answer[i][1][2] = Error;
		}else if( SubLoopCount == 7 ){
			Answer[i][2][0] = Error;
		}else if( SubLoopCount == 8 ){
			Answer[i][2][1] = Error;
		}else if( SubLoopCount == 9 ){
			Answer[i][2][2] = Error;
		}
	}else{
		if( SubLoopCount == 1 ){
			Answer[i][0][1] = Error;
		}else if( SubLoopCount == 2 ){
			Answer[i][1][0] = Error;
		}else if( SubLoopCount == 3 ){
			Answer[i][1][1] = Error;
		}else if( SubLoopCount == 4 ){
			Answer[i][1][2] = Error;
		}else if( SubLoopCount == 5 ){
			Answer[i][2][1] = Error;
		}
	}


}

void CBisectionMethod02::init(){
	LoopCount = 1;
	SubLoopCount = 0;
	area = 0;
	Answer[0][0][0] = -1000.0;
	Answer[0][0][1] = -1000.0;
	Answer[0][0][2] = -1000.0;
	Answer[0][1][0] = 1000.0;
	Answer[0][1][1] = 1000.0;
	Answer[0][1][2] = 1000.0;
	Answer[0][2][0] = 1000.0;
	Answer[0][2][1] = 1000.0;
	Answer[0][2][2] = 1000.0;
	Answer[1][0][0] = -1000.0;
	Answer[1][0][1] = -1000.0;
	Answer[1][0][2] = -1000.0;
	Answer[1][1][0] = 1000.0;
	Answer[1][1][1] = 1000.0;
	Answer[1][1][2] = 1000.0;
	Answer[1][2][0] = 1000.0;
	Answer[1][2][1] = 1000.0;
	Answer[1][2][2] = 1000.0;
}

bool CBisectionMethod02::check2_sub( double a1 , double a2 , double b1 , double b2 ){

	if( a1 < 0.0 && a2 < 0.0 && b1 >= 0.0 && b2 >= 0.0 ){
		return true;
	}else if( a1 >= 0.0 && a2 >= 0.0 && b1 < 0.0 && b2 < 0.0 ){
		return true;
	}else{
		return false;
	}

}


bool CBisectionMethod02::check2(){

	bool flag = false;

	flag = check2_sub( Answer[0][0][0] , Answer[1][0][0] , Answer[0][1][1] , Answer[1][1][1] );
	if( flag == true ){
		area = 1;
	}
	flag = check2_sub( Answer[0][0][1] , Answer[1][0][1] , Answer[0][1][0] , Answer[1][1][0] );
	if( flag == true ){
		area = 1;
	}

	flag = check2_sub( Answer[0][0][1] , Answer[1][0][1] , Answer[0][1][2] , Answer[1][1][2] );
	if( flag == true ){
		area = 2;
	}
	flag = check2_sub( Answer[0][0][2] , Answer[1][0][2] , Answer[0][1][1] , Answer[1][1][1] );
	if( flag == true ){
		area = 2;
	}

	flag = check2_sub( Answer[0][1][0] , Answer[1][1][0] , Answer[0][2][1] , Answer[1][2][1] );
	if( flag == true ){
		area = 3;
	}
	flag = check2_sub( Answer[0][1][1] , Answer[1][1][1] , Answer[0][2][0] , Answer[1][2][0] );
	if( flag == true ){
		area = 3;
	}

	flag = check2_sub( Answer[0][1][1] , Answer[1][1][1] , Answer[0][2][2] , Answer[1][2][2] );
	if( flag == true ){
		area = 4;
	}
	flag = check2_sub( Answer[0][1][2] , Answer[1][1][2] , Answer[0][2][1] , Answer[1][2][1] );
	if( flag == true ){
		area = 4;
	}

	if( area == 1 ){
		Answer[0][2][2] = Answer[0][1][1]; 
		Answer[1][2][2] = Answer[1][1][1]; 

		Answer[0][0][2] = Answer[0][0][1]; 
		Answer[1][0][2] = Answer[1][0][1]; 

		Answer[0][2][0] = Answer[0][1][0]; 
		Answer[1][2][0] = Answer[1][1][0];

		Value[0][2] = Value[0][1];
		Value[1][2] = Value[1][1];
	}
	if( area == 2 ){
		Answer[0][2][2] = Answer[0][1][2]; 
		Answer[1][2][2] = Answer[1][1][2]; 

		Answer[0][0][0] = Answer[0][0][1]; 
		Answer[1][0][0] = Answer[1][0][1]; 

		Answer[0][2][0] = Answer[0][1][1]; 
		Answer[1][2][0] = Answer[1][1][1]; 

		Value[0][0] = Value[0][1];
		Value[1][2] = Value[1][1];
	}
	if( area == 3 ){
		Answer[0][0][2] = Answer[0][1][1]; 
		Answer[1][0][2] = Answer[1][1][1]; 

		Answer[0][0][0] = Answer[0][1][0]; 
		Answer[1][0][0] = Answer[1][1][0]; 

		Answer[0][2][2] = Answer[0][2][1]; 
		Answer[1][2][2] = Answer[1][2][1]; 

		Value[0][2] = Value[0][1];
		Value[1][0] = Value[1][1];
	}
	if( area == 4 ){
		Answer[0][0][2] = Answer[0][1][2]; 
		Answer[1][0][2] = Answer[1][1][2]; 

		Answer[0][0][0] = Answer[0][1][1]; 
		Answer[1][0][0] = Answer[1][1][1]; 

		Answer[0][2][0] = Answer[0][2][1]; 
		Answer[1][2][0] = Answer[1][2][1]; 

		Value[0][0] = Value[0][1];
		Value[1][0] = Value[1][1];
	}

	Value[0][1] = ( Value[0][0] + Value[0][2] ) / 2.0;
	Value[1][1] = ( Value[1][0] + Value[1][2] ) / 2.0;


	return true;

}

bool CBisectionMethod02::check(){

	if( LoopCount == 1 && SubLoopCount == 9 ){
		LoopCount++;
		SubLoopCount = 0;
		check2();
	}
	if( LoopCount > 1 && SubLoopCount == 5 ){
		LoopCount++;
		SubLoopCount = 0;
		check2();
	}



	SubLoopCount++;


	if( LoopCount == 1 ){
		if( SubLoopCount == 1 ){
			ReturnValue[0] = Value[0][0];
			ReturnValue[1] = Value[1][0];
		}else if( SubLoopCount == 2 ){
			ReturnValue[0] = Value[0][0];
			ReturnValue[1] = Value[1][1];
		}else if( SubLoopCount == 3 ){
			ReturnValue[0] = Value[0][0];
			ReturnValue[1] = Value[1][2];
		}else if( SubLoopCount == 4 ){
			ReturnValue[0] = Value[0][1];
			ReturnValue[1] = Value[1][0];
		}else if( SubLoopCount == 5 ){
			ReturnValue[0] = Value[0][1];
			ReturnValue[1] = Value[1][1];
		}else if( SubLoopCount == 6 ){
			ReturnValue[0] = Value[0][1];
			ReturnValue[1] = Value[1][2];
		}else if( SubLoopCount == 7 ){
			ReturnValue[0] = Value[0][2];
			ReturnValue[1] = Value[1][0];
		}else if( SubLoopCount == 8 ){
			ReturnValue[0] = Value[0][2];
			ReturnValue[1] = Value[1][1];
		}else if( SubLoopCount == 9 ){
			ReturnValue[0] = Value[0][2];
			ReturnValue[1] = Value[1][2];
		}


	}
	if( LoopCount > 1 ){
		if( SubLoopCount == 1 ){
			ReturnValue[0] = Value[0][0];
			ReturnValue[1] = Value[1][1];
		}else if( SubLoopCount == 2 ){
			ReturnValue[0] = Value[0][1];
			ReturnValue[1] = Value[1][0];
		}else if( SubLoopCount == 3 ){
			ReturnValue[0] = Value[0][1];
			ReturnValue[1] = Value[1][1];
		}else if( SubLoopCount == 4 ){
			ReturnValue[0] = Value[0][1];
			ReturnValue[1] = Value[1][2];
		}else if( SubLoopCount == 5 ){
			ReturnValue[0] = Value[0][2];
			ReturnValue[1] = Value[1][1];
		}

	}


	if( fabs( Answer[0][1][1] ) < Error && fabs( Answer[1][1][1] ) < Error ){
		return false;
	}
	return true;
}

void CBisectionMethod02::prt(){

	if( (LoopCount == 1 && SubLoopCount == 9) || (LoopCount > 1 && SubLoopCount == 5) ){

		std::cout << " Loop = " << LoopCount;
		std::cout << " Error01 = " << Answer[0][1][1];
		std::cout << " Error02 = " << Answer[1][1][1];
		std::cout << std::endl;
	}

}


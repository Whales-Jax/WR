

#ifndef __CMASTER222_H__
#define __CMASTER222_H__

#include <iostream>
#include <sstream>
#include <fstream>

/*!
2008/06/08　5108C015大野
xに対応するyの値，yに対応するxの値を返します．

使い方
以下のようなcsvファイルを用意してください．
//ここから
0,5
1,2
2,4
3,7
4,2
5,8
6,1
end
//ここまで
一列目がxの値，二列目がyの値として認識されます．

以下使い方の例
#include <iostream>
#include "CMaster222.h"
int main(void){
	CMaster222 cma;
	cma.setup("sample.csv");//セットアップ
	double y;
	y = cma.xtoy( 3.5 );//xが3.5のときのyの値が返される．yの値は比例補間で計算される．
	return 1;
}
*/
class CMaster222
{
public:

	CMaster222(){}
	~CMaster222(){}

	void setup(std::string FileName);
	double xtoy(double a);
	double **mt;
	int num;

};

#endif
#include <iostream>
#include <sstream>
#include <fstream>
#include "CMaster222.h"


void CMaster222::setup(std::string FileName){
	std::string buf;
	int i;
	int j;

	//テキストファイル読み込み
	std::fstream desp( FileName.c_str() );
	if( desp.fail() ){
		std::cout << "file open error" << std::endl; exit(1);
	}

	//行数取得
	i = 0;
	while( desp && getline( desp , buf ) ){
		i++;
		if( buf == "end," || buf == "end" ){
			break;
		}
	}
	num = i;

	//配列宣言
	mt = new double *[2];
	for( j = 0 ; j < 2 ; j++ ){
		mt[j] = new double [i];
	}

	desp.close();


	std::fstream desp2( FileName.c_str() );
	if( desp2.fail() ){
		std::cout << "file open error" << std::endl; exit(1);
	}



	//値読み込み
	desp2.seekg( 0 , std::ios::beg );
	for( j = 0 ; j < i ; j++ ){
		getline( desp2 , buf );
		mt[0][j] = atof( buf.substr( 0 , buf.find(",") ).c_str() );
		mt[1][j] = atof( buf.substr( buf.find(",")+1 , buf.length() ).c_str() );

	}
}


double CMaster222::xtoy( double a ){

	int i;
	int j;
	int k;
	double ans;

	i = 0;
	j = 1;


	for( k = 0 ; k < num-1 ; k++ ){
		if( mt[i][k] <= a && a < mt[i][k+1] ){
			break;
		}
	}

	ans = ( mt[j][k+1] - mt[j][k] ) / ( mt[i][k+1] - mt[i][k] ) * ( a - mt[i][k] ) + mt[j][k] ;
	return ans;
}


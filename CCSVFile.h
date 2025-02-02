#ifndef __CCSVFILE_H__
#define __CCSVFILE_H__

#include <windows.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>


#define CCSVFILE_is_sjis1(c) (((((unsigned char)(c))>=0x81)&&(((unsigned char)(c))<=0x9F))||((((unsigned char)(c))>=0xE0)&&(((unsigned char)(c))<=0xFC)))
#define CCSVFILE_is_sjis2(c) ((((unsigned char)(c))!=0x7F)&&(((unsigned char)(c))>=0x40)&&(((unsigned char)(c))<=0xFC))

#define CCSVFILE_NOBUFFER() (rowsize == 0L && colsize == 0L && datumsize == 0L && buffer == (char***)NULL)

class CCSVFile
{
	FILE			*fp;

	unsigned long	rowsize;	// 行数
	unsigned long	colsize;	// 列数
	unsigned long	datumsize;	// 一データの最大長

	char			***buffer;

	bool			CheckCSV();
	bool			ReadData();

public:
	CCSVFile();
	~CCSVFile();

	// ファイルからCSVデータ読み出し
//	bool			LoadBuffer(char *filename);
	bool			LoadBuffer(std::string filename);
	// 読み出したデータ解放
	void			ReleaseBuffer();

/*	bool			WriteBuffer(char *filename, double **data, unsigned long rows, unsigned long cols);
	bool			WriteBuffer(char *filename, long **data, unsigned long rows, unsigned long cols);
	bool			WriteBuffer(char *filename, int **data, unsigned long rows, unsigned long cols);*/
	// バッファの内容をファイルへCSVデータとして書き込み
	bool			WriteBuffer(char *filename);
	// 指定データをファイルへCSVデータとして書き込み
	bool			WriteBuffer(char *filename, char ***data, unsigned long rows, unsigned long cols, unsigned long datumlen);

	// 各プライベート変数読み出し
	char			***GetCSVBuffer(){ return buffer; };
	char			**GetCSVBuffer(unsigned long row){ return(row < rowsize)?buffer[row]:(char**)NULL; };
	char			*GetCSVBuffer(unsigned long row, unsigned long col){ return(row < rowsize && col < colsize)?buffer[row][col]:(char*)NULL; };
	unsigned long	GetRowSize(){ return rowsize; };
	unsigned long	GetColSize(){ return colsize; };
	unsigned long	GetDatumSize(){ return datumsize; };
};

#endif
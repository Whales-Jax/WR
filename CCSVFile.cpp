#include "CCSVFile.h"
using namespace std;


//csv解析用
CCSVFile::CCSVFile()
{
	rowsize = 0L;
	colsize = 0L;
	datumsize = 0L;
	buffer = (char***)NULL;
};

CCSVFile::~CCSVFile()
{
	ReleaseBuffer();
};

// CSVファイル読み込み
bool CCSVFile::LoadBuffer(std::string filename)
//bool CCSVFile::LoadBuffer(char *filename)
{
	unsigned long i, j;

	// 既に何らかのデータを読み込んであり、解放されていない
	if(!CCSVFILE_NOBUFFER())
		return false;

	// ファイルオープン
	if( fopen_s( &fp, filename.c_str(), "r") != 0)
		return false;

	// プロパティ値取得
	if(!CheckCSV())
		return false;

	// バッファ確保
	buffer = new char**[rowsize];
	for(i=0;i<rowsize;i++)
	{
		buffer[i] = new char*[colsize];
		for(j=0;j<colsize;j++)
			buffer[i][j] = new char[datumsize];
	}

	// データ読み込み
	if(!ReadData())
		return false;

	//ファイルクローズ
	fclose(fp);

	return true;
}

// 全メモリ解放
void CCSVFile::ReleaseBuffer()
{
	if(CCSVFILE_NOBUFFER())
		return;

	unsigned long i, j;

	// バッファの解放
	for(i=0;i<rowsize;i++)
	{
		for(j=0;j<colsize;j++)
			delete buffer[i][j];
		delete buffer[i];
	}
	delete buffer;

	// 値のリセット
	rowsize = colsize = datumsize = 0L;
	buffer = (char***)NULL;
	return;
}

bool CCSVFile::CheckCSV()
{
	bool quated = false;
	bool first = true;
	unsigned long rn = 0, cn = 0, dn = 0;
	int state = 0;	 // { 0, 1, 2 } = { 1byte文字, 2byte文字の第1byte, 2byte文字の第2byte }
	int c;	// 文字

	// データが既に読み込んである場合、解放する
	ReleaseBuffer();

	do
	{
		c = getc(fp);

		// 文字ステートメントの変更
		if      ( ( state == 0 ) && ( CCSVFILE_is_sjis1( c ) ) )	state = 1; // 0 -> 1
		else if ( ( state == 1 ) && ( CCSVFILE_is_sjis2( c ) ) )	state = 2; // 1 -> 2
		else if ( ( state == 2 ) && ( CCSVFILE_is_sjis1( c ) ) )	state = 1; // 2 -> 1
		else										state = 0; // 2 -> 0, その他

		// CSVでのフォーマットは1byte文字に依る
		if(state == 0)
		{
			if(c == '"')
				quated = !quated;
			if(quated)
			{
				if(c != '"' || (quated && !first))
					if(datumsize < ++dn)
						datumsize = dn;
				first = false;
			}
			else
			{
				if(c == ',')
				{
					if(colsize < ++cn)
						colsize = cn;
					dn = 0;
					first = true;
				}
				else if(c == '\n')
				{
					if(colsize < ++cn)
						colsize = cn;
					if(rowsize < ++rn)
						rowsize = rn;
					cn = 0;
					dn = 0;
					first = true;
				}
				else
				{
					if(c != '"' || quated)
						if(datumsize < ++dn)
							datumsize = dn;
					first = false;
				}
			}
		}
	}while(!feof(fp));

	rewind(fp);

	// 終末ヌル文字'\0'冗長
	datumsize++;

	return true;
}

bool CCSVFile::ReadData()
{
	bool quated = false;
	bool first = true;
	unsigned long rn = 0, cn = 0, dn = 0;
	int state = 0;	 // { 0, 1, 2 } = { 1byte文字, 2byte文字の第1byte, 2byte文字の第2byte }
	int c;	// 文字

	do
	{
		c = getc(fp);
		if(c == EOF)
			break;

		// 文字ステートメントの変更
		if      ( ( state == 0 ) && ( CCSVFILE_is_sjis1( c ) ) )	state = 1; // 0 -> 1
		else if ( ( state == 1 ) && ( CCSVFILE_is_sjis2( c ) ) )	state = 2; // 1 -> 2
		else if ( ( state == 2 ) && ( CCSVFILE_is_sjis1( c ) ) )	state = 1; // 2 -> 1
		else										state = 0;	// 2 -> 0, その他

		// CSVでのフォーマットは1byte文字に依る
		if(state == 0)
		{
			if(c == '"')
				quated = !quated;
			if(quated)
			{
				if(c != '"' || quated)
					if(c != '"' || !first)
						buffer[rn][cn][dn++] = (char)c;
				first = false;
			}
			else
			{
				if(c == ',')
				{
					buffer[rn][cn++][dn] = '\0';
					dn = 0;
					first = true;
				}
				else if(c == '\n')
				{
					buffer[rn++][cn][dn] = '\0';
					cn = 0;
					dn = 0;
					first = true;
				}
				else
				{
					if(c != '"' || quated)
						buffer[rn][cn][dn++] = (char)c;

					first = false;
				}
			}
		}
	}while(!feof(fp));

	return true;
}

bool CCSVFile::WriteBuffer(char *filename)
{
	unsigned long rn, cn, dn;
	char tmp;
	bool errChar;
	int state;
 
	// バッファがない
	if(CCSVFILE_NOBUFFER())
		return false;

	if( fopen_s( &fp , filename, "a") != 0 )
		return false;

	for(rn=0;rn<rowsize;rn++)
	{
		for(cn=0;cn<colsize;cn++)
		{
			// 違反文字が含まれるチェック
			errChar = false;
			state = 0;
			for(dn=0;dn<datumsize && buffer[rn][cn][dn] != '\0';dn++)
			{
				// 文字ステートメントの変更
				if      ( ( state == 0 ) && ( CCSVFILE_is_sjis1( buffer[rn][cn][dn] ) ) )	state = 1; // 0 -> 1
				else if ( ( state == 1 ) && ( CCSVFILE_is_sjis2( buffer[rn][cn][dn] ) ) )	state = 2; // 1 -> 2
				else if ( ( state == 2 ) && ( CCSVFILE_is_sjis1( buffer[rn][cn][dn] ) ) )	state = 1; // 2 -> 1
				else															state = 0; // 2 -> 0, その他

				if(state == 0 && ( (tmp = buffer[rn][cn][dn]) == '\n' || tmp == '\r' || tmp == ',' || tmp == '"' ) )
				{
					errChar = true;
					break;
				}
			}

			// 違反文字が含まれる⇒クォートする
			errChar && putc('"', fp);

			state = 0;
			// ダブルクォーテーションは二つ重ねる処理を加えながら文字列書き出し
			for(dn=0;dn<datumsize && buffer[rn][cn][dn] != '\0';dn++)
			{
				// 文字ステートメントの変更
				if      ( ( state == 0 ) && ( CCSVFILE_is_sjis1( buffer[rn][cn][dn] ) ) )	state = 1; // 0 -> 1
				else if ( ( state == 1 ) && ( CCSVFILE_is_sjis2( buffer[rn][cn][dn] ) ) )	state = 2; // 1 -> 2
				else if ( ( state == 2 ) && ( CCSVFILE_is_sjis1( buffer[rn][cn][dn] ) ) )	state = 1; // 2 -> 1
				else															state = 0; // 2 -> 0, その他

				if(state == 0 && buffer[rn][cn][dn] == '"')
					putc('"', fp);
				putc(buffer[rn][cn][dn], fp);
			}

			// 違反文字が含まれる⇒クォートする
			errChar && putc('"', fp);
		}
	}
	fclose(fp);

	return true;
}

bool CCSVFile::WriteBuffer(char *filename, char ***data, unsigned long rows, unsigned long cols, unsigned long datumlen)
{
	unsigned long rn, cn, dn;
	if(rows == 0L || cols == 0L || datumlen == 0L || data == (char***)NULL)
		return false;

	// 既にデータがあった場合、解放する
	ReleaseBuffer();

	// 各設定値設定
	rowsize = rows;
	colsize = cols;
	datumsize = datumlen;
	buffer = new char**[rows];
	for(rn=0;rn<rowsize;rn++)
	{
		buffer[rn] = new char*[cols];
		for(cn=0;cn<colsize;cn++)
		{
			buffer[rn][cn] = new char[datumlen];
			for(dn=0;dn<datumsize;dn++)
				buffer[rn][cn][dn] = data[rn][cn][dn];
		}
	}

	if(filename != NULL)
		if(!WriteBuffer(filename))
			return false;
	return true;
}

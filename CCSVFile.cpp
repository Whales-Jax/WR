#include "CCSVFile.h"
using namespace std;


//csv��͗p
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

// CSV�t�@�C���ǂݍ���
bool CCSVFile::LoadBuffer(std::string filename)
//bool CCSVFile::LoadBuffer(char *filename)
{
	unsigned long i, j;

	// ���ɉ��炩�̃f�[�^��ǂݍ���ł���A�������Ă��Ȃ�
	if(!CCSVFILE_NOBUFFER())
		return false;

	// �t�@�C���I�[�v��
	if( fopen_s( &fp, filename.c_str(), "r") != 0)
		return false;

	// �v���p�e�B�l�擾
	if(!CheckCSV())
		return false;

	// �o�b�t�@�m��
	buffer = new char**[rowsize];
	for(i=0;i<rowsize;i++)
	{
		buffer[i] = new char*[colsize];
		for(j=0;j<colsize;j++)
			buffer[i][j] = new char[datumsize];
	}

	// �f�[�^�ǂݍ���
	if(!ReadData())
		return false;

	//�t�@�C���N���[�Y
	fclose(fp);

	return true;
}

// �S���������
void CCSVFile::ReleaseBuffer()
{
	if(CCSVFILE_NOBUFFER())
		return;

	unsigned long i, j;

	// �o�b�t�@�̉��
	for(i=0;i<rowsize;i++)
	{
		for(j=0;j<colsize;j++)
			delete buffer[i][j];
		delete buffer[i];
	}
	delete buffer;

	// �l�̃��Z�b�g
	rowsize = colsize = datumsize = 0L;
	buffer = (char***)NULL;
	return;
}

bool CCSVFile::CheckCSV()
{
	bool quated = false;
	bool first = true;
	unsigned long rn = 0, cn = 0, dn = 0;
	int state = 0;	 // { 0, 1, 2 } = { 1byte����, 2byte�����̑�1byte, 2byte�����̑�2byte }
	int c;	// ����

	// �f�[�^�����ɓǂݍ���ł���ꍇ�A�������
	ReleaseBuffer();

	do
	{
		c = getc(fp);

		// �����X�e�[�g�����g�̕ύX
		if      ( ( state == 0 ) && ( CCSVFILE_is_sjis1( c ) ) )	state = 1; // 0 -> 1
		else if ( ( state == 1 ) && ( CCSVFILE_is_sjis2( c ) ) )	state = 2; // 1 -> 2
		else if ( ( state == 2 ) && ( CCSVFILE_is_sjis1( c ) ) )	state = 1; // 2 -> 1
		else										state = 0; // 2 -> 0, ���̑�

		// CSV�ł̃t�H�[�}�b�g��1byte�����Ɉ˂�
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

	// �I���k������'\0'�璷
	datumsize++;

	return true;
}

bool CCSVFile::ReadData()
{
	bool quated = false;
	bool first = true;
	unsigned long rn = 0, cn = 0, dn = 0;
	int state = 0;	 // { 0, 1, 2 } = { 1byte����, 2byte�����̑�1byte, 2byte�����̑�2byte }
	int c;	// ����

	do
	{
		c = getc(fp);
		if(c == EOF)
			break;

		// �����X�e�[�g�����g�̕ύX
		if      ( ( state == 0 ) && ( CCSVFILE_is_sjis1( c ) ) )	state = 1; // 0 -> 1
		else if ( ( state == 1 ) && ( CCSVFILE_is_sjis2( c ) ) )	state = 2; // 1 -> 2
		else if ( ( state == 2 ) && ( CCSVFILE_is_sjis1( c ) ) )	state = 1; // 2 -> 1
		else										state = 0;	// 2 -> 0, ���̑�

		// CSV�ł̃t�H�[�}�b�g��1byte�����Ɉ˂�
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
 
	// �o�b�t�@���Ȃ�
	if(CCSVFILE_NOBUFFER())
		return false;

	if( fopen_s( &fp , filename, "a") != 0 )
		return false;

	for(rn=0;rn<rowsize;rn++)
	{
		for(cn=0;cn<colsize;cn++)
		{
			// �ᔽ�������܂܂��`�F�b�N
			errChar = false;
			state = 0;
			for(dn=0;dn<datumsize && buffer[rn][cn][dn] != '\0';dn++)
			{
				// �����X�e�[�g�����g�̕ύX
				if      ( ( state == 0 ) && ( CCSVFILE_is_sjis1( buffer[rn][cn][dn] ) ) )	state = 1; // 0 -> 1
				else if ( ( state == 1 ) && ( CCSVFILE_is_sjis2( buffer[rn][cn][dn] ) ) )	state = 2; // 1 -> 2
				else if ( ( state == 2 ) && ( CCSVFILE_is_sjis1( buffer[rn][cn][dn] ) ) )	state = 1; // 2 -> 1
				else															state = 0; // 2 -> 0, ���̑�

				if(state == 0 && ( (tmp = buffer[rn][cn][dn]) == '\n' || tmp == '\r' || tmp == ',' || tmp == '"' ) )
				{
					errChar = true;
					break;
				}
			}

			// �ᔽ�������܂܂��˃N�H�[�g����
			errChar && putc('"', fp);

			state = 0;
			// �_�u���N�H�[�e�[�V�����͓�d�˂鏈���������Ȃ��當���񏑂��o��
			for(dn=0;dn<datumsize && buffer[rn][cn][dn] != '\0';dn++)
			{
				// �����X�e�[�g�����g�̕ύX
				if      ( ( state == 0 ) && ( CCSVFILE_is_sjis1( buffer[rn][cn][dn] ) ) )	state = 1; // 0 -> 1
				else if ( ( state == 1 ) && ( CCSVFILE_is_sjis2( buffer[rn][cn][dn] ) ) )	state = 2; // 1 -> 2
				else if ( ( state == 2 ) && ( CCSVFILE_is_sjis1( buffer[rn][cn][dn] ) ) )	state = 1; // 2 -> 1
				else															state = 0; // 2 -> 0, ���̑�

				if(state == 0 && buffer[rn][cn][dn] == '"')
					putc('"', fp);
				putc(buffer[rn][cn][dn], fp);
			}

			// �ᔽ�������܂܂��˃N�H�[�g����
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

	// ���Ƀf�[�^���������ꍇ�A�������
	ReleaseBuffer();

	// �e�ݒ�l�ݒ�
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

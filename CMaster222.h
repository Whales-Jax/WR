

#ifndef __CMASTER222_H__
#define __CMASTER222_H__

#include <iostream>
#include <sstream>
#include <fstream>

/*!
2008/06/08�@5108C015���
x�ɑΉ�����y�̒l�Cy�ɑΉ�����x�̒l��Ԃ��܂��D

�g����
�ȉ��̂悤��csv�t�@�C����p�ӂ��Ă��������D
//��������
0,5
1,2
2,4
3,7
4,2
5,8
6,1
end
//�����܂�
���ڂ�x�̒l�C���ڂ�y�̒l�Ƃ��ĔF������܂��D

�ȉ��g�����̗�
#include <iostream>
#include "CMaster222.h"
int main(void){
	CMaster222 cma;
	cma.setup("sample.csv");//�Z�b�g�A�b�v
	double y;
	y = cma.xtoy( 3.5 );//x��3.5�̂Ƃ���y�̒l���Ԃ����Dy�̒l�͔���ԂŌv�Z�����D
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
#ifndef _PropertyLithiumBromide_ver05_h_
#define _PropertyLithiumBromide_ver05_h_

//////////////////////////////////////////////////////////////////////////////////
//																				//
//  LiBr �����l�֐�(SI�P��)														//
//																				//
//------------------------------------------------------------------------------//
//	�֐���	   |����				|�P��		|���̓p�����[�^			|�o�W	//
//------------------------------------------------------------------------------//
//	sc_h_XT	   |�n�t�̔�G���^���s	|[kJ/kg]	|(�Z�x�C���x)			|*1		//
//		       (�g�p�͈�:40%�`70%,15���`165��)									//
//	sc_T_Xh	   |�n�t�̉��x			|[��]		|(�Z�x�C��G���^���s)	|*1		//
//		       (�g�p�͈�:40%�`70%,15���`165��)									//
//	sc_T_XTsat |�n�t�̉��x			|[��]		|(�Z�x�C�O�a���x)		|*1		//
//		       (�g�p�͈�:40%�`70%,5���`175��)									//
//	sc_X_TTsat |�n�t�̔Z�x	�@�@���A��@|		|(���x,�O�a���x)		|*1'	//
//																				//
//	sc_Tsat_XT |�n�t�̖O�a���x	���A��@|		|(�Z�x, ���x)			|*1'	//
//																				//
//	sc_rho_XT  |�n�t�̖��x			|[kg/m^3]	|(�Z�x�C���x)			|*4		//
//																				//
//	sc_X_Trho  |�n�t�̔Z�x ���A��@ |[kg/m^3]	|(���x�C���x)			|*4'	//
//																				//
//	sc_visc_XT |�n�t�̔S���W��		|[Pa�s]		|(�Z�x�C���x)			|*2		//
//																				//
//	sc_visc_XT2|�n�t�̔S���W��		|[Pa�s]		|(�Z�x�C���x)			|*3		//
//																				//
// 	sc_thc_XT  |�n�t�̔M�`����		|[W/m�K]	|(�Z�x�C���x)			|*6		//
//																				//
// 	sc_thc_XT2 |�n�t�̔M�`����		|[W/m�K]	|(�Z�x�C���x)			|*3		//
//																				//
//	sc_cp_XT   |�n�t�̒舳��M		|[kJ/kg�K]	|(�Z�x�C���x)			|*2		//
//																				//
//	sc_cp_XT2  |�n�t�̒舳��M		|[kJ/kg�K]	|(�Z�x�C���x)			|*3		//
//																				//
//	sc_st_XT   |�n�t�̕\�ʒ���		|[N/m]		|(�Z�x�C���x)			|*5-1	//
//																				//
//	sc_st_XT2  |�n�t�̕\�ʒ���		|[N/m]		|(�Z�x�C���x)			|*3		//	��������
//																				//
//	sc_d_X     |�n�t�̕����g�U�W��	|[m^2/s]	|(�Z�x)					|*3		//
//																				//
//	sc_d_XT    |�n�t�̕����g�U�W��	|[m^2/s]	|(�Z�x�C���x)			|*4		//
//																				//
//////////////////////////////////////////////////////////////////////////////////

// *1	2005 ASHRAE Handbook - Fundamentals (SI)
// *2	�V�ŁE��6��  �Ⓚ�󒲕֗��T��, (2010)
// *3	���{�����w���, �V�� �M�����n���h�u�b�N, �{����, (2008)  �o�W�L�q�Ȃ�
//		cp�@�g��E�A��, �Ⓚ35,397(1960) 815���Ǝv����D 
// *4	�����K�X�̎���
// *5	Patterson, M.R. and Perez-Blanco, H., "Numerical fits of the properties of lithium-bromide water solutions", ASHRAE annual meeting, June 1988.
// *5-1	�A��, �z�����Ⓚ�@�p��}-�z���܌n�̕���, �Ⓚ, Vol.52, No.600, 1975.
// *6	PROPERTIES OF LITHIUM BROMIDE-WATER SOLUTIONS AT HIGH TEMPERATURES AND CONDENSATIONS - PART I. Thermal Conductivity", ASHRAE Trans. 1990 Vol.96
//ASHRAE Trans. 1990 Vol.96 PROPERTIES OF LITHIUM BROMIDE-WATER SOLUTIONS AT HIGH TEMPERATURES AND CONDENSATIONS - PART II.  Density and Viscosity
//ASHRAE Trans. 1992 Vol.98 PROPERTIES OF LITHIUM BROMIDE-WATER SOLUTIONS AT HIGH TEMPERATURES AND CONDENSATIONS - PART III. Specific Heat


#include <cmath>

class PropertyLithiumBromide{
public:

	PropertyLithiumBromide();
	~PropertyLithiumBromide();

	double sc_h_XT    ( double _X , double _T    );
	double sc_T_Xh    ( double _X , double _h    );
	double sc_T_XTsat ( double _X , double _Tsat );
	double sc_X_TTsat ( double _T , double _Tsat );
	double sc_Tsat_XT ( double _X , double _T    );
	double sc_rho_XT  ( double _X , double _T    );
	double sc_X_Trho  ( double _T , double _rho  );
	double sc_visc_XT ( double _X , double _T    );
	double sc_visc_XT2( double _X , double _T    );
	double sc_thc_XT  ( double _X , double _T    );
	double sc_thc_XT2 ( double _X , double _T    );
	double sc_cp_XT   ( double _X , double _T    );
	double sc_cp_XT2  ( double _X , double _T    );
	double sc_st_XT   ( double _X , double _T    );
	double sc_st_XT2  ( double _X , double _T    );
	double sc_d_X	  ( double _X );
	double sc_d_XT    ( double _X , double _T    );

	double A[5];
	double B[5];
	double C[5];
	double D[4];
	double E[4];
	double F[2];
	double G[2];

	double An;
	double Bn;
	double Cn;
	double Dn;
	double En;

	double a;
	double b;

	double h;
	double T;
	double Tsat;
	double X;
	double rho;
	double thc;
	double visc;
	double cp;
	double st;
	double d;

	double t10;
	double a10;
	double b10;
	double c10;

	double dec;
	double eps;

	double E0;
	double E1;
	double E2;

};

#endif
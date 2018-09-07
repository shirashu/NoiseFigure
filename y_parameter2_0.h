/* --------------------------------------------------------------------
	Y�p�����[�^����ъe�함���萔�̌v�Z�ɗp�����\���̈ꗗ
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (y_parameter2_0.h)
-------------------------------------------------------------------- */

/* ----- �d���l�̍\���� -----*/
struct Current{
	double _[3][2];
		/*--�d���l [0][*]:Gate,[1][*]:Drain,[2][*]:Source*/
		/*--       [*][0]:Gate Step�d�����,[*][1]:Drain Step�d�����*/
};



/* ----- �A�h�~�^���X�̍\���� -----*/
struct Admittance{
	double g[2][2];
	double b[2][2];
};



/* ----- ������H�p�����[�^�̍\���� -----*/
struct CircuitPara{
	double cgd, cgs, cds;
	double gm0, gds;
	double ri;
	double tau;
	double ft, gm;
};

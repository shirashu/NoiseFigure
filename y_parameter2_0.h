/* --------------------------------------------------------------------
	Yパラメータおよび各種物理定数の計算に用いた構造体一覧
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (y_parameter2_0.h)
-------------------------------------------------------------------- */

/* ----- 電流値の構造体 -----*/
struct Current{
	double _[3][2];
		/*--電流値 [0][*]:Gate,[1][*]:Drain,[2][*]:Source*/
		/*--       [*][0]:Gate Step電圧印加,[*][1]:Drain Step電圧印加*/
};



/* ----- アドミタンスの構造体 -----*/
struct Admittance{
	double g[2][2];
	double b[2][2];
};



/* ----- 等価回路パラメータの構造体 -----*/
struct CircuitPara{
	double cgd, cgs, cds;
	double gm0, gds;
	double ri;
	double tau;
	double ft, gm;
};

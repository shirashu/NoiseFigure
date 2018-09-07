/* ----------------------------------------------------
	ノイズ指数用ヘッダファイル
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (noise_figure1_1.h)
------------------------------------------------------- */
void NoiseFigure(double sig[], double sid[], struct ComplexNumbers sigd[],
				 struct ComplexNumbers y11[], struct ComplexNumbers y12[],
				 struct ComplexNumbers y21[], struct ComplexNumbers y22[], 
				 double fr[], char *fp_name);

void MinimumNF_Gonzalez(double sig[], double sid[], struct ComplexNumbers sigd[],
				struct ComplexNumbers y11[], struct ComplexNumbers y12[],
				struct ComplexNumbers y21[], struct ComplexNumbers y22[], 
				double p[], double r[], double c[], double fmin[],
				double fr[], double gn[], struct ComplexNumbers zopt[]);

void MinimumNF_Cappy(double sig[], double sid[], struct ComplexNumbers sigd[],
				struct ComplexNumbers y11[], struct ComplexNumbers y12[],
				struct ComplexNumbers y21[], struct ComplexNumbers y22[], 
				double p[], double r[], double c[], double fmin[],
				double fr[], double gn[], struct ComplexNumbers zopt[]);

void NoiseFigureCircle(double fmin[], double fr[], double gn[], 
					   struct ComplexNumbers zopt[]);
void WriteNoise(double p[], double r[], double c[], double fmin[],
				double fr[], double gn[], struct ComplexNumbers zopt[], 
				char *fp_name);

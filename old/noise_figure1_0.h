/* ----------------------------------------------------
	ノイズ指数用ヘッダファイル
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (noise_figure1_0.h)
------------------------------------------------------- */
void NoiseFigure(double sig[], double sid[], struct ComplexNumbers sigd[],
				 struct ComplexNumbers y11[], struct ComplexNumbers y21[],
				 double fr[], char *fp_name);
void WriteNoise(double p[], double r[], double c[], double fmin[],
				double fr[], char *fp_name);

/* ----------------------------------------------------
	Yパラメータおよび各種物理定数の計算
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (y_parameter1_0.h)
------------------------------------------------------- */
void YParameter(double cur1[][JT][IN0], double ss[IN0][IN0], double rex[][IN0],
				struct ComplexNumbers y11[], struct ComplexNumbers y21[],
				double fr[], char *fn_yparam, char *fn_out);

void CalculateYPara_inte(double cur1[][JT][IN0], double g[][IN0][DATA],
						 double b[][IN0][DATA], double rex[][IN0],
						 double ss[][IN0], double dv[], double fr[]);
void FlatYPara_inte(double g[][IN0][DATA], double b[][IN0][DATA], double fr[]);

void CurrentFFT(double cur1[][JT][IN0],
				struct ComplexNumbers cur_f[][DATA][IN0]);
void ConverterToComp(struct ComplexNumbers temp[], double cur1[][JT][IN0],
					 int i, int ii);
void ConverterFromComp(struct ComplexNumbers temp[], 
					   struct ComplexNumbers cur_f[][DATA][IN0], int i, int ii);
void CalculateYPara_FFT(struct ComplexNumbers cur_f[][DATA][IN0],
						double g[][IN0][DATA], double b[][IN0][DATA],
						double rex[][IN0], double ss[][IN0],
						double dv[], double fr[]);

void CalculateOtherPara(double g[][IN0][DATA], double b[][IN0][DATA],
						double cgd[], double cgs[],	double ri[],
						double gm0[], double tau[], double cds[],
						double gds[], double fr[]);
void WriteParameter(double g[][IN0][DATA], double b[][IN0][DATA],
					double cgd[], double cgs[], double ri[],
					double gm0[], double tau[], double cds[],
					double gds[], double fr[], char *fn_yparam, char *fn_out);

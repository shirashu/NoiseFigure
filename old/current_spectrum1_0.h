/* ----------------------------------------------------
	スペクトル計算用ヘッダファイル
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (correlation_function1_0.h)
------------------------------------------------------- */
void CurrentSpectrum(double cur1[][JT][IN0], double corr[][NUMS],
					 double sig[], double sid[], struct ComplexNumbers sigd[],
					 double fr[], char *fp_name_spec);

void CalculateCorrelation_inte(double cur1[][JT][IN0], double corr[][NUMS]);
void CalculateSpectrum_inte(double corr[][NUMS], double sig[], double sid[],
							struct ComplexNumbers sigd[], double fr[]);
void WriteSpectrum_inte(double corr[][NUMS], double sig[], double sid[],
						struct ComplexNumbers sigd[], double fr[], char *fp_name);

void ConvertCompNumber(struct ComplexNumbers tg[],
					   struct ComplexNumbers td[], double cur1[][JT][IN0]);
void SpectrumWindows(struct ComplexNumbers t[]);
void CalculateSpectrum_FFT(double sig[], double sid[], struct ComplexNumbers sigd[],
						   struct ComplexNumbers tg[], struct ComplexNumbers td[]);
void FlatSpectrum_FFT(double sig[], double sid[], struct ComplexNumbers sigd[]);
void CalculateCorrelation_FFT(double sig[], double sid[], 
								struct ComplexNumbers sigd[], double corr[][NUMS]);
void WriteSpectrum_FFT(double corr[][NUMS], double sig[], double sid[], 
				   struct ComplexNumbers sigd[], char *fp_name);


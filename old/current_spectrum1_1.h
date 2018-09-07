/* ----------------------------------------------------
	スペクトル計算用ヘッダファイル
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (correlation_function1_1.h)
------------------------------------------------------- */
void CurrentSpectrum(double cur[][JTN0], double corr[][JTN0],
					 double sig[], double sid[], struct ComplexNumbers sigd[],
					 double fr[], char *fp_name_spec);

void LowPassFilter(double cur[][JTN0]);
void CalculateCorrelation_inte(double cur[][JTN0], double corr[][JTN0]);
void CalculateSpectrum_inte(double corr[][JTN0], double sig[], double sid[],
							struct ComplexNumbers sigd[], double fr[]);
void SpectrumWindow_inte(double sig[], double sid[], 
					   struct ComplexNumbers sigd[], double fr[]);
void SpectrumAppro_FFT(double sig[], double sid[], 
				  struct ComplexNumbers sigd[], double fr[]);
void WriteSpectrum_inte(double corr[][JTN0], double sig[], double sid[],
						struct ComplexNumbers sigd[], double fr[], char *fp_name);

void ConvertCompNumber(struct ComplexNumbers tg[],
					   struct ComplexNumbers td[], double cur[][JTN0]);
void TimeWindow(struct ComplexNumbers t[]);
void CalculateSpectrum_FFT(double sig[], double sid[], struct ComplexNumbers sigd[],
						   struct ComplexNumbers tg[], struct ComplexNumbers td[]);
void SpectrumWindow_FFT(double sig[], double sid[], 
					   struct ComplexNumbers sigd[], double fr[]);
void CalculateCorrelation_FFT(double sig[], double sid[], 
								struct ComplexNumbers sigd[], double corr[][JTN0]);
void WriteSpectrum_FFT(double corr[][JTN0], double sig[], double sid[], 
				   struct ComplexNumbers sigd[], double fr[], char *fp_name);


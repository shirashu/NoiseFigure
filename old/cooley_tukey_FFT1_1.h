/* ----------------------------------------------
	FFT関数用ヘッダファイル
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (Cooley_Tukey_FFT1_1.h)
---------------------------------------------- */
int FFTCalculate(struct ComplexNumbers sif[], int bit_num);
void BitRev(int bit_r[], int bit_num);
void BitRIndex(int bit_r[], struct ComplexNumbers x[], int bit_num);
void CooleyTukeyFFT(struct ComplexNumbers x[], int bit_num);

double PowReW(double y, int data_num);
double PowImW(double y, int data_num);

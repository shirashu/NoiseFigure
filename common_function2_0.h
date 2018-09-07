/* ---------------------------------------------------------------
	Main関数で使われる一般的な関数ヘッダファイル
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (common_function2_0.h)
----------------------------------------------------------------- */

void FrequencyIndex(double *fr, struct SetPara para);

FILE* OpenFile(char *argv, char mode);
int ReadLine(char *buffer, FILE *fp);

void Sigma(double x,double *sum, double *reg);
void TriangleWindow(double x[], int nmax, int avenum, int linenum);
void HanningWindow_f(double x[], int nmax);
void LowPassFilter(double x[], int nmax, double cutfreq, int fwin, struct SetPara para);

void StopPG();

int error(size_t st);


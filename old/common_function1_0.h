/* ---------------------------------------------------------------
	Main関数で使われる一般的な関数ヘッダファイル
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (common_function1_0.h)
----------------------------------------------------------------- */
void ReadCurrent(double cur[][JT][IN0], int ii, char *fp_name);
void FrequencyIndex(double fr[]);
void CurInit(double cur[][JT][IN0], double cur1[][JT][IN0], double ss[][IN0], 
			 double rex[][IN0], char *fp_name);

FILE* OpenFile(char *argv, char mode);
int ReadLine(char *buffer, FILE *fp);

void Sigma(double x,double *sum, double *reg);


/* ---------------------------------------------------------------
	Main関数で使われる一般的な関数ヘッダファイル
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (common_function1_1.h)
----------------------------------------------------------------- */
void FrequencyIndex(double fr[]);

void ReadCurrent(double cur[][JT][IN0], int ii, char *fp_name);
void CurInit(double cur[][JT][IN0], double cur1[][JT][IN0],
			 double ss[][IN0], double rex[][IN0]);

void ReadCurrentVC(double cur[][JTN0], char *fp_name);
void CurInitVC(double cur[][JTN0], char *fp_name);

FILE* OpenFile(char *argv, char mode);
int ReadLine(char *buffer, FILE *fp);

void Sigma(double x,double *sum, double *reg);


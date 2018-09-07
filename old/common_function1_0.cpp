/* ---------------------------------------------------------------
	Main関数で使われる一般的な関数
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (common_function1_0.cpp)
----------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_0.h" 
#include "common_function1_0.h" 

/* ************************************************ */
/* ----- FILEから電流値を読み込み、変数に代入 ----- */
void ReadCurrent(double cur[][JT][IN0], int ii, char *fp_name) {
	char buffer[128][64];
	int jt;
    FILE *fp;

    fp = OpenFile(fp_name, 'r');

	/* ----- 1行ずつ変数に代入 ----- */
	for(jt=0; jt<=JT-1; jt++) {
		ReadLine(&buffer[0][0], fp);
		sscanf(&buffer[0][0], " %lf %lf %lf ", &cur[1][jt][ii],
			&cur[0][jt][ii], &cur[2][jt][ii]);
	}
	fclose(fp);
	printf("'%s' read\n", fp_name);

}

/* ********************************************** */
/* ----- 求めたい周波数のインデックスを生成 ----- */
void FrequencyIndex(double fr[]){
	int jf;

	/* ------------------------------ */
	/* --- 数値積分法を用いた場合 --- */
	/* ------------------------------ */
	if(WAY==0){
		int fnum;
		double df, fmin, fmax;
		df   = 5.0e+9;		//周波数間隔
		fmin = -2.0e+10;	//最小周波数
		fmax = 1.2e+11;		//最大周波数
		fnum = (int)((fmax-fmin)/df); //周波数インデックスの個数

		for(jf=0; jf<=fnum; jf++){
			fr[jf] = df * double(jf) + fmin;
		}

		fr[fnum+1] = 32168; //配列の終わりを示す
	}

	/* ------------------------- */
	/* --- FFT法を用いた場合 --- */
	/* ------------------------- */
	if(WAY==1){
		for(jf=0; jf<=DATA/2-1; jf++){
			fr[jf] = jf/DT/DATA;
		}
		fr[DATA/2] = 32168; //配列の終わりを示す
	}
}


/* *********************************************** */
/* 電流のふらつき、定常状態の電流値を計算 ------- */
void CurInit(double cur[][JT][IN0], double cur1[][JT][IN0], double ss[][IN0], 
			 double rex[][IN0], char *fp_name){
	int i, ii, jt;
	double rs1, rs2, rr1, rr2;

    FILE *fp;
    fp = OpenFile(fp_name, 'w');
	fprintf(fp, "Ig-vg-bef	Id-vg-bef	Ig-vg-aft	Id-vg-aft	");
	fprintf(fp, "Ig-vd-bef	Id-vd-bef	Ig-vd-aft	Id-vd-aft\n");

	// --- 定常状態における電流値の計算
	for(ii=0; ii<=1; ii++){
		 ss[0][ii] = 0.0;	 ss[1][ii] = 0.0;
		rex[0][ii] = 0.0;	rex[1][ii] = 0.0;
		rs1 = 0.0;	rs2 = 0.0;
		rr1 = 0.0;	rr2 = 0.0;
		for(jt=JTP0-NUMSBEF; jt<=JTP0-1; jt++){
			Sigma(cur[0][jt][ii], &ss[0][ii], &rs1);
			Sigma(cur[1][jt][ii], &ss[1][ii], &rs2);
		}
		for(jt=JT-NUMS; jt<=JT-1; jt++){
			Sigma(cur[0][jt][ii], &rex[0][ii], &rr1);
			Sigma(cur[1][jt][ii], &rex[1][ii], &rr2);
		}
		 ss[0][ii] =  ss[0][ii]/double(NUMSBEF);
		 ss[1][ii] =  ss[1][ii]/double(NUMSBEF);
		rex[0][ii] = rex[0][ii]/double(NUMS);
		rex[1][ii] = rex[1][ii]/double(NUMS);
		fprintf(fp, "%lf	%lf	",  ss[0][ii],  ss[1][ii]);
		fprintf(fp, "%lf	%lf	", rex[0][ii], rex[1][ii]);
	}
	fprintf(fp, "\n");
	fclose(fp);


	/* -----------------------------------------------------*/
	/* --- 単純にIg=Is-Idにしたがって計算*/
	/* -----------------------------------------------------*/
	if(CUR0==1){
		// --- 電流のふらつきを計算
		for(i=0; i<=IN0-1; i++){
			for(ii=0; ii<=1; ii++){
				for(jt=0; jt<=JTL0; jt++){
					cur1[i][jt][ii] = cur[i][jt+JTP0][ii] - rex[i][ii];
				}
			}
		}
	}


	/* -----------------------------------------------------*/
	/* --- Ig=0、つまりIs=Idと仮定して計算*/
	/* -----------------------------------------------------*/
	if(CUR0==0){
		// --- 定常状態における電流値の計算
		for(ii=0; ii<=1; ii++){
			ss[1][ii] = 0.0;	rs1 = 0.0;
			for(jt=JTP0-NUMSBEF; jt<=JTP0-1; jt++){
				Sigma(cur[1][jt][ii]-cur[2][jt][ii], &ss[1][ii], &rs1);
			}
			ss[0][ii] = 0.0;
			ss[1][ii] =  ss[1][ii]/double(NUMSBEF)/2;
			ss[2][ii] = -ss[1][ii];
		}

		for(ii=0; ii<=1; ii++){
			rex[1][ii] = 0.0;	rex[0][ii] = 0.0;
			rr1 = 0.0;	rr2 = 0.0;
			for(jt=JT-NUMS; jt<=JT-1; jt++){
				Sigma(cur[0][jt][ii], &rex[0][ii], &rr1);
				Sigma(cur[1][jt][ii]-cur[2][jt][ii], &rex[1][ii], &rr2);
			}
			rex[0][ii] =  rex[0][ii]/double(NUMS);
			rex[1][ii] =  rex[1][ii]/double(NUMS)/2;
			rex[2][ii] = -rex[1][ii];
		}

		// --- 電流のふらつきを計算
		for(i=0; i<=IN0-1; i++){
			for(ii=0; ii<=1; ii++){
				for(jt=0; jt<=JTL0; jt++){
					cur1[i][jt][ii] = cur[i][jt+JTP0][ii] - rex[i][ii];
				}
			}
		}
		rex[0][0] = 0.0;
		rex[0][1] = 0.0;
	}


}

/* ****************************** */
/* ----- ファイルをオープン ----- */
FILE* OpenFile(char *argv, char mode) {
   FILE *fp;
   if((fp = fopen(argv, &mode)) == NULL) {
           fprintf(stderr, "Error: File Not Found (%s)\n", argv);
           return 0;
   }
   return fp;
}

/* ************************ */
/* ----- 一行読み込み ----- */
int ReadLine(char *buffer, FILE *fp) {
   int i=0;
   if(fgets(buffer, 256, fp) == NULL) {
           if(ferror(fp)) {
                   fprintf(stderr, "Error: Cannot Read File");
                   return 0;
           }
   if(feof(fp)) i = 1;
   }
   return i;
}

/* ********************************** */
/* ----- Siguma計算(数値計算用) ----- */
void Sigma(double x,double *sum, double *reg)
{
	double temp;
	*(reg)  = *(reg) + x;
	temp = *(sum);
	*sum  = *(sum) + *(reg);
	temp = *(sum) - temp;
	*(reg)  = *(reg) - temp;
}

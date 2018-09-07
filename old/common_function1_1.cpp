/* ---------------------------------------------------------------
	Main関数で使われる一般的な関数
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (common_function1_1.cpp)
----------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_1.h" 
#include "common_function1_1.h" 

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
		fmax = 1.5e+11;		//最大周波数
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
		double fmax;
		fmax = 1.50e+11;		//最大周波数
		for(jf=0; jf<=DATA/2-1; jf++){
			fr[jf] = jf/DT/DATA;
			if(fr[jf] > fmax){break;}
		}
		fr[jf+1] = 32168; //配列の終わりを示す
	}
}


/* **************************************************** */
/* ----- FILEから電流値を読み込み，等価回路特性用 ----- */
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

/* ************************************************************* */
/* 電流のふらつき、定常状態の電流値を計算、等価回路特性用------- */
void CurInit(double cur[][JT][IN0], double cur1[][JT][IN0],
			 double ss[][IN0], double rex[][IN0]){
	int i, ii, jt;
	double rs1, rs2, rr1, rr2;

	/* -----------------------------------------------------*/
	/* --- 単純にIg=Is-Idにしたがって計算*/
	/* -----------------------------------------------------*/
	if(CUR0==1){
		// --- 定常状態における電流値の計算
		for(ii=0; ii<=1; ii++){
			 ss[0][ii] = 0.0;	 ss[1][ii] = 0.0;
			rex[0][ii] = 0.0;	rex[1][ii] = 0.0;
			rs1 = 0.0;	rs2 = 0.0;
			rr1 = 0.0;	rr2 = 0.0;
			for(jt=JTP0-NUMP; jt<=JTP0-1; jt++){
				Sigma(cur[0][jt][ii], &ss[0][ii], &rs1);
				Sigma(cur[1][jt][ii], &ss[1][ii], &rs2);
			}
			for(jt=JT-NUML; jt<=JT-1; jt++){
				Sigma(cur[0][jt][ii], &rex[0][ii], &rr1);
				Sigma(cur[1][jt][ii], &rex[1][ii], &rr2);
			}
			 ss[0][ii] =  ss[0][ii]/double(NUMP);
			 ss[1][ii] =  ss[1][ii]/double(NUMP);
			rex[0][ii] = rex[0][ii]/double(NUML);
			rex[1][ii] = rex[1][ii]/double(NUML);
		}
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
			for(jt=JTP0-NUMP; jt<=JTP0-1; jt++){
				Sigma(cur[1][jt][ii]-cur[2][jt][ii], &ss[1][ii], &rs1);
			}
			ss[0][ii] = 0.0;
			ss[1][ii] =  ss[1][ii]/double(NUMP)/2;
			ss[2][ii] = -ss[1][ii];
		}

		for(ii=0; ii<=1; ii++){
			rex[1][ii] = 0.0;	rex[0][ii] = 0.0;
			rr1 = 0.0;	rr2 = 0.0;
			for(jt=JT-NUML; jt<=JT-1; jt++){
				Sigma(cur[0][jt][ii], &rex[0][ii], &rr1);
				Sigma(cur[1][jt][ii]-cur[2][jt][ii], &rex[1][ii], &rr2);
			}
			rex[0][ii] =  rex[0][ii]/double(NUML);
			rex[1][ii] =  rex[1][ii]/double(NUML)/2;
			rex[2][ii] = -rex[1][ii];
		}

		// --- 電流のふらつきを計算
		for(i=0; i<=IN0-1; i++){
			for(ii=0; ii<=1; ii++){
				for(jt=0; jt<=JTL0-1; jt++){
					cur1[i][jt][ii] = cur[i][jt+JTP0][ii] - rex[i][ii];
				}
			}
		}
		rex[0][0] = 0.0;
		rex[0][1] = 0.0;
	}

}
/* **************************************************** */
/* ----- FILEから電流値を読み込み、雑音特性解析用 ----- */
void ReadCurrentVC(double cur[][JTN0], char *fp_name) {
	char buffer[128][64];
	int jt;
    FILE *fp;

    fp = OpenFile(fp_name, 'r');

	for(jt=0; jt<=JTN0-1; jt++) {
		ReadLine(&buffer[0][0], fp);
		sscanf(&buffer[0][0], " %lf %lf %lf ", &cur[1][jt],
			&cur[0][jt], &cur[2][jt]);
	}
	fclose(fp);
	printf("'%s' read\n", fp_name);

}


/* ************************************************************* */
/* 電流のふらつき、定常状態の電流値を計算、雑音特性解析用------- */
void CurInitVC(double cur[][JTN0], char *fp_name){
	int jt;
	double rex[IN0], rr1, rr2, rr3;

    FILE *fp;
    fp = OpenFile(fp_name, 'w');
	fprintf(fp, "Ig	Id	Is\n");

	// --- 定常状態における電流値の計算
	rex[0] = 0.0;	rex[1] = 0.0;	rex[2] = 0.0;
	rr1 = 0.0;	rr2 = 0.0;	rr3 = 0.0;
	for(jt=JTN0-NUMN; jt<=JTN0-1; jt++){
		Sigma(cur[0][jt], &rex[0], &rr1);
		Sigma(cur[1][jt], &rex[1], &rr2);
		Sigma(cur[2][jt], &rex[2], &rr3);
	}
	rex[0] = rex[0]/double(NUMN);
	rex[1] = rex[1]/double(NUMN);
	rex[2] = rex[2]/double(NUMN);
	fprintf(fp, "%lf	%lf	%lf\n", rex[0], rex[1], rex[2]);
	fclose(fp);


	/* -----------------------------------------------------*/
	/* --- 単純にIg=Is-Idにしたがって計算*/
	/* -----------------------------------------------------*/
	if(CUR1==1){
		// --- 電流のふらつきを計算
		for(jt=0; jt<=JTN0-1; jt++){
			cur[0][jt] = cur[0][jt] - rex[0];
			cur[1][jt] = cur[1][jt] - rex[1];
			cur[2][jt] = cur[2][jt] - rex[2];
		}
	}

	/* -----------------------------------------------------*/
	/* --- Ig=0、つまりIs=Idと仮定して計算*/
	/* -----------------------------------------------------*/
	if(CUR1==0){
		// --- 定常状態における電流値の計算
		rex[0] = 0.0;	rex[1] = 0.0;
		rr1 = 0.0;	rr2 = 0.0;
		for(jt=JTN0-NUMN; jt<=JTN0-1; jt++){
			Sigma(cur[0][jt], &rex[0], &rr1);
			Sigma(cur[1][jt]-cur[2][jt], &rex[1], &rr2);
		}
		rex[0] =  rex[0]/double(NUMN);
		rex[1] =  rex[1]/double(NUMN)/2;
		rex[2] = -rex[1];

		// --- 電流のふらつきを計算
		for(jt=0; jt<=JTN0-1; jt++){
			cur[0][jt] = cur[0][jt] - rex[0];
			cur[1][jt] = cur[1][jt] - rex[1];
			cur[2][jt] = cur[2][jt] - rex[2];
		}
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

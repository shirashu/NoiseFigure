/* ---------------------------------------------------------------
	Main関数で使われる一般的な関数
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (common_function2_0.cpp)
----------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "new.h"

#include "parameter2_0.h" 
#include "common_function2_0.h" 


/* ********************************************** */
/* ----- 求めたい周波数のインデックスを生成 ----- */
void FrequencyIndex(double *fr, struct SetPara para){
	int jf;

	/* ------------------------------ */
	/* --- 数値積分法を用いた場合 --- */
	/* ------------------------------ */
	if(para.way==0){
		int fnum;
		fnum = (int)((para.fmax_y-para.fmin_y)/para.fstep_y); //周波数Indexの数

		for(jf=0;para.data/2-1 ; jf++){
			fr[jf] = para.fstep_y * double(jf) + para.fmin_y;
			if(fr[jf] > para.fmax_y){break;}
		}

//		fr[fnum+1] = 32168; // 配列の終わりを示す
	}

	/* ------------------------- */
	/* --- FFT法を用いた場合 --- */
	/* ------------------------- */
	else if(para.way==1){
		for(jf=0; jf<=para.data/2-1; jf++){
			fr[jf] = jf/para.dt/para.data;
//			if(fr[jf] > (para.fmax_y /*+ para.yave/para.dt/para.data)*/)){break;}
		}

//		fr[jf+1] = 32168; // 配列の終わりを示す
	}

	/* ------------------------------- */
	/* --- Yパラのみを計算する場合 --- */
	/* ------------------------------- */
	if(para.culc==1){
		for(jf=0;jf<=para.data/2-1; jf++){
			fr[jf] = jf/para.dt/para.data;
			if(fr[jf] > (para.fmax_y /*+ para.yave/para.dt/para.data)*/)){break;}
		}
//		fr[jf+1] = 32168; // 配列の終わりを示す
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
   if(fgets(buffer, 20000, fp) == NULL) {
           if(ferror(fp)) {
                   fprintf(stderr, "Error: Cannot Read File");
                   return 0;
           }
   if(feof(fp)) i = 1;
   }
   return i;
}


/* ********************************************** */
/* ----- 動的メモリ確保に失敗した場合の処理 ----- */
int error(size_t st){
	printf("メモリの確保に失敗しました。強制終了します。\n");
	exit(1);
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



/* **************************************************** */
/* ----- 三角形ウィンドウによるデータの平滑化処理 ----- */
void TriangleWindow(double x[], int nmax, int avenum, int linenum){
	int k, jf, *l;
	double *temp;
	double rs;

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	l = new int [linenum];
	temp = new double [linenum];

	for(jf=0; jf<=nmax; jf++){
		temp[jf] = 0;		rs = 0;

		l[jf] = avenum;
		if(jf-l[jf]+1<0)    {do{l[jf]--;}while(jf-l[jf]+1<0);}
		if(jf+l[jf]-1>nmax)    {do{l[jf]--;}while(jf+l[jf]-1>nmax);}

		for(k=-l[jf]+1; k<=l[jf]-1; k++){
			Sigma((l[jf]-abs(k))*x[jf-k], &temp[jf], &rs);
		}
	}
	for(jf=0; jf<=nmax; jf++){
		x[jf] = temp[jf]/l[jf]/l[jf];
	}
}
	
void HanningWindow_f(double x[], int nmax){

	int jf;
	double *y;
	y = new double [nmax+1];

	for(jf=0; jf<=nmax; jf++){
		if(jf==0){
			y[jf] = -0.25*x[nmax]+0.50*x[jf]-0.25*x[jf+1];
		}
		if(jf==nmax){
			y[jf] = -0.25*x[jf-1]+0.50*x[jf]-0.25*x[0];
		}
		if(jf!=0 && jf!=nmax){
			y[jf] = -0.25*x[jf-1]+0.50*x[jf]-0.25*x[jf+1];
		}
	}
			for(jf=0; jf<=nmax; jf++){
				x[jf] = y[jf];
			}
}


/* ***************************************** */
/* ----- 低域通過フィルタ ------------------ */
/* ------- 同時に4つの行に対して計算可能 --- */
void LowPassFilter(double x[], int nmax, double cutfreq, int fwin, struct SetPara para){

	/* ----- 変数宣言 */
	int jt, tau;
	double *temp0, *temp1, *temp2, *temp3, a, rr[4];

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	temp0 = new double [nmax];
	temp1 = new double [nmax];
	temp2 = new double [nmax];
	temp3 = new double [nmax];

	/* ----- フィルタの計算 */
	for(jt=0; jt<=nmax-1; jt++){
		temp0[jt]=0.0;		rr[0]=0.0;
		temp1[jt]=0.0;		rr[1]=0.0;
		temp2[jt]=0.0;		rr[2]=0.0;
		temp3[jt]=0.0;		rr[3]=0.0;
		for(tau=0; tau<=fwin; tau++){
			if(jt-tau<0){;}
			else if(jt-tau>nmax-1){;}
			else {
				if(tau==0){ a = 2*cutfreq*para.dt;}
				else{
					a = sin(2*para.pai*cutfreq*para.dt*tau)/(2*para.pai*tau)*
						(1+cos(para.pai*tau/fwin));
				}
				Sigma(x[jt-tau]*a, &temp0[jt], &rr[0]);
				Sigma(x[nmax*1 + jt-tau]*a, &temp1[jt], &rr[1]);
				Sigma(x[nmax*2 + jt-tau]*a, &temp2[jt], &rr[2]);
				Sigma(x[nmax*3 + jt-tau]*a, &temp3[jt], &rr[3]);
			}
		}
	}

	for(jt=0; jt<=nmax-1; jt++){
		x[jt] = temp0[jt];
		x[nmax*1 + jt] = temp1[jt];
		x[nmax*2 + jt] = temp2[jt];
		x[nmax*3 + jt] = temp3[jt];
	}


	/* ----- メモリの開放 */
	delete [] temp0;
	delete [] temp1;
	delete [] temp2;
	delete [] temp3;

}



/* **************************************** */
/* ----- プログラムをいったん停止する ----- */
void StopPG(){
	char a;
	printf("Press any key to continue.");
	scanf("%c", &a);
}

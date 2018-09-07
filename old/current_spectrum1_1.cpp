/* ----------------------------------------------------
	スペクトル計算用ヘッダファイル
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (correlation_function1_1.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_1.h" 
#include "common_function1_1.h" 
#include "current_spectrum1_1.h" 
#include "cooley_tukey_FFT1_1.h"


/* ***************************+************** */
/* ----- 電流のスペクトル計算のMian関数 ----- */
void CurrentSpectrum(double cur[][JTN0], double corr[][JTN0],
					 double sig[], double sid[], struct ComplexNumbers sigd[],
					 double fr[], char *fp_name_spec){

	struct ComplexNumbers tg[NUMN], td[NUMN];

	// 低域通過フィルタ
	if(LOWPASS==1){
		LowPassFilter(cur);
	}

	/* -----------------------------------------------------*/
	/* --- 相関関数からフーリエ変換によりスペクトルを計算   */
	/* -----------------------------------------------------*/
	if(WAY==0){
		// --- 電流のふらつきから相関関数を計算
		CalculateCorrelation_inte(cur, corr);
		// --- 相関関数からスペクトルを計算
		CalculateSpectrum_inte(corr, sig, sid, sigd, fr);
		// --- スペクトルを三角波ウィンドウにより補正
		if(SAVE!=0){
			SpectrumWindow_inte(sig, sid, sigd, fr);
		}
		// ----- FILEへスペクトルと相関関数を書き込む
		WriteSpectrum_inte(corr, sig, sid, sigd, fr, fp_name_spec);
	}

	/* -----------------------------------------------------*/
	/* --- FFTを用いたスペクトル計算*/
	/* -----------------------------------------------------*/
	if(WAY==1){
		// --- ふらつき電流を複素数の値に変換
		ConvertCompNumber(tg, td, cur);
		// --- タイムウィンドウをかける
		if(TWIN==1){
			TimeWindow(tg);
			TimeWindow(td);
		}
		// --- 電流ふらつきをFFTにかける
		FFTCalculate(tg, JBIT);
		FFTCalculate(td, JBIT);
		// --- スペクトルを計算する
		CalculateSpectrum_FFT(sig, sid, sigd, tg, td);
		// --- スペクトルウィンドウにより補正
		if(SAVE!=0){
			SpectrumWindow_FFT(sig, sid, sigd, fr);
		}
		// --- スペクトルから相関関数を計算する
		CalculateCorrelation_FFT(sig, sid, sigd, corr);
		// ----- FILEへスペクトルと相関関数を書き込む
		WriteSpectrum_FFT(corr, sig, sid, sigd, fr, fp_name_spec);
	}

}


/* ********************************************************************** 
/* -------------------------------------------------------------------------
  以下各種関数
------------------------------------------------------------------------- */


/* ****************************************************** */
/* 数値積分法(相関関数を求めて計算) --------------------- */
/* ****************************************************** */

/* ****************************************************** */
/* ----- 低域通過フィルタ ----- */
void LowPassFilter(double cur[][JTN0]){
	int jt, tau;
	double temp[IN0][JTN0], rr0, rr1, rr2;
	double a;

	for(jt=0; jt<=JTN0-1; jt++){
		temp[0][jt]=0.0;	temp[1][jt]=0.0;	temp[2][jt]=0.0;
		rr0=0.0;			rr1=0.0;			rr2=0.0;
		for(tau=-FWIN; tau<=FWIN; tau++){
			if(tau==0){ a = 2*CUTFREQ*DT;}
			else{
				a = sin(2*PAI*CUTFREQ*DT*tau)/(PAI*tau)/2*(1+cos(PAI*tau/(FWIN+1)));
			}
			if(jt-tau<0){a=0;}
			else if(jt-tau>JTN0-1){a=0;}

			Sigma(cur[0][jt-tau]*a, &temp[0][jt], &rr0);
			Sigma(cur[1][jt-tau]*a, &temp[1][jt], &rr1);
			Sigma(cur[2][jt-tau]*a, &temp[2][jt], &rr2);
		}
	}

	for(jt=0; jt<=JTN0-1; jt++){
		cur[0][jt] = temp[0][jt];
		cur[1][jt] = temp[1][jt];
		cur[2][jt] = temp[2][jt];
	}

}

/* ****************************************************** */
/* ----- 電流のふらつきから直接に相関関数を計算する ----- */
void CalculateCorrelation_inte(double cur[][JTN0], double corr[][JTN0]){
	int dt, jt;
	double temp1, temp2, temp3, temp4, rs1, rs2, rs3, rs4;

	for(jt=0; jt<=MSMALLS; jt++){
		temp1 = 0;	temp2 = 0;		rs1 = 0;	rs2 = 0;
		for(dt=(JTN0-NUMN)+1; dt<=(JTN0-NUMN)+MLARGE-MSMALLS; dt++){
			Sigma(cur[0][dt]*cur[0][dt+jt], &temp1, &rs1);//自己相関
			Sigma(cur[1][dt]*cur[1][dt+jt], &temp2, &rs2);//自己相関
		}
		corr[0][jt] = temp1/(double)(MLARGE-MSMALLS)/EPP;
		corr[1][jt] = temp2/(double)(MLARGE-MSMALLS)/EPP;
	}

	for(jt=0; jt<=MSMALLC; jt++){
		temp3=0;	temp4=0;	rs3 = 0;	rs4 = 0;
		for(dt=(JTN0-NUMN)+1; dt<=(JTN0-NUMN)+MLARGE-MSMALLC; dt++){
			Sigma(cur[0][dt]*cur[1][dt+jt], &temp3, &rs3);//相互相関
			Sigma(cur[1][dt]*cur[0][dt+jt], &temp4, &rs4);//相互相関
		}
		corr[2][jt] = temp3/(double)(MLARGE-MSMALLC)/EPP;
		corr[3][jt] = temp4/(double)(MLARGE-MSMALLC)/EPP;
	}
}

/* ********************************************************** */
/* --- 相関関数からフーリエ変換によりスペクトルを計算する --- */
void CalculateSpectrum_inte(double corr[][JTN0], double sig[], double sid[],
							struct ComplexNumbers sigd[], double fr[]){
	int jt, fmax, jf, tau;
	double tg[NUMN], td[NUMN];
	struct ComplexNumbers tgd[NUMN];
	double rs1, rs2;
	double cmax, cgd, cdg, a[NUMN], b[NUMN];


	/* ---------------------------------------*/
	/* --- 自己相関関数からのスペクトル計算   */
	/* ---------------------------------------*/
	// --- 相関関数からのスペクトル計算
	for(jf=0; ; jf++){
		if(fr[jf]==32168){fmax=jf-1; break;} //最大周波数で計算をストップ
		rs1 = 0;	rs2 = 0;
		tg[jf] = corr[0][0] + corr[0][MSMALLS]*cos(2*PAI*fr[jf]*MSMALLS*DT);
		td[jf] = corr[1][0] + corr[1][MSMALLS]*cos(2*PAI*fr[jf]*MSMALLS*DT);
		for(jt=1; jt<=MSMALLS-1; jt++){
			Sigma(2*corr[0][jt]*cos(2*PAI*fr[jf]*jt*DT), &tg[jf], &rs1);
			Sigma(2*corr[1][jt]*cos(2*PAI*fr[jf]*jt*DT), &td[jf], &rs2);
		}
		tg[jf] = tg[jf]*DT;
		td[jf] = td[jf]*DT;
	}

	for(jf=0; jf<=fmax; jf++){
		sig[jf] = tg[jf];
		sid[jf] = td[jf];
	}


	/* -----------------------------------------------------*/
	/* --- 相互相関関数からのスペクトル計算   */
	/* -----------------------------------------------------*/
	cmax = corr[2][0];	tau = 0;
	for(jt=1; jt<=MSMALLC; jt++){
		if(corr[2][jt]>cmax){tau=jt; cmax=corr[2][jt];}
	}

	for(jt=0; jt<=MSMALLC-tau; jt++){
		cgd=corr[2][jt+tau];
		if(jt-tau>=0){cdg=corr[3][jt-tau];}
		else {cdg=corr[2][-jt+tau];}
		a[jt] = (cgd + cdg)/2;
		b[jt] = (cgd - cdg)/2;
	}

	// --- 相関関数からのスペクトル計算
	for(jf=0; jf<=fmax; jf++){
		rs1 = 0;	rs2 = 0;
		tgd[jf].re = a[0] + a[MSMALLC-tau]*cos(2*PAI*fr[jf]*(MSMALLC-tau)*DT);
		tgd[jf].im =        b[MSMALLC-tau]*sin(2*PAI*fr[jf]*(MSMALLC-tau)*DT);
		for(jt=1; jt<=MSMALLC-tau-1; jt++){
			Sigma(2*a[jt]*cos(2*PAI*fr[jf]*jt*DT), &tgd[jf].re, &rs1);
			Sigma(2*b[jt]*sin(2*PAI*fr[jf]*jt*DT), &tgd[jf].im, &rs2);
		}
		tgd[jf].re = tgd[jf].re*DT;
		tgd[jf].im = tgd[jf].im*DT;
	}

	for(jf=0; jf<=fmax; jf++){
		a[jf] = tgd[jf].re;
		b[jf] = tgd[jf].im;
	}

	for(jf=0; jf<=fmax-1; jf++){
		sigd[jf].re = (a[jf]*cos(2*PAI*fr[jf]*tau*DT) -
					   b[jf]*sin(2*PAI*fr[jf]*tau*DT));
		sigd[jf].im = (a[jf]*sin(2*PAI*fr[jf]*tau*DT) +
					   b[jf]*cos(2*PAI*fr[jf]*tau*DT));
	}
}


/* ************************************************************ */
/* ----- 三角形ウィンドによるスペクトルの周波数平滑を行う ----- */
void SpectrumWindow_inte(double sig[], double sid[], 
				  struct ComplexNumbers sigd[], double fr[]){
	int k, jf, fmax, fmin, l1[NUMN], l2[NUMN];
	double t1[NUMN], t2[NUMN], t3[NUMN], t4[NUMN];
	double rs1, rs2, rs3, rs4;

	for(jf=0; ; jf++){
		if(fr[jf]==0){fmin=jf;} //周波数0Hzのインデックスを記録
		else if(fr[jf]==32168){fmax=jf-1; break;} //最大周波数で計算をストップ
	}

	// --- Sig，Sigdを求める
	for(jf=0; jf<=fmax; jf++){
		t1[jf] = 0; t2[jf] = 0; t3[jf] = 0; t4[jf] = 0;
		rs1 = 0;	rs2 = 0;	rs3 = 0;	rs4 = 0;

		l1[jf] = SAVE; l2[jf] = SAVE;
		if(jf-l1[jf]+1<0)    {do{l1[jf]--;}while(jf-l1[jf]+1<0);}
		if(jf+l2[jf]-1>fmax)    {do{l2[jf]--;}while(jf+l2[jf]-1>fmax);}

		for(k=-l2[jf]+1; k<=l1[jf]-1; k++){
			Sigma((SAVE-abs(k))*sig[jf-k], &t1[jf], &rs1);
			Sigma((SAVE-abs(k))*sigd[jf-k].re, &t3[jf], &rs3);
			Sigma((SAVE-abs(k))*sigd[jf-k].im, &t4[jf], &rs4);
		}
	}

	for(jf=0; jf<=fmax; jf++){
		sig[jf] = t1[jf]/SAVE/(l1[jf]+l2[jf])*2;
		sigd[jf].re = t3[jf]/SAVE/(l1[jf]+l2[jf])*2;
		sigd[jf].im = t4[jf]/SAVE/(l1[jf]+l2[jf])*2;
	}

	// --- Sidを求める
	for(jf=1; jf<=fmax; jf++){
		t2[jf] = 0;	rs2 = 0;
		l1[jf] = SAVE; l2[jf] = SAVE;
		if(jf-l1[jf]+1<1)    {do{l1[jf]--;}while(jf-l1[jf]+1<1);}
		if(jf+l2[jf]-1>fmax)    {do{l2[jf]--;}while(jf+l2[jf]-1>fmax);}

		for(k=-l2[jf]+1; k<=l1[jf]-1; k++){
			Sigma((SAVE-abs(k))*sid[jf-k], &t2[jf], &rs2);
		}
	}

	for(jf=1; jf<=fmax; jf++){
		sid[jf] = t2[jf]/SAVE/(l1[jf]+l2[jf])*2;
	}


	// ---  周波数0Hzでスペクトルsig，sigdを0として補正
	for(jf=fmax; jf>=0; jf--){
		sig[jf] = sig[jf] - sig[fmin];
		sigd[jf].im = sigd[jf].im - sigd[fmin].im;
		sigd[jf].re = sigd[jf].re - sigd[fmin].re;

	}

}

/* ************************************************ */
/* ----- FILEへスペクトルと相関関数を書き込む ----- */
void WriteSpectrum_inte(double corr[][JTN0], double sig[], double sid[],
						struct ComplexNumbers sigd[], double fr[], char *fp_name){
	int k, fmax, fmin, ms;
	char fn_spe[]=SPE;
    FILE *fp;

	// --- ファイルオープン
	fp = OpenFile(fn_spe, 'a');

	fprintf(fp, "ΔI_G	ΔI_D	");
	fprintf(fp, "NUMN	MSMALL_Self	MSMALL_Cross	MLARGE	DATA\n");
	fprintf(fp, "%lf	%lf	",sqrt(corr[0][0]), sqrt(corr[1][0]));
	fprintf(fp, "%d	%d	%d	%d	%d\n", NUMN, MSMALLS, MSMALLC, MLARGE, DATA);

	fprintf(fp, "t[sec]	Cor_G(t)	Cor_D(t)	Cor_GD(t)	");
	fprintf(fp, "f[Hz]	Si_G(f)	Si_D(f)	Si_gd(f)[re]	Si_gd(f)[im]\n");

	for(k=0; ; k++){
		if(fr[k]==0){fmin=k; continue;}
		else if(fr[k]>1.05e+11){fmax=k-1; break;} //最大周波数で計算をストップ
	}

	for(k=0; k<=fmax-fmin; k++){
		fprintf(fp, "%1.4e	", k*DT);
		fprintf(fp, "%1.4e	%1.4e	%1.4e	", corr[0][k], corr[1][k], corr[2][k]);
		fprintf(fp, "%1.4e	%1.4e	%1.4e	%1.4e	%1.4e\n", fr[k+fmin] , sig[k+fmin], sid[k+fmin], sigd[k+fmin].re, sigd[k+fmin].im);
	}

	if(MSMALLS > MSMALLC){ ms = MSMALLC;}
	else { ms = MSMALLS;}
	for(k=fmax-fmin+1; k<=ms; k++){
		fprintf(fp, "%1.4e	", k*DT);
		fprintf(fp, "%1.4e	%1.4e	%1.4e\n", corr[0][k], corr[1][k], corr[2][k]);
	}

	// --- ファイルクローズ
	fclose(fp);
	printf("'%s' write\n", fn_spe);
}








/* ****************************************************** */
/* FFT法(スペクトルを直接計算) -------------- */
/* ****************************************************** */
/* ****************************************** */
/* ----- ふらつき電流を複素数の値に変換 ----- */
void ConvertCompNumber(struct ComplexNumbers tg[],
					   struct ComplexNumbers td[], double cur[][JTN0]){
	int k;

	for(k=0; k<=NUMN-1; k++){
		tg[k].re = cur[0][JTN0-NUMN+k];
		tg[k].im = 0.0;
		td[k].re = cur[1][JTN0-NUMN+k];
		td[k].im = 0.0;
	}
}

/* **************************************** */
/* ----- タイムウィンドウを掛ける ----- */
void TimeWindow(struct ComplexNumbers t[]){
	int i;
	// 0～T/10 までのウィンド
	for(i=0; i<=(int)(DATA/10); i++){
		t[i].re = t[i].re * (1-cos(PAI*10*i/DATA))/2;
	}

	// 9T/10～T までのウィンド
	for(i=(int)(DATA*9/10); i<=DATA-1; i++){
		t[i].re = t[i].re * (1+cos(PAI*10*((i-(int)(DATA*9/10))/DATA)))/2;
	}

}

/* ******************************** */
/* ----- スペクトルを計算する ----- */
void CalculateSpectrum_FFT(double sig[], double sid[], struct ComplexNumbers sigd[],
						   struct ComplexNumbers tg[], struct ComplexNumbers td[]){
	int jf;
	double st;

	if(TWIN==0){st = 1.0;}
	else if(TWIN==1){st = 0.875;}//タイムウィンドウを掛けた際の強度調整

	for(jf=0; jf<=DATA/2-1; jf++){
		sig[jf]    = (tg[jf].re*tg[jf].re + tg[jf].im*tg[jf].im)/st/EPP;
		sid[jf]    = (td[jf].re*td[jf].re + td[jf].im*td[jf].im)/st/EPP;
		sigd[jf].re= (tg[jf].re*td[jf].re + tg[jf].im*td[jf].im)/st/EPP;
		sigd[jf].im=-(tg[jf].im*td[jf].re - tg[jf].re*td[jf].im)/st/EPP;
	}
}

/* ************************************************************ */
/* ----- スペクトルウィンドによる周波数平滑を行う ----- */
void SpectrumWindow_FFT(double sig[], double sid[], 
				  struct ComplexNumbers sigd[], double fr[]){
	int k, jf, fmax, l1[NUMN], l2[NUMN];
	double t1[NUMN], t2[NUMN], t3[NUMN], t4[NUMN];
	double rs1, rs2, rs3, rs4;

	// 最大周波数の設定
	fmax = DATA/2-1;

/* -------------------------------------------*/
/* ----- 三角形ウィンドによる周波数平滑 ----- */
	for(jf=0; jf<=fmax; jf++){
		t1[jf] = 0; t2[jf] = 0; t3[jf] = 0; t4[jf] = 0;
		rs1 = 0;	rs2 = 0;	rs3 = 0;	rs4 = 0;

		l1[jf] = SAVE; l2[jf] = SAVE;
		if(jf-l1[jf]+1<0)    {do{l1[jf]--;}while(jf-l1[jf]+1<0);}
		if(jf+l2[jf]-1>fmax)    {do{l2[jf]--;}while(jf+l2[jf]-1>fmax);}

		for(k=-l2[jf]+1; k<=l1[jf]-1; k++){
			Sigma((SAVE-abs(k))*sig[jf-k], &t1[jf], &rs1);
			Sigma((SAVE-abs(k))*sid[jf-k], &t2[jf], &rs2);
			Sigma((SAVE-abs(k))*sigd[jf-k].re, &t3[jf], &rs3);
			Sigma((SAVE-abs(k))*sigd[jf-k].im, &t4[jf], &rs4);
		}
	}
	for(jf=0; jf<=fmax; jf++){
		sig[jf] = t1[jf]/SAVE/(l1[jf]+l2[jf])*2 *DT/DATA;
		sid[jf] = t2[jf]/SAVE/(l1[jf]+l2[jf])*2 *DT/DATA;
		sigd[jf].re = t3[jf]/SAVE/(l1[jf]+l2[jf])*2 *DT/DATA;
		sigd[jf].im = t4[jf]/SAVE/(l1[jf]+l2[jf])*2 *DT/DATA;
	}


/* -------------------------------------------------------------------*/
/* ----- 周波数0Hzでスペクトルsig，sigdを0として2次曲線にて近似 ----- */
	if(FETAPP==1){
		int l, fmin_ave_index, fmax_ave_index;
		double a, rs, fmin_ave, fmax_ave;
		fmin_ave = FMINAVE;
		fmax_ave = FMAXAVE;

		// 近似を行うための周波数帯域の周波数インデックスを設定
		for(jf=0; jf<=fmax; jf++){
			if(fr[jf] >= fmin_ave){fmin_ave_index=jf;	break;}
		}
		l=0;
		for(jf=fmin_ave_index; jf<=fmax; jf++){
			if(fr[jf] >= fmax_ave){fmax_ave_index=jf-1;	break;}
			l = l+1;
		}

		/* --- Sid : 一定値に設定
		a = 0;	rs = 0;
		for(jf=fmin_ave_index; jf<=fmax_ave_index; jf++){
			Sigma(sid[jf], &a, &rs);
		}
		a = a / ((float)l);
		for(jf=0; jf<=fmax; jf++){
			sid[jf] = a;
			if(fr[jf]>1.05e+11){break;}
		}

		/* --- Sig : 周波数の2乗に比例するように設定 */
		a = 0;	rs = 0;
		for(jf=fmin_ave_index; jf<=fmax_ave_index; jf++){
			Sigma(sig[jf]/fr[jf]/fr[jf], &a, &rs);
		}
		a = a / ((float)l);
		for(jf=0; jf<=fmax; jf++){
			sig[jf] = a * fr[jf]*fr[jf];
			if(fr[jf]>1.05e+11){break;}
		}

		/* --- Im[Sigd] : 周波数に比例するように設定 */
		a = 0;	rs = 0;
		for(jf=fmin_ave_index; jf<=fmax_ave_index; jf++){
			Sigma(sigd[jf].im/fr[jf], &a, &rs);
		}
		a = a / ((float)l);
		for(jf=0; jf<=fmax; jf++){
			sigd[jf].im = a * fr[jf];
			if(fr[jf]>1.05e+11){break;}
		}
	}

}

/* ******************************************** */
/* ----- スペクトルから相関関数を計算する ----- */
void CalculateCorrelation_FFT(double sig[], double sid[], 
								struct ComplexNumbers sigd[], double corr[][JTN0]){
	int i;
	struct ComplexNumbers temp1[DATA], temp2[DATA], temp3[DATA], temp4[DATA];

	// スペクトルを一時変数配列に格納
	for(i=0; i<=DATA/2-1; i++){
		temp1[i].re= sig[i]; //自己相関関数 sig
		temp1[i].im= 0;
		temp2[i].re= sid[i]; //自己相関関数 sig
		temp2[i].im= 0;
		temp3[i].re= sigd[i].re; //相互相関関数 sigd
		temp3[i].im= 0;
		temp4[i].re= sigd[i].im; //相互相関関数 sigd
		temp4[i].im= 0;
	}

	// スペクトルをFFTにかける
	FFTCalculate(&temp1[0], JBIT-1);
	FFTCalculate(&temp2[0], JBIT-1);
	FFTCalculate(&temp3[0], JBIT-1);
	FFTCalculate(&temp4[0], JBIT-1);

	// スペクトルのFFTから相関関数を計算する
	for(i=0; i<=DATA/4-1; i++) {
		corr[0][i] = 2*temp1[i].re/DATA/DT;
		corr[1][i] = 2*temp2[i].re/DATA/DT;
		corr[2][i] = 2*(temp3[i].re-temp4[i].im)/DATA/DT;
		corr[3][i] = 2*(temp3[i].re+temp4[i].im)/DATA/DT;
	}

}

/* ************************************************ */
/* ----- FILEへスペクトルと相関関数を書き込む ----- */
void WriteSpectrum_FFT(double corr[][JTN0], double sig[], double sid[], 
				   struct ComplexNumbers sigd[], double fr[], char *fp_name) {
	int k, fmax;
    FILE *fp;

	fp = OpenFile(fp_name, 'a');

	fprintf(fp, "ΔI_G	ΔI_D	JTN0	DATA\n");
	fprintf(fp, "%lf	%lf	%d	%d\n",sqrt(corr[0][0]), sqrt(corr[1][0]), JTN0, DATA);

	fprintf(fp, "t[sec]	Cor_G(t)	Cor_D(t)	CrossCor(t)	");
	fprintf(fp, "f[Hz]	Si_G(f)	Si_D(f)	Si_G-D(f).re	Si_G-D(f).im");

	fprintf(fp, "\n");


	for(k=0; ; k++){
		if(fr[k]==32168){fmax=k-1; break;}
	}

	for(k=0; k<=DATA/4-1; k++){
		if(k <= fmax){
//		if(k <= 60000){
			// --- 相関関数の書き込み
			fprintf(fp, "%1.4e	", k*DT);
			fprintf(fp, "%1.4e	%1.4e	%1.4e	", corr[0][k], corr[1][k], corr[2][k]);
			// --- スペクトルの書き込み
			fprintf(fp, "%1.4e	", k/DT/DATA);
			fprintf(fp, "%1.4e	" , sig[k]);
			fprintf(fp, "%1.4e	" , sid[k]);
			fprintf(fp, "%1.4e	%1.4e", sigd[k].re, sigd[k].im);
			fprintf(fp, "\n");
		}
	}

	fclose(fp);
	printf("'%s' write\n", fp_name);
}


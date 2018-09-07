/* ----------------------------------------------------
	スペクトル計算 --- 実際の数値計算関数
	（FFTを用いたスペクトル計算）
	2015/ T.takahashi 
	Ver. 3.0 (spe_FFT2_0.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "new.h"

#include "parameter2_0.h" 
#include "common_function2_0.h" 
#include "current_spectrum2_0.h" 

void FFTCalculate2(struct ComplexNumbers t[], struct SetPara para);
void FFTCalculate3(struct ComplexNumbers t[], struct SetPara para);

/* ***************************+********************** */
/* ----- FFTによるスペクトル計算のMian関数 ----- */
void Spectrum_FFT(struct Current cur[], struct Correlation corr[],
				  struct Spectrum si[], double fr[], struct SetPara para){


	/* ================================================================= */
	/*** --- 変数宣言 --- ***/
	struct ComplexNumbers *tg, *td;

	/*** --- プロトタイプ宣言 --- ***/
	void SpectrumAppro_FFT(struct Spectrum si[], double fr[], 
						   struct SetPara para);
	void ConvertCompNumber(struct ComplexNumbers tg[],struct ComplexNumbers td[], 
						   struct Current cur[], struct SetPara para);
	void TimeWindow(struct ComplexNumbers t[], struct SetPara para);
	void HanningWindow_f(double x[], int nmax);
	void CalculateSpectrum_FFT(struct Spectrum si[], struct ComplexNumbers tg[], 
							   struct ComplexNumbers td[], struct SetPara para);
	void SpectrumWindow_FFT(struct Spectrum si[], 
							struct ComplexNumbers tg[],struct ComplexNumbers td[], 
							double fr[], struct SetPara para);
	void CalculateCorrelation_FFT(struct Spectrum si[], struct Correlation corr[],
								  struct SetPara para);
	void WriteSpectrum_FFT(struct Correlation corr[], struct Spectrum si[], 
						   struct ComplexNumbers tg[],struct ComplexNumbers td[], 
						   double fr[], struct SetPara para);

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	tg = new ComplexNumbers[para.data];
	td = new ComplexNumbers[para.data];

	/* ================================================================= */

	/* --- ふらつき電流を複素数の値に変換 */
	ConvertCompNumber(tg, td, cur, para);

	/* --- タイムウィンドウをかける */
	if(para.twin==1){
		TimeWindow(tg, para);
		TimeWindow(td, para);
	}
	/* --- 電流ふらつきをFFTにかける */
	FFTCalculate2(td, para);
	FFTCalculate2(tg, para);

	/* --- スペクトルを計算する */
	CalculateSpectrum_FFT(si, tg, td, para);

	/* --- スペクトルウィンドウにより補正 */
	SpectrumWindow_FFT(si, tg, td, fr, para);

	/* --- スペクトルから相関関数を計算する */
	if(para.cor_out == 1){
		CalculateCorrelation_FFT(si, corr, para);
	}

	/* --- [Sig ∝ f^2]　[Sid = Const]　[Sigid ∝ f] の計算 */
	SpectrumAppro_FFT(si, fr, para);

	/* ----- FILEへスペクトルと相関関数を書き込む */
	WriteSpectrum_FFT(corr, si, tg, td, fr, para);

	/* ----- メモリの開放 */
	delete [] tg;
	delete [] td;

}


/* ********************************************************************** 
/* -------------------------------------------------------------------------
  以下各種関数
------------------------------------------------------------------------- */

/* ****************************************** */
/* ----- ふらつき電流を複素数の値に変換 ----- */
void ConvertCompNumber(struct ComplexNumbers tg[], struct ComplexNumbers td[], 
					   struct Current cur[], struct SetPara para){
	int k;
//------------------初期化------------------//
	for(k=0; k<=para.data-1; k++){
		tg[k].re = 0.0;
		tg[k].im = 0.0;
		td[k].re = 0.0;
		td[k].im = 0.0;
	}

//---------電流揺らぎのデータの代入---------//
	for(k=0; k<=para.numn-1; k++){
		tg[k].re = cur[para.jtn0-para.numn+k]._[0];
		td[k].re = cur[para.jtn0-para.numn+k]._[1];			
	}
}

/* ************************************ */
/* ----- タイムウィンドウを掛ける ----- */
void TimeWindow(struct ComplexNumbers t[], struct SetPara para){
	int i;
	/* 0～T/10 までのウィンド */
	for(i=0; i<=(int)(para.data/10); i++){
		t[i].re = t[i].re * (1-cos(para.pai*10*i/para.data))/2;
	}

	/* 9T/10～T までのウィンド */
	for(i=(int)(para.data*9/10); i<=para.data-1; i++){
//		t[i].re = t[i].re * (1+cos(para.pai*10*((i-(int)(para.data*9/10))/para.data)))/2;
		t[i].re = t[i].re * (cos(para.pai*10*((i-(int)(para.data*9/10))/para.data)))/2;
	}

}

/* ******************************** */
/* ----- スペクトルを計算する ----- */
void CalculateSpectrum_FFT(struct Spectrum si[], struct ComplexNumbers tg[], 
						   struct ComplexNumbers td[], struct SetPara para){
	int jf;
	double st;

	if(para.twin==0){st = 1.0;}
	else if(para.twin==1){st = 0.875;}/*タイムウィンドウを掛けた際の強度調整 */

	for(jf=0; jf<=para.data/2-1; jf++){
		si[jf].g    = 2*fabs(tg[jf].re*tg[jf].re + tg[jf].im*tg[jf].im)/st;
		si[jf].d    = 2*fabs(td[jf].re*td[jf].re + td[jf].im*td[jf].im)/st;
		si[jf].gd.re= 2*fabs(tg[jf].re*td[jf].re + tg[jf].im*td[jf].im)/st;
		si[jf].gd.im= 2*fabs(-tg[jf].im*td[jf].re + tg[jf].re*td[jf].im)/st;

		td[jf].re = td[jf].re * sqrt(2*para.dt/para.numn/st/para.epp);
		td[jf].im = td[jf].im * sqrt(2*para.dt/para.numn/st/para.epp);
		tg[jf].re = tg[jf].re * sqrt(2*para.dt/para.numn/st/para.epp);
		tg[jf].im = tg[jf].im * sqrt(2*para.dt/para.numn/st/para.epp);
	}
	
}

/* **************************************************** */
/* ----- スペクトルウィンドによる周波数平滑を行う ----- */
void SpectrumWindow_FFT(struct Spectrum si[], 
					    struct ComplexNumbers tg[],struct ComplexNumbers td[], 
						double fr[], struct SetPara para){

	int jf, fmax;	
	double *x, *y, *z, *w;

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	x = new double [para.data];
	y = new double [para.data];
	z = new double [para.data];
	w = new double [para.data];

	/* 最大周波数の設定 */
//	i = 1/para.dt/2.56;
//	fmax = i / (1/para.dt*para.data);
	fmax = para.data/2-1;
//	fmax = para.fmax_y;

	/* ----- 周波数ウィンドによる周波数平滑 ----- */
	for(jf=0; jf<=fmax; jf++){
//		if(fr[jf]>para.fmax_y){break;} //最大周波数で計算をストップ
		x[jf] = si[jf].g;
		y[jf] = si[jf].d;
		z[jf] = si[jf].gd.re;
		w[jf] = si[jf].gd.im;
	}

//	if(para.saveg!=0){HanningWindow_f(x, fmax);}
//	if(para.saved!=0){HanningWindow_f(y, fmax);}
//	if(para.savegd!=0){
//		HanningWindow_f(z, fmax);
//		HanningWindow_f(w, fmax);

	if(para.saveg!=0){TriangleWindow(x, fmax, para.saveg, para.data);}
	if(para.saved!=0){TriangleWindow(y, fmax, para.saved, para.data);}
	if(para.savegd!=0){
		TriangleWindow(z, fmax, para.savegd, para.data);
		TriangleWindow(w, fmax, para.savegd, para.data);
	}
	for(jf=0; jf<=fmax; jf++){
		si[jf].g= x[jf] *para.dt/para.numn/para.epp;
		si[jf].d= y[jf] *para.dt/para.numn/para.epp;
		si[jf].gd.re= z[jf] *para.dt/para.numn/para.epp;
		si[jf].gd.im= w[jf] *para.dt/para.numn/para.epp;
	}

	/* ----- メモリの開放 */
	delete [] x;
	delete [] y;
	delete [] z;
	delete [] w;

}

/* ****************************************************************** */
/* ----- [Sig ∝ f^2]　[Sid = Const]　[Sigid ∝ f] の計算 */

void SpectrumAppro_FFT(struct Spectrum si[], double fr[], struct SetPara para){

	int jf, fmax;
	int l, fmin_index, fmax_index;
	double a, rs, b, rsb;
	double freq;

	/* 最大周波数の設定 */
	fmax = para.data/2-1;


	/* --- Sid : 一定値の計算 */
	if(para.appr_sd==1){
		/* 近似を行うための周波数帯域の周波数インデックスを設定 */
		for(jf=0; jf<=fmax; jf++){
			freq = jf/para.dt/para.data;
			if(freq >= para.fmin_sd){fmin_index=jf;	break;}
		}
		l=0;
		for(jf=fmin_index; jf<=fmax; jf++){
			freq = jf/para.dt/para.data;
			if(freq >= para.fmax_sd){fmax_index=jf-1;	break;}
			l = l+1;
		}
		/* 平滑化の計算 */
		a = 0;	rs = 0;
		for(jf=fmin_index; jf<=fmax_index; jf++){
			Sigma(si[jf].d, &a, &rs);
		}
		a = a / ((float)l);
		for(jf=0; jf<=fmax; jf++){
			freq = jf/para.dt/para.data;
			si[jf].d = a;
			if(freq > para.fmax_sd2){break;}
		}
	}


	/* --- Sig ∝ f^2 の計算 */
	if(para.appr_sg==1){
		/* 近似を行うための周波数帯域の周波数インデックスを設定 */
		for(jf=0; jf<=fmax; jf++){
			freq = jf/para.dt/para.data;
			if(freq >= para.fmin_sg){fmin_index=jf;	break;}
		}
		l=0;
		for(jf=fmin_index; jf<=fmax; jf++){
			freq = jf/para.dt/para.data;
			if(freq >= para.fmax_sg){fmax_index=jf-1;	break;}
			l = l+1;
		}
		/* 平滑化の計算 */
		a = 0;	rs = 0;
		for(jf=fmin_index; jf<=fmax_index; jf++){
			freq = jf/para.dt/para.data;
			Sigma(si[jf].g/pow(freq, 2.0), &a, &rs);
		}
		a = a / ((float)l);
		for(jf=0; jf<=fmax; jf++){
			freq = jf/para.dt/para.data;
			si[jf].g = a * pow(freq, 2.0);
			if(freq > para.fmax_sg2){break;}
		}

	}


	/* --- Sigd ∝ f の計算 */
	if(para.appr_sgd==1){
		/* 近似を行うための周波数帯域の周波数インデックスを設定 */
		for(jf=0; jf<=fmax; jf++){
			freq = jf/para.dt/para.data;
			if(freq >= para.fmin_sgd){fmin_index=jf;	break;}
		}
		l=0;
		for(jf=fmin_index; jf<=fmax; jf++){
			freq = jf/para.dt/para.data;
			if(freq >= para.fmax_sgd){fmax_index=jf-1;	break;}
			l = l+1;
		}
		/* 平滑化の計算 */
		a = 0;	rs = 0;
		b = 0;	rsb = 0;
		for(jf=fmin_index; jf<=fmax_index; jf++){
			freq = jf/para.dt/para.data;
			Sigma(si[jf].gd.im/freq, &a, &rs);
			Sigma(si[jf].gd.re/freq, &b, &rsb);
		}
		a = a / ((float)l);
		b = b / ((float)l);
		for(jf=0; jf<=fmax; jf++){
			freq = jf/para.dt/para.data;
			si[jf].gd.im = a * freq;
			si[jf].gd.re = b * freq;
			if(freq > para.fmax_sgd2){break;}
		}

	}
}

/* ******************************************** */
/* ----- スペクトルから相関関数を計算する ----- */
void CalculateCorrelation_FFT(struct Spectrum si[], struct Correlation corr[],
							  struct SetPara para){
	int i;
	struct ComplexNumbers *temp1, *temp2, *temp3, *temp4;

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	temp1 = new ComplexNumbers [para.data];
	temp2 = new ComplexNumbers [para.data];
	temp3 = new ComplexNumbers [para.data];
	temp4 = new ComplexNumbers [para.data];

	/* スペクトルを一時変数配列に格納 */
	for(i=0; i<=para.data/2-1; i++){
		temp1[i].re= si[i].g; /*自己相関関数 sig */
		temp1[i].im= 0.0;
		temp2[i].re= si[i].d; /*自己相関関数 sig */
		temp2[i].im= 0.0;
		temp3[i].re= si[i].gd.re; /*相互相関関数 sigd */
		temp3[i].im= 0.0;
		temp4[i].re= si[i].gd.im; /*相互相関関数 sigd */
		temp4[i].im= 0.0;
	}
	for(i=para.data/2; i<=para.data-1; i++){
		temp1[i].re= 0.0;	temp1[i].im= 0.0;
		temp2[i].re= 0.0;	temp2[i].im= 0.0;
		temp3[i].re= 0.0;	temp3[i].im= 0.0;
		temp4[i].re= 0.0;	temp4[i].im= 0.0;
	}

	/* スペクトルをFFTにかける */
	FFTCalculate3(temp1, para);
	FFTCalculate3(temp2, para);
	FFTCalculate3(temp3, para);
	FFTCalculate3(temp4, para);

	/* スペクトルのFFTから相関関数を計算する */
	for(i=0; i<=para.data/4-1; i++) {
		corr[i]._[0] = 2*temp1[i].re/para.data/para.dt;
		corr[i]._[1] = 2*temp2[i].re/para.data/para.dt;
		corr[i]._[2] = 2*(temp3[i].re-temp4[i].im)/para.data/para.dt;
		corr[i]._[3] = 2*(temp3[i].re+temp4[i].im)/para.data/para.dt;
	}




	/* ----- メモリの開放 */
	delete [] temp1;
	delete [] temp2;
	delete [] temp3;
	delete [] temp4;

}


/* ********************************************** */
/* ----- FFT計算（"fftsg.cpp"を用いた関数） ----- */
void FFTCalculate2(struct ComplexNumbers t[], struct SetPara para){

	void cdft(int, int, double *, int *, double *);
	int i, *ip;
	double *a, *w;

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	ip = new int [para.data];
	a = new double [2*para.data];
	w = new double [(int)(para.data/2)];

	ip[0] = 0;
	for(i=0; i<=para.data-1; i++) {
		a[2*i] = t[i].re;
		a[2*i+1] = t[i].im;
	}

	cdft(para.data*2, -1, a, ip, w);

	for(i=0; i<=para.data/2-1; i++) {//FFT後の出力数はこれでよいのか
		t[i].re = a[2*i];
		t[i].im = a[2*i+1];
	}

	/* ----- メモリの開放 */
	delete [] ip;
	delete [] a;
	delete [] w;

}


/* ********************************************** */
/* ----- FFT計算（"fftsg.cpp"を用いた関数） ----- */
void FFTCalculate3(struct ComplexNumbers t[], struct SetPara para){

	void cdft(int, int, double *, int *, double *);
	int i, *ip;
	double *a, *w;

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	ip = new int [para.data/2];
	a = new double [para.data];
	w = new double [para.data * 5 / 8];

	ip[0] = 0;
	for(i=0; i<=para.data/2-1; i++) {
		a[2*i] = t[i].re;
		a[2*i+1] = t[i].im;
	}

	cdft(para.data, -1, a, ip, w);

	for(i=0; i<=para.data/2-1; i++) {
		t[i].re = a[2*i];
		t[i].im = a[2*i+1];
	}
	
	/* ----- メモリの開放 */
	delete [] ip;
	delete [] a;
	delete [] w;

}



/* ************************************************ */
/* ----- FILEへスペクトルと相関関数を書き込む ----- */
void WriteSpectrum_FFT(struct Correlation corr[], struct Spectrum si[], 
					   struct ComplexNumbers tg[],struct ComplexNumbers td[], 
					   double fr[], struct SetPara para) {
	int k;
    FILE *fp_spe, *fp_cor, *fp_current, *fp;

	/* --- 平均電流値の書き込み */
	if(para.current_out ==1){
		fp_current = OpenFile(para.current, 'a');
		fprintf(fp_current, "ΔI_G	ΔI_D\n");
		fprintf(fp_current, "%1.4e	%1.4e\n",sqrt(corr[0]._[0]), sqrt(corr[0]._[1]));
		fclose(fp_current);
		printf("'%s' write\n", para.current);
	}

	/* --- 相関関数の書き込み */
	if(para.cor_out ==1){
		fp_cor = OpenFile(para.cor, 'w');
		fprintf(fp_cor, "t[sec]	Cor_d(t)	Cor_g(t)	CrossCor(t)\n");
		for(k=0; k<=para.data/4-1; k++){
			if(k > para.fnum){break;}
			fprintf(fp_cor, "%1.4e	", 2*k*para.dt);
			fprintf(fp_cor, "%1.4e	%1.4e	%1.4e", corr[k]._[1], corr[k]._[0], corr[k]._[2]);
			fprintf(fp_cor, "\n");
		}
		fclose(fp_cor);
		printf("'%s' write\n", para.cor);
	}


	/* --- フーリエ変換成分の書き込み */
	if(para.fouricomp_out ==1){
		fp = OpenFile(para.fouricomp, 'w');
		fprintf(fp, "f[Hz]	Xd(f).re	Xd(f).im	Xg(f).re	Xg(f).im\n");
		for(k=0; k<=para.data/4-1; k++){
			if(k > para.fnum){break;}

			fprintf(fp, "%1.4e	", k/para.dt/para.data);
			fprintf(fp, "%1.4e	" , td[k].re);
			fprintf(fp, "%1.4e	" , td[k].im);
			fprintf(fp, "%1.4e	" , tg[k].re);
			fprintf(fp, "%1.4e" , tg[k].im);
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf("'%s' write\n", para.fouricomp);
	}


	/* --- スペクトルの書き込み */
	fp_spe = OpenFile(para.spe, 'w');
	fprintf(fp_spe, "f[Hz]	Si_D(f)	Si_G(f)	Si_G-D(f).re	Si_G-D(f).im\n");
	for(k=0; k<=para.data/4-1; k++){
		if(k > para.fnum){break;}
		fprintf(fp_spe, "%1.4e	", k/para.dt/para.data);
		fprintf(fp_spe, "%1.4e	" , si[k].d);
		fprintf(fp_spe, "%1.4e	" , si[k].g);
		fprintf(fp_spe, "%1.4e	" , si[k].gd.re);
		fprintf(fp_spe, "%1.4e"   , si[k].gd.im);

		fprintf(fp_spe, "\n");
	}
	fclose(fp_spe);
	printf("'%s' write\n", para.spe);


}
